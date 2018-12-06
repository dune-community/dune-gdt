// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2016 - 2017)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
#define DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH

#include <algorithm>
#include <cmath>
#include <list>
#include <memory>

#include <boost/align/aligned_allocator.hpp>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/container/pattern.hh>

#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/operators/fv/reconstruction/internal.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#if HAVE_CLP
#include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

namespace Dune {
namespace GDT {


template <class KeyVectorType, class ValueVectorType>
class EntropyLocalCache
{
public:
  using MapType = typename std::map<KeyVectorType, ValueVectorType, XT::Common::VectorLess>;
  using KAbs = XT::Common::VectorAbstraction<KeyVectorType>;
  using VAbs = XT::Common::VectorAbstraction<ValueVectorType>;
  using IteratorType = typename MapType::iterator;
  using ConstIteratorType = typename MapType::const_iterator;
  using RangeFieldType = typename KAbs::ScalarType;

  EntropyLocalCache(const size_t capacity)
    : capacity_(capacity)
  {}

  void insert(const KeyVectorType& u, const ValueVectorType& alpha)
  {
    cache_.insert(std::make_pair(u, alpha));
    keys_.push_back(u);
    if (cache_.size() > capacity_) {
      cache_.erase(keys_.front());
      keys_.pop_front();
    }
  }

  void keep(const KeyVectorType& u)
  {
    keys_.remove(u);
    keys_.push_back(u);
  }

  ConstIteratorType find_closest(const KeyVectorType& u) const
  {
    ConstIteratorType ret = cache_.begin();
    if (ret == end())
      return ret;
    auto diff = u - ret->first;
    // use infinity_norm as distance
    RangeFieldType distance = std::abs(diff[0]);
    for (size_t ii = 1; ii < diff.size(); ++ii)
      distance = std::max(distance, diff[ii]);
    auto it = ret;
    while (++it != end()) {
      diff = u - it->first;
      RangeFieldType new_distance = std::abs(diff[0]);
      for (size_t ii = 1; ii < diff.size(); ++ii)
        new_distance = std::max(new_distance, diff[ii]);
      if (new_distance < distance) {
        distance = new_distance;
        ret = it;
      }
    }
    return ret;
  }

  IteratorType begin()
  {
    return cache_.begin();
  }

  ConstIteratorType begin() const
  {
    return cache_.begin();
  }

  IteratorType end()
  {
    return cache_.end();
  }

  ConstIteratorType end() const
  {
    return cache_.end();
  }

  void set_capacity(const size_t new_capacity)
  {
    capacity_ = new_capacity;
  }

private:
  std::size_t capacity_;
  MapType cache_;
  std::list<KeyVectorType> keys_;
};

template <class BasisfunctionImp, class GridLayerImp, class U>
class EntropyBasedLocalFlux;

/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 */
template <class BasisfunctionImp, class GridLayerImp, class U>
class EntropyBasedLocalFlux
  : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                           typename BasisfunctionImp::DomainFieldType,
                                                           BasisfunctionImp::dimFlux,
                                                           U,
                                                           0,
                                                           typename BasisfunctionImp::RangeFieldType,
                                                           BasisfunctionImp::dimRange,
                                                           BasisfunctionImp::dimFlux>
{
  using BaseType =
      typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                               typename BasisfunctionImp::DomainFieldType,
                                                               BasisfunctionImp::dimFlux,
                                                               U,
                                                               0,
                                                               typename BasisfunctionImp::RangeFieldType,
                                                               BasisfunctionImp::dimRange,
                                                               BasisfunctionImp::dimFlux>;
  using ThisType = EntropyBasedLocalFlux;

public:
  using BasisfunctionType = BasisfunctionImp;
  using GridLayerType = GridLayerImp;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using typename BaseType::StateType;
  // make matrices a little larger to align to 64 byte boundary
  static constexpr size_t matrix_num_cols = dimRange % 8 ? dimRange : dimRange + (8 - dimRange % 8);
  using MatrixType = XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using VectorType = XT::Common::FieldVector<RangeFieldType, dimRange>;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using BasisValuesMatrixType = XT::LA::CommonDenseMatrix<RangeFieldType>;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, dimDomain>;
  using AlphaReturnType = std::pair<VectorType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, VectorType>;
  using AlphaStorageType = std::map<DomainType, StateRangeType, XT::Common::VectorFloatLess>;
  static const size_t cache_size = 4 * dimDomain + 2;

  // get permutation instead of sorting directly to be able to sort two vectors the same way
  // see
  // https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  template <typename T, typename Compare>
  std::vector<std::size_t> get_sort_permutation(const std::vector<T>& vec, const Compare& compare)
  {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
    return p;
  }

  template <typename T>
  void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p)
  {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
      if (done[i]) {
        continue;
      }
      done[i] = true;
      std::size_t prev_j = i;
      std::size_t j = p[i];
      while (i != j) {
        std::swap(vec[prev_j], vec[j]);
        done[j] = true;
        prev_j = j;
        j = p[j];
      }
    }
  }

  // Joins duplicate quadpoints, vectors have to be sorted!
  void join_duplicate_quadpoints(std::vector<DomainType>& quad_points, std::vector<RangeFieldType>& quad_weights)
  {
    // Index of first quad_point of several quad_points with the same position
    size_t curr_index = 0;
    std::vector<size_t> indices_to_remove;
    for (size_t ll = 1; ll < quad_weights.size(); ++ll) {
      if (XT::Common::FloatCmp::eq(quad_points[curr_index], quad_points[ll])) {
        quad_weights[curr_index] += quad_weights[ll];
        indices_to_remove.push_back(ll);
      } else {
        curr_index = ll;
      }
    } // ll
    assert(indices_to_remove.size() < std::numeric_limits<int>::max());
    // remove duplicate points, from back to front to avoid invalidating indices
    for (int ll = static_cast<int>(indices_to_remove.size()) - 1; ll >= 0; --ll) {
      quad_points.erase(quad_points.begin() + indices_to_remove[ll]);
      quad_weights.erase(quad_weights.begin() + indices_to_remove[ll]);
    }
  }

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , quad_points_(basis_functions_.quadratures().merged().size())
    , quad_weights_(quad_points_.size())
    , M_(quad_points_.size(), matrix_num_cols, 0., 0)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , T_minus_one_(std::make_unique<MatrixType>())
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , alpha_storage_(index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {
    XT::LA::eye_matrix(*T_minus_one_);
    size_t ll = 0;
    for (const auto& quad_point : basis_functions_.quadratures().merged()) {
      quad_points_[ll] = quad_point.position();
      quad_weights_[ll] = quad_point.weight();
      ++ll;
    }
    // Join duplicate quad_points. For that purpose, first sort the vectors
    const auto permutation = get_sort_permutation(quad_points_, XT::Common::VectorFloatLess{});
    apply_permutation_in_place(quad_points_, permutation);
    apply_permutation_in_place(quad_weights_, permutation);
    // Now join duplicate quad_points by removing all quad_points with the same position except one and adding the
    // weights of the removed points to the remaining point
    join_duplicate_quadpoints(quad_points_, quad_weights_);
    assert(quad_points_.size() == quad_weights_.size());
    // evaluate basis functions and store in matrix
    M_.resize(quad_points_.size(), matrix_num_cols);
    for (ll = 0; ll < quad_points_.size(); ++ll) {
      const auto val = basis_functions_.evaluate(quad_points_[ll]);
      for (size_t ii = 0; ii < dimRange; ++ii)
        M_.set_entry(ll, ii, val[ii]);
    }
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using LocalfunctionType::dimRange;
    using typename LocalfunctionType::ColPartialURangeType;
    using typename LocalfunctionType::ColRangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const std::vector<DomainType>& quad_points,
                  const std::vector<RangeFieldType>& quad_weights,
                  const BasisValuesMatrixType& M,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const MatrixType& T_minus_one,
                  LocalCacheType& cache,
                  AlphaStorageType& alpha_storage,
                  std::mutex& mutex
#if HAVE_CLP
                  ,
                  XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp)
#else
                  )
#endif
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , quad_points_(quad_points)
      , quad_weights_(quad_weights)
      , M_(M)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , T_minus_one_(T_minus_one)
      , cache_(cache)
      , alpha_storage_(alpha_storage)
      , mutex_(mutex)
#if HAVE_CLP
      , realizability_helper_(basis_functions_, quad_points_, lp)
#else
      , realizability_helper_(basis_functions_, quad_points_)
#endif
    {}

    template <class BasisFuncImp = BasisfunctionType, bool quadrature_contains_vertices = true, bool anything = true>
    struct RealizabilityHelper;

#if HAVE_CLP
    template <class BasisFuncImp, bool quadrature_contains_vertices, bool anything>
    struct RealizabilityHelper
    {
      static_assert(std::is_same<BasisFuncImp, BasisfunctionType>::value, "BasisFuncImp has to be BasisfunctionType!");

      RealizabilityHelper(const BasisfunctionType& basis_functions,
                          const std::vector<DomainType>& quad_points,
                          XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp)
        : basis_functions_(basis_functions)
        , quad_points_(quad_points)
        , lp_(lp)
      {}

      // The ClpSimplex structure seems to get corrupted sometimes (maybe some problems with infs/NaNs?), so we
      // reinitialize it if the stopping conditions is always false
      void setup_linear_program(const bool reinitialize) const
      {
        if (!*lp_ || reinitialize) {
          // We start with creating a model with dimRange rows and num_quad_points columns */
          constexpr int num_rows = static_cast<int>(dimRange);
          assert(quad_points_.size() < std::numeric_limits<int>::max());
          int num_cols = static_cast<int>(quad_points_.size()); /* variables are x_1, ..., x_{num_quad_points} */
          *lp_ = std::make_unique<ClpSimplex>(false);
          auto& lp = **lp_;
          // set number of rows
          lp.resize(num_rows, 0);

          // Clp wants the row indices that are non-zero in each column. We have a dense matrix, so provide all indices
          // 0..num_rows
          std::array<int, num_rows> row_indices;
          for (int ii = 0; ii < num_rows; ++ii)
            row_indices[static_cast<size_t>(ii)] = ii;

          // set columns for quadrature points
          for (int ii = 0; ii < num_cols; ++ii) {
            const auto v_i = basis_functions_.evaluate(quad_points_[static_cast<size_t>(ii)]);
            // First argument: number of elements in column
            // Second/Third argument: indices/values of column entries
            // Fourth/Fifth argument: lower/upper column bound, i.e. lower/upper bound for x_i. As all x_i should be
            // positive, set to 0/inf, which is the default.
            // Sixth argument: Prefactor in objective for x_i, this is 0 for all x_i, which is also the default;
            lp.addColumn(num_rows, row_indices.data(), &(v_i[0]));
          }

          // silence lp
          lp.setLogLevel(0);
        } // if (!lp_)
      }

      bool is_realizable(const StateRangeType& u, const bool reinitialize) const
      {
        const auto density = basis_functions_.density(u);
        if (!(density > 0.) || std::isinf(density))
          return false;
        const auto u_prime = u / density;
        setup_linear_program(reinitialize);
        auto& lp = **lp_;
        constexpr int num_rows = static_cast<int>(dimRange);
        // set rhs (equality constraints, so set both bounds equal
        for (int ii = 0; ii < num_rows; ++ii) {
          size_t uii = static_cast<size_t>(ii);
          lp.setRowLower(ii, u_prime[uii]);
          lp.setRowUpper(ii, u_prime[uii]);
        }
        // set maximal wall time. If this is not set, in rare cases the primal method never returns
        lp.setMaximumWallSeconds(60);
        // Now check solvability
        lp.primal();
        return lp.primalFeasible();
      }

    private:
      const BasisfunctionType& basis_functions_;
      const std::vector<DomainType>& quad_points_;
      XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp_;
    }; // struct RealizabilityHelper<...>
#else // HAVE_CLP
    template <class BasisFuncImp, bool quadrature_contains_vertices, bool anything>
    struct RealizabilityHelper
    {
      RealizabilityHelper(const BasisfunctionType& /*basis_functions*/, const std::vector<DomainType>& /*quad_points*/)
      {
        DUNE_THROW(Dune::NotImplemented, "You are missing Clp!");
      }

      bool is_realizable(const StateRangeType& /*u*/, const bool /*reinitialize*/) const
      {
        DUNE_THROW(Dune::NotImplemented, "You are missing Clp!");
        return false;
      }
    }; // struct RealizabilityHelper<...>
#endif // HAVE_CLP

    // specialization for hatfunctions
    template <size_t dimRange_or_refinements, bool anything>
    struct RealizabilityHelper<
        HatFunctionMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange_or_refinements, 1, dimDomain>,
        true,
        anything>
    {
      RealizabilityHelper(const BasisfunctionType& /*basis_functions*/,
                          const std::vector<DomainType>& /*quad_points*/
#if HAVE_CLP
                          ,
                          XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& /*lp*/)
#else
      )
#endif
      {}

      static bool is_realizable(const StateRangeType& u, const bool /*reinitialize*/)
      {
        for (const auto& u_i : u)
          if (!(u_i > 0.) || std::isinf(u_i))
            return false;
        return true;
      }
    }; // struct RealizabilityHelper<Hatfunctions, ...>

    void keep(const StateRangeType& u)
    {
      cache_.keep(u);
    }

    using LocalfunctionType::entity;

    // temporary vectors to store inner products and exponentials
    std::vector<RangeFieldType>& working_storage() const
    {
      thread_local std::vector<RangeFieldType> work_vec;
      work_vec.resize(quad_points_.size());
      return work_vec;
    }

    void calculate_scalar_products(const VectorType& beta_in,
                                   const BasisValuesMatrixType& M,
                                   std::vector<RangeFieldType>& scalar_products) const
    {
#if HAVE_MKL || HAVE_CBLAS
      XT::Common::Blas::dgemv(XT::Common::Blas::row_major(),
                              XT::Common::Blas::no_trans(),
                              static_cast<int>(quad_points_.size()),
                              dimRange,
                              1.,
                              M.data(),
                              matrix_num_cols,
                              &(beta_in[0]),
                              1,
                              0.,
                              scalar_products.data(),
                              1);
#else
      const size_t num_quad_points = quad_points_.size();
      std::fill(scalar_products.begin(), scalar_products.end(), 0.);
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto* basis_ll = M.get_ptr(ll);
        scalar_products[ll] = std::inner_product(beta_in.begin(), beta_in.end(), basis_ll, 0.);
      }
#endif
    }

    void apply_exponential(std::vector<RangeFieldType>& values) const
    {
      assert(values.size() < std::numeric_limits<int>::max());
      XT::Common::Mkl::exp(static_cast<int>(values.size()), values.data(), values.data());
    }

    // calculate ret = \int (exp(beta_in * m))
    RangeFieldType calculate_scalar_integral(const VectorType& beta_in, const BasisValuesMatrixType& M) const
    {
      auto& work_vec = working_storage();
      calculate_scalar_products(beta_in, M, work_vec);
      apply_exponential(work_vec);
      return std::inner_product(quad_weights_.begin(), quad_weights_.end(), work_vec.begin(), RangeFieldType(0.));
    }

    // calculate ret = \int (m1 exp(beta_in * m2))
    void calculate_vector_integral(const VectorType& beta_in,
                                   const BasisValuesMatrixType& M1,
                                   const BasisValuesMatrixType& M2,
                                   VectorType& ret,
                                   bool same_beta = false,
                                   bool only_first_component = false) const
    {
      auto& work_vec = working_storage();
      if (!same_beta) {
        calculate_scalar_products(beta_in, M2, work_vec);
        apply_exponential(work_vec);
      }
      std::fill(ret.begin(), ret.end(), 0.);
      const size_t num_quad_points = quad_weights_.size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto factor_ll = work_vec[ll] * quad_weights_[ll];
        const auto* basis_ll = M1.get_ptr(ll);
        for (size_t ii = 0; ii < (only_first_component ? 1 : dimRange); ++ii)
          ret[ii] += basis_ll[ii] * factor_ll;
      } // ll
    }

    void copy_transposed(const MatrixType& T_k, MatrixType& T_k_trans) const
    {
      for (size_t ii = 0; ii < dimRange; ++ii)
        for (size_t jj = 0; jj <= ii; ++jj)
          T_k_trans[jj][ii] = T_k[ii][jj];
    }

    // For each basis evaluation b, calculates T_k^{-1} b. As the basis evaluations are the rows of M, we want to
    // calculate (T_k^{-1} M^T)^T = M T_k^{-T}
    void apply_inverse_matrix(const MatrixType& T_k, BasisValuesMatrixType& M) const
    {
#if HAVE_MKL || HAVE_CBLAS
      // Calculate the transpose here first as this is much faster than passing the matrix to dtrsm and using CblasTrans
      thread_local auto T_k_trans = std::make_unique<MatrixType>(0.);
      copy_transposed(T_k, *T_k_trans);
      assert(quad_points_.size() < std::numeric_limits<int>::max());
      XT::Common::Blas::dtrsm(XT::Common::Blas::row_major(),
                              XT::Common::Blas::right(),
                              XT::Common::Blas::upper(),
                              XT::Common::Blas::no_trans(),
                              XT::Common::Blas::non_unit(),
                              static_cast<int>(quad_points_.size()),
                              dimRange,
                              1.,
                              &((*T_k_trans)[0][0]),
                              dimRange,
                              M.data(),
                              matrix_num_cols);
#else
      assert(quad_points_.size() == M.rows());
      VectorType tmp_vec, tmp_vec2;
      for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
        std::copy_n(M.get_ptr(ll), dimRange, tmp_vec.begin());
        XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
        std::copy_n(tmp_vec2.begin(), dimRange, M.get_ptr(ll));
      }
#endif
    }

    void store_alpha(const DomainType& x_local, const StateRangeType& alpha)
    {
      alpha_storage_[x_local] = alpha;
    }

    StateRangeType get_stored_alpha(const DomainType& x_local) const
    {
      return alpha_storage_.at(x_local);
    }

    template <class GridLayerType>
    void center_results_to_intersections(const GridLayerType& grid_layer)
    {
      const auto center = entity().geometry().local(entity().geometry().center());
      const auto center_alpha = get_stored_alpha(center);
      for (const auto& intersection : Dune::intersections(grid_layer, entity()))
        store_alpha(entity().geometry().local(intersection.geometry().center()), center_alpha);
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& x_local,
                                               const StateRangeType& u,
                                               const XT::Common::Parameter& param,
                                               const bool regularize) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      auto ret = std::make_unique<AlphaReturnType>();
      mutex_.lock();
      if (boundary)
        cache_.set_capacity(cache_size + dimDomain);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      if (!(density > 0.) || std::isinf(density)) {
        mutex_.unlock();
        DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
      }
      VectorType u_prime = u / density;
      auto alpha_iso_dyn = basis_functions_.alpha_iso();
      VectorType alpha_iso;
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_iso[ii] = alpha_iso_dyn[ii];
      VectorType v, u_eps_diff, beta_in;
      RangeFieldType first_error_cond, second_error_cond, tau_prime;

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u_prime, 1e-14, 1e-14)) {
        const auto alpha_prime = cache_iterator->second;
        ret->first = alpha_prime + alpha_iso * std::log(density);
        ret->second = 0.;
        alpha_storage_[x_local] = ret->first;
        mutex_.unlock();
        return ret;
      } else {
        auto u_iso = basis_functions_.u_iso();
        const RangeFieldType dim_factor = is_full_moment_basis<BasisfunctionType>::value ? 1. : std::sqrt(dimDomain);
        tau_prime = std::min(tau_ / ((1 + dim_factor * u_prime.two_norm()) * density + dim_factor * tau_), tau_);

        // define further variables
        VectorType g_k, beta_out;
        beta_in = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;
        thread_local auto T_k = XT::Common::make_unique<MatrixType>();

        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence) {
          // regularize u
          v = u_prime;
          if (r > 0) {
            beta_in = alpha_iso;
            DynamicRangeType r_times_u_iso = u_iso;
            r_times_u_iso *= r;
            v *= 1 - r;
            v += r_times_u_iso;
          }
          *T_k = T_minus_one_;
          // calculate T_k u
          VectorType v_k = v;
          // calculate values of basis p = S_k m
          thread_local BasisValuesMatrixType P_k(M_.backend(), false, 0., 0);
          std::copy_n(M_.data(), M_.rows() * M_.cols(), P_k.data());
          // calculate f_0
          RangeFieldType f_k = calculate_scalar_integral(beta_in, P_k);
          f_k -= beta_in * v_k;

          thread_local auto H = XT::Common::make_unique<MatrixType>(0.);

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used or cholesky decomposition fails
            if (kk > k_0_ && r < r_max)
              break;
            try {
              change_basis(beta_in, v_k, P_k, *T_k, g_k, beta_out, *H);
            } catch (const Dune::MathError&) {
              if (r < r_max)
                break;
              mutex_.unlock();
              const std::string err_msg =
                  "Failed to converge for " + XT::Common::to_string(u) + " with density "
                  + XT::Common::to_string(density) + " and multiplier " + XT::Common::to_string(beta_in)
                  + " at position " + XT::Common::to_string(entity().geometry().center())
                  + " due to errors in change_basis! Last u_eps_diff = " + XT::Common::to_string(u_eps_diff)
                  + ", first_error_cond = " + XT::Common::to_string(first_error_cond) + ", second_error_cond = "
                  + XT::Common::to_string(second_error_cond) + ", tau_prime = " + XT::Common::to_string(tau_prime);
              DUNE_THROW(MathError, err_msg);
            }
            // calculate descent direction d_k;
            VectorType d_k = g_k;
            d_k *= -1;
            // Calculate stopping criteria (in original basis). Variables with _k are in current basis, without k in
            // original basis.
            VectorType alpha_tilde;
            XT::LA::solve_lower_triangular_transposed(*T_k, alpha_tilde, beta_out);
            VectorType u_alpha_tilde;
            calculate_vector_integral(alpha_tilde, M_, M_, u_alpha_tilde);
            VectorType g_alpha_tilde = u_alpha_tilde - v;
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            if (!(density_tilde > 0.) || std::isinf(density_tilde))
              break;
            const auto alpha_prime = alpha_tilde - alpha_iso * std::log(density_tilde);
            VectorType u_alpha_prime;
            calculate_vector_integral(alpha_prime, M_, M_, u_alpha_prime);
            u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
            VectorType d_alpha_tilde;
            XT::LA::solve_lower_triangular_transposed(*T_k, d_alpha_tilde, d_k);
            first_error_cond = g_alpha_tilde.two_norm();
            second_error_cond = std::exp(d_alpha_tilde.one_norm() + std::abs(std::log(density_tilde)));
            if (first_error_cond < tau_prime && 1 - epsilon_gamma_ < second_error_cond
                && realizability_helper_.is_realizable(u_eps_diff, kk == static_cast<size_t>(0.8 * k_0_))) {
              ret->first = alpha_prime + alpha_iso * std::log(density);
              ret->second = r;
              cache_.insert(v, alpha_prime);
              alpha_storage_[x_local] = ret->first;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              beta_in = beta_out;
              // backtracking line search
              // while (pure_newton >= 2 || zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm() * 100.) {
              while (pure_newton >= 2 || zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm()) {
                VectorType beta_new = d_k;
                beta_new *= zeta_k;
                beta_new += beta_out;
                RangeFieldType f = calculate_scalar_integral(beta_new, P_k);
                f -= beta_new * v_k;
                if (pure_newton >= 2 || XT::Common::FloatCmp::le(f, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  beta_in = beta_new;
                  f_k = f;
                  pure_newton = 0;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              if (zeta_k <= epsilon_ * beta_out.two_norm() / d_k.two_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)
        mutex_.unlock();
        const std::string err_msg =
            "Failed to converge for " + XT::Common::to_string(u) + " with density " + XT::Common::to_string(density)
            + " and multiplier " + XT::Common::to_string(beta_in) + " at position "
            + XT::Common::to_string(entity().geometry().center()) + " due to too many iterations! Last u_eps_diff = "
            + XT::Common::to_string(u_eps_diff) + ", first_error_cond = " + XT::Common::to_string(first_error_cond)
            + ", second_error_cond = " + XT::Common::to_string(second_error_cond)
            + ", tau_prime = " + XT::Common::to_string(tau_prime);
        DUNE_THROW(MathError, err_msg);
      } // else ( value has not been calculated before )

    outside_all_loops:
      mutex_.unlock();
      return ret;
    }

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      ColRangeType col_ret;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        evaluate_col(dd, x_local, u, col_ret, param);
        for (size_t ii = 0; ii < dimRange; ++ii)
          helper<dimDomain>::get_ref(ret, ii, dd) = col_ret[ii];
      } // dd
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      std::fill(ret.begin(), ret.end(), 0.);
      const auto alpha = get_alpha(x_local, u, param, true)->first;
      auto& work_vecs = working_storage();
      calculate_scalar_products(alpha, M_, work_vecs);
      apply_exponential(work_vecs);
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      for (size_t ll = 0; ll < quad_weights_.size(); ++ll) {
        const auto factor = work_vecs[ll] * quad_weights_[ll] * quad_points_[ll][col];
        for (size_t ii = 0; ii < dimRange; ++ii)
          ret[ii] += M_.get_entry(ll, ii) * factor;
      } // ll
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& /*u*/,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local auto H = XT::Common::make_unique<MatrixType>();
      calculate_hessian(alpha, M_, *H);
      helper<dimDomain>::partial_u(M_, *H, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& /*u*/,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local auto H = XT::Common::make_unique<MatrixType>();
      calculate_hessian(alpha, M_, *H);
      partial_u_col_helper(col, M_, *H, ret);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    template <size_t domainDim = dimDomain, class anything = void>
    struct helper
    {
      static void partial_u(const BasisValuesMatrixType& M,
                            MatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        for (size_t dd = 0; dd < domainDim; ++dd)
          entropy_flux->partial_u_col_helper(dd, M, H, ret[dd], dd > 0);
      } // void partial_u(...)

      static RangeFieldType& get_ref(RangeType& ret, const size_t rr, const size_t cc)
      {
        return ret[rr][cc];
      }
    }; // class helper<...>

    template <class anything>
    struct helper<1, anything>
    {
      static void partial_u(const BasisValuesMatrixType& M,
                            MatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        entropy_flux->partial_u_col_helper(0, M, H, ret, false);
      } // void partial_u(...)

      static RangeFieldType& get_ref(RangeType& ret, const size_t rr, const size_t DXTC_DEBUG_ONLY(cc))
      {
        assert(cc == 0);
        return ret[rr];
      }
    }; // class helper<1, ...>

    void partial_u_col_helper(const size_t col,
                              const BasisValuesMatrixType& M,
                              MatrixType& H,
                              ColPartialURangeType& ret,
                              bool L_calculated = false) const
    {
      assert(col < dimDomain);
      calculate_J(M, ret, col);
      calculate_A_Binv(ret, H, L_calculated);
    } // void partial_u_col(...)

    // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
    static void calculate_A_Binv(ColPartialURangeType& A, MatrixType& B, bool L_calculated = false)
    {
      // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
      // calculate B = LL^T first
      if (!L_calculated)
        XT::LA::cholesky(B);
      VectorType tmp_vec;
      for (size_t ii = 0; ii < dimRange; ++ii) {
        // calculate C = A (L^T)^{-1} and store in B
        XT::LA::solve_lower_triangular(B, tmp_vec, A[ii]);
        // calculate ret = C L^{-1}
        XT::LA::solve_lower_triangular_transposed(B, A[ii], tmp_vec);
      } // ii
    } // void calculate_A_Binv(...)

    void calculate_hessian(const VectorType& alpha, const BasisValuesMatrixType& M, MatrixType& H) const
    {
      std::fill(H.begin(), H.end(), 0.);
      auto& work_vec = working_storage();
      calculate_scalar_products(alpha, M, work_vec);
      apply_exponential(work_vec);
      const size_t num_quad_points = quad_weights_.size();
      // matrix is symmetric, we only use lower triangular part
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        auto factor_ll = work_vec[ll] * quad_weights_[ll];
        const auto* basis_ll = M.get_ptr(ll);
        for (size_t ii = 0; ii < dimRange; ++ii) {
          auto* H_row = &(H[ii][0]);
          const auto factor_ll_ii = basis_ll[ii] * factor_ll;
          if (!XT::Common::is_zero(factor_ll_ii)) {
            for (size_t kk = 0; kk <= ii; ++kk) {
              H_row[kk] += basis_ll[kk] * factor_ll_ii;
            } // kk
          }
        } // ii
      } // ll
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) values
    void calculate_J(const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<RangeFieldType, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimRangeCols);
      const auto& work_vecs = working_storage();
      std::fill(J_dd.begin(), J_dd.end(), 0);
      const size_t num_quad_points = quad_points_.size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto factor_ll = work_vecs[ll] * quad_points_[ll][dd] * quad_weights_[ll];
        const auto* basis_ll = M.get_ptr(ll);
        for (size_t ii = 0; ii < dimRange; ++ii) {
          const auto factor_ll_ii = factor_ll * basis_ll[ii];
          if (!XT::Common::is_zero(factor_ll_ii)) {
            for (size_t kk = 0; kk <= ii; ++kk)
              J_dd[ii][kk] += basis_ll[kk] * factor_ll_ii;
          }
        } // ii
      } // ll
      // symmetric update for upper triangular part of J
      for (size_t mm = 0; mm < dimRange; ++mm)
        for (size_t nn = mm + 1; nn < dimRange; ++nn)
          J_dd[mm][nn] = J_dd[nn][mm];
    } // void calculate_J(...)

    void change_basis(const VectorType& beta_in,
                      VectorType& v_k,
                      BasisValuesMatrixType& P_k,
                      MatrixType& T_k,
                      VectorType& g_k,
                      VectorType& beta_out,
                      MatrixType& H) const
    {
      calculate_hessian(beta_in, P_k, H);
      XT::LA::cholesky(H);
      const auto& L = H;
      thread_local std::unique_ptr<MatrixType> tmp_mat = std::make_unique<MatrixType>();
      *tmp_mat = T_k;
      rightmultiply(T_k, *tmp_mat, L);
      L.mtv(beta_in, beta_out);
      StateRangeType tmp_vec;
      XT::LA::solve_lower_triangular(L, tmp_vec, v_k);
      v_k = tmp_vec;
      apply_inverse_matrix(L, P_k);
      calculate_vector_integral(beta_out, P_k, P_k, g_k, false);
      g_k -= v_k;
    } // void change_basis(...)

    const BasisfunctionType& basis_functions_;
    const std::vector<DomainType>& quad_points_;
    const std::vector<RangeFieldType>& quad_weights_;
    const BasisValuesMatrixType& M_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const MatrixType& T_minus_one_;
    const std::string name_;
    // constructor)
    LocalCacheType& cache_;
    AlphaStorageType& alpha_storage_;
    std::mutex& mutex_;
    RealizabilityHelper<BasisfunctionType, true, true> realizability_helper_;
  }; // class Localfunction

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           quad_points_,
                                           quad_weights_,
                                           M_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           *T_minus_one_,
                                           cache_[index],
                                           alpha_storage_[index],
                                           mutexes_[index]
#if HAVE_CLP
                                           ,
                                           lp_);
#else
    );
#endif
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& /*u_i*/,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType& u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& /*param*/,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    const bool boundary = static_cast<bool>(param_neighbor.get("boundary")[0]);
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_stored_alpha(x_local_entity);
    StateRangeType alpha_j;
    if (boundary)
      alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, true)->first;
    else
      alpha_j = local_function_neighbor->get_stored_alpha(x_local_neighbor);
    thread_local FieldVector<std::vector<RangeFieldType>, 2> work_vecs;
    work_vecs[0].resize(quad_points_.size());
    work_vecs[1].resize(quad_points_.size());
    local_function_entity->calculate_scalar_products(alpha_i, M_, work_vecs[0]);
    local_function_entity->calculate_scalar_products(alpha_j, M_, work_vecs[1]);
    StateRangeType ret(0);
    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      const auto position = quad_points_[ll][dd];
      RangeFieldType factor = position * n_ij[dd] > 0. ? std::exp(work_vecs[0][ll]) : std::exp(work_vecs[1][ll]);
      factor *= quad_weights_[ll] * position;
      const auto* basis_ll = M_.get_ptr(ll);
      for (size_t ii = 0; ii < dimRange; ++ii)
        ret[ii] += basis_ll[ii] * factor;
    } // ll
    ret *= n_ij[dd];
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

  const BasisfunctionType& basis_functions() const
  {
    return basis_functions_;
  }

private:
  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  std::vector<DomainType> quad_points_;
  std::vector<RangeFieldType> quad_weights_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const std::unique_ptr<MatrixType> T_minus_one_;
  const std::string name_;
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see
  // constructor)
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<AlphaStorageType> alpha_storage_;
  mutable std::vector<std::mutex> mutexes_;
#if HAVE_CLP
  mutable XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>> lp_;
#endif
};


#if 1
/**
 * Specialization for DG basis
 */
template <class GridLayerImp, class U, size_t domainDim, size_t dimRange_or_refinements>
class EntropyBasedLocalFlux<
    PartialMomentBasis<typename U::DomainFieldType, domainDim, typename U::RangeFieldType, dimRange_or_refinements, 1>,
    GridLayerImp,
    U>
  : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                           typename U::DomainFieldType,
                                                           GridLayerImp::dimension,
                                                           U,
                                                           0,
                                                           typename U::RangeFieldType,
                                                           U::dimRange,
                                                           GridLayerImp::dimension>
{
  using BaseType =
      typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                               typename U::DomainFieldType,
                                                               GridLayerImp::dimension,
                                                               U,
                                                               0,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               GridLayerImp::dimension>;
  using ThisType = EntropyBasedLocalFlux;

public:
  using GridLayerType = GridLayerImp;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using typename BaseType::StateType;
  static const size_t block_size = (dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  using BlockMatrixType = XT::Common::BlockedFieldMatrix<RangeFieldType, num_blocks, block_size>;
  using LocalMatrixType = typename BlockMatrixType::BlockType;
  using BlockVectorType = XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>;
  using LocalVectorType = typename BlockVectorType::BlockType;
  using BasisValuesMatrixType = FieldVector<XT::LA::CommonDenseMatrix<RangeFieldType>, num_blocks>;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, dimDomain>;
  using QuadraturePointsType =
      FieldVector<std::vector<DomainType, boost::alignment::aligned_allocator<DomainType, 64>>, num_blocks>;
  using QuadratureWeightsType =
      FieldVector<std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>, num_blocks>;
  using BasisfunctionType =
      PartialMomentBasis<DomainFieldType, dimDomain, RangeFieldType, dimRange_or_refinements, 1, dimDomain>;
  using AlphaReturnType = std::pair<BlockVectorType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, BlockVectorType>;
  using AlphaStorageType = std::map<DomainType, BlockVectorType, XT::Common::VectorFloatLess>;
  using TemporaryVectorType = std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>;
  using TemporaryVectorsType = FieldVector<TemporaryVectorType, num_blocks>;
  static const size_t cache_size = 4 * dimDomain + 2;

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using LocalfunctionType::dimRange;
    using typename LocalfunctionType::ColPartialURangeType;
    using typename LocalfunctionType::ColRangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const QuadraturePointsType& quad_points,
                  const QuadratureWeightsType& quad_weights,
                  const BasisValuesMatrixType& M,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const LocalMatrixType& T_minus_one,
                  LocalCacheType& cache,
                  AlphaStorageType& alpha_storage,
                  std::mutex& mutex)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , quad_points_(quad_points)
      , quad_weights_(quad_weights)
      , M_(M)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , T_minus_one_(T_minus_one)
      , cache_(cache)
      , alpha_storage_(alpha_storage)
      , mutex_(mutex)
    {}

    using LocalfunctionType::entity;

    // temporary vectors to store inner products and exponentials
    TemporaryVectorsType& working_storage() const
    {
      thread_local TemporaryVectorsType work_vecs;
      for (size_t jj = 0; jj < num_blocks; ++jj)
        work_vecs[jj].resize(quad_points_[jj].size());
      return work_vecs;
    }

    void copy_basis_matrix(const BasisValuesMatrixType& source_mat, BasisValuesMatrixType& range_mat) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        range_mat[jj].backend() = source_mat[jj].backend();
    }

    void calculate_scalar_products_block(const size_t jj,
                                         const LocalVectorType& beta_in,
                                         const XT::LA::CommonDenseMatrix<RangeFieldType>& M,
                                         TemporaryVectorType& scalar_products) const
    {
      const size_t num_quad_points = quad_points_[jj].size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto* basis_ll = M.get_ptr(ll);
        scalar_products[ll] = std::inner_product(beta_in.begin(), beta_in.end(), basis_ll, 0.);
      } // ll
    }

    void calculate_scalar_products(const BlockVectorType& beta_in,
                                   const BasisValuesMatrixType& M,
                                   TemporaryVectorsType& scalar_products) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        calculate_scalar_products_block(jj, beta_in.block(jj), M[jj], scalar_products[jj]);
    }

    void apply_exponential(TemporaryVectorType& values) const
    {
      assert(values.size() < std::numeric_limits<int>::max());
      XT::Common::Mkl::exp(static_cast<int>(values.size()), values.data(), values.data());
    }

    void apply_exponential(TemporaryVectorsType& values) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        apply_exponential(values[jj]);
    }

    // calculate ret = \int (exp(beta_in * m))
    RangeFieldType calculate_scalar_integral(const BlockVectorType& beta_in, const BasisValuesMatrixType& M) const
    {
      auto& work_vecs = working_storage();
      calculate_scalar_products(beta_in, M, work_vecs);
      apply_exponential(work_vecs);
      RangeFieldType ret(0.);
      for (size_t jj = 0; jj < num_blocks; ++jj)
        ret += std::inner_product(
            quad_weights_[jj].begin(), quad_weights_[jj].end(), work_vecs[jj].begin(), RangeFieldType(0.));
      return ret;
    }

    // calculate ret = \int (m1 exp(beta_in * m2))
    void calculate_vector_integral_block(const size_t jj,
                                         const LocalVectorType& beta_in,
                                         const XT::LA::CommonDenseMatrix<RangeFieldType>& M1,
                                         const XT::LA::CommonDenseMatrix<RangeFieldType>& M2,
                                         LocalVectorType& ret) const
    {
      auto& work_vec = working_storage()[jj];
      calculate_scalar_products_block(jj, beta_in, M2, work_vec);
      apply_exponential(work_vec);
      std::fill(ret.begin(), ret.end(), 0.);
      const size_t num_quad_points = quad_weights_[jj].size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto factor = work_vec[ll] * quad_weights_[jj][ll];
        const auto* basis_ll = M1.get_ptr(ll);
        for (size_t ii = 0; ii < block_size; ++ii)
          ret[ii] += basis_ll[ii] * factor;
      } // ll
    }

    // calculate ret = \int (m1 exp(beta_in * m2))
    void calculate_vector_integral(const BlockVectorType& beta_in,
                                   const BasisValuesMatrixType& M1,
                                   const BasisValuesMatrixType& M2,
                                   BlockVectorType& ret) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        calculate_vector_integral_block(jj, beta_in.block(jj), M1[jj], M2[jj], ret.block(jj));
    }

    void copy_transposed(const LocalMatrixType& T_k, LocalMatrixType& T_k_trans) const
    {
      for (size_t ii = 0; ii < block_size; ++ii)
        for (size_t kk = 0; kk <= ii; ++kk)
          T_k_trans[kk][ii] = T_k[ii][kk];
    }

    void apply_inverse_matrix_block(const size_t jj,
                                    const LocalMatrixType& T_k,
                                    XT::LA::CommonDenseMatrix<RangeFieldType>& M) const
    {
      const size_t num_quad_points = quad_points_[jj].size();
      if (block_size == 2) {
        const auto T_00_inv = 1 / T_k[0][0];
        const auto T_11_inv = 1 / T_k[1][1];
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          auto* basis_ll = M.get_ptr(ll);
          basis_ll[0] *= T_00_inv;
          basis_ll[1] = (basis_ll[1] - T_k[1][0] * basis_ll[0]) * T_11_inv;
        }
      } else if (block_size == 4) {
        FieldVector<RangeFieldType, 4> diag_inv;
        for (size_t ii = 0; ii < 4; ++ii)
          diag_inv[ii] = 1. / T_k[ii][ii];
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          auto* basis_ll = M.get_ptr(ll);
          basis_ll[0] *= diag_inv[0];
          basis_ll[1] = (basis_ll[1] - T_k[1][0] * basis_ll[0]) * diag_inv[1];
          basis_ll[2] = (basis_ll[2] - T_k[2][0] * basis_ll[0] - T_k[2][1] * basis_ll[1]) * diag_inv[2];
          basis_ll[3] =
              (basis_ll[3] - T_k[3][0] * basis_ll[0] - T_k[3][1] * basis_ll[1] - T_k[3][2] * basis_ll[2]) * diag_inv[3];
        }
      } else {
#if HAVE_MKL || HAVE_CBLAS
        thread_local LocalMatrixType T_k_trans(0.);
        assert(num_quad_points < std::numeric_limits<int>::max());
        // Calculate the transpose here first as this is much faster than passing the matrix to dtrsm and using
        // CblasTrans
        copy_transposed(T_k, T_k_trans);
        XT::Common::Blas::dtrsm(XT::Common::Blas::row_major(),
                                XT::Common::Blas::right(),
                                XT::Common::Blas::upper(),
                                XT::Common::Blas::no_trans(),
                                XT::Common::Blas::non_unit(),
                                static_cast<int>(num_quad_points),
                                block_size,
                                1.,
                                &(T_k_trans[0][0]),
                                block_size,
                                M.data(),
                                block_size);
#else
        LocalVectorType tmp_vec, tmp_vec2;
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          std::copy_n(M.get_ptr(ll), block_size, tmp_vec.begin());
          XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
          std::copy_n(tmp_vec2.begin(), block_size, M.get_ptr(ll));
        }
#endif
      }
    }

    void apply_inverse_matrix(const BlockMatrixType& T_k, BasisValuesMatrixType& M) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        apply_inverse_matrix_block(jj, T_k.block(jj), M[jj]);
    }

    void store_alpha(const DomainType& x_local, const BlockVectorType& alpha)
    {
      alpha_storage_[x_local] = alpha;
    }

    BlockVectorType get_stored_alpha(const DomainType& x_local) const
    {
      return alpha_storage_.at(x_local);
    }

    template <class GridLayerType>
    void center_results_to_intersections(const GridLayerType& grid_layer)
    {
      const auto center = entity().geometry().local(entity().geometry().center());
      const auto center_alpha = get_stored_alpha(center);
      for (const auto& intersection : Dune::intersections(grid_layer, entity()))
        store_alpha(entity().geometry().local(intersection.geometry().center()), center_alpha);
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& x_local,
                                               const StateRangeType& u_in,
                                               const XT::Common::Parameter& param,
                                               const bool regularize) const
    {
      const bool boundary = static_cast<bool>(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      auto ret = std::make_unique<AlphaReturnType>();
      auto v_in = std::make_unique<StateRangeType>();
      mutex_.lock();
      if (boundary)
        cache_.set_capacity(cache_size + dimDomain);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u_in);
      if (!(density > 0. || !(basis_functions_.min_density(u_in) > 0.)) || std::isinf(density)) {
        mutex_.unlock();
        DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
      }
      auto u_prime_in = std::make_unique<StateRangeType>(u_in / density);
      auto alpha_iso_in = std::make_unique<StateRangeType>(basis_functions_.alpha_iso());
      auto alpha_iso = std::make_unique<BlockVectorType>(*alpha_iso_in);

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(*u_prime_in);
      if (cache_iterator != cache_.end()
          && XT::Common::FloatCmp::eq(cache_iterator->first, *u_prime_in, 1e-14, 1e-14)) {
        const auto& alpha_prime = cache_iterator->second;
        ret->first = *alpha_iso;
        ret->first *= std::log(density);
        ret->first += alpha_prime;
        ret->second = 0.;
        alpha_storage_[x_local] = ret->first;
        mutex_.unlock();
        return ret;
      } else {
        RangeFieldType tau_prime = std::min(
            tau_ / ((1 + std::sqrt(dimRange) * u_prime_in->two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);

        // calculate moment vector for isotropic distribution
        auto u_iso_dyn = basis_functions_.u_iso();
        auto u_iso_in = std::make_unique<StateRangeType>();
        for (size_t ii = 0; ii < dimRange; ++ii)
          (*u_iso_in)[ii] = u_iso_dyn[ii];

        // define further variables
        auto g_k = std::make_unique<BlockVectorType>();
        auto beta_out = std::make_unique<BlockVectorType>();
        auto v = std::make_unique<BlockVectorType>();
        thread_local auto T_k = XT::Common::make_unique<BlockMatrixType>();
        auto beta_in =
            std::make_unique<BlockVectorType>(cache_iterator != cache_.end() ? cache_iterator->second : *alpha_iso);

        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence) {
          // regularize u
          *v_in = *u_prime_in;
          if (r > 0.) {
            *beta_in = *alpha_iso;
            // calculate v = (1-r) u + r u_iso
            // use alpha_iso_in as storage for u_iso_in * r
            *v_in *= (1 - r);
            *alpha_iso_in = *u_iso_in;
            *alpha_iso_in *= r;
            *v_in += *alpha_iso_in;
          }
          *v = *v_in;
          for (size_t jj = 0; jj < num_blocks; ++jj)
            T_k->block(jj) = T_minus_one_;
          // calculate T_k u
          auto v_k = std::make_unique<BlockVectorType>(*v);
          // calculate values of basis p = S_k m
          thread_local BasisValuesMatrixType P_k(XT::LA::CommonDenseMatrix<RangeFieldType>(0, 0, 0., 0));
          copy_basis_matrix(M_, P_k);
          // calculate f_0
          RangeFieldType f_k = calculate_scalar_integral(*beta_in, P_k) - *beta_in * *v_k;

          thread_local auto H = XT::Common::make_unique<BlockMatrixType>(0.);

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used or cholesky decomposition fails
            if (kk > k_0_ && r < r_max)
              break;
            try {
              change_basis(*beta_in, *v_k, P_k, *T_k, *g_k, *beta_out, *H);
            } catch (const Dune::MathError&) {
              if (r < r_max)
                break;
              mutex_.unlock();
              //              std::cerr << "Failed to converge for " << XT::Common::to_string(u_in, 15) << " with
              //              density "
              //                        << XT::Common::to_string(density, 15) << " at position "
              //                        << XT::Common::to_string(entity().geometry().center(), 15) << " due to too many
              //                        iterations!"
              //                        << std::endl;
              DUNE_THROW(Dune::MathError, "Failure to converge!");
            }

            // calculate descent direction d_k;
            thread_local auto d_k = std::make_unique<BlockVectorType>();
            *d_k = *g_k;
            *d_k *= -1;

            // Calculate stopping criteria (in original basis). Variables with _k are in current basis, without k in
            // original basis.
            thread_local auto alpha_tilde = std::make_unique<BlockVectorType>();
            thread_local auto alpha_prime = std::make_unique<BlockVectorType>();
            thread_local auto u_alpha_tilde = std::make_unique<BlockVectorType>();
            thread_local auto u_alpha_prime = std::make_unique<BlockVectorType>();
            thread_local auto d_alpha_tilde = std::make_unique<BlockVectorType>();
            thread_local auto g_alpha_tilde = std::make_unique<BlockVectorType>();
            thread_local auto u_eps_diff = std::make_unique<BlockVectorType>();
            // convert everything to original basis
            for (size_t jj = 0; jj < num_blocks; ++jj) {
              XT::LA::solve_lower_triangular_transposed(T_k->block(jj), alpha_tilde->block(jj), beta_out->block(jj));
              XT::LA::solve_lower_triangular_transposed(T_k->block(jj), d_alpha_tilde->block(jj), d_k->block(jj));
            } // jj
            calculate_vector_integral(*alpha_tilde, M_, M_, *u_alpha_tilde);
            *g_alpha_tilde = *u_alpha_tilde;
            *g_alpha_tilde -= *v;
            auto density_tilde = basis_functions_.density(*u_alpha_tilde);
            if (!(density_tilde > 0.) || !(basis_functions_.min_density(*u_alpha_tilde) > 0.)
                || std::isinf(density_tilde))
              break;
            *alpha_prime = *alpha_iso;
            *alpha_prime *= -std::log(density_tilde);
            *alpha_prime += *alpha_tilde;
            calculate_vector_integral(*alpha_prime, M_, M_, *u_alpha_prime);
            *u_eps_diff = *u_alpha_prime;
            *u_eps_diff *= -(1 - epsilon_gamma_);
            *u_eps_diff += *v;
            if (g_alpha_tilde->two_norm() < tau_prime
                && 1 - epsilon_gamma_ < std::exp(d_alpha_tilde->one_norm() + std::abs(std::log(density_tilde)))
                && helper<dimDomain>::is_realizable(*u_eps_diff, basis_functions_)) {
              ret->first = *alpha_iso;
              ret->first *= std::log(density);
              ret->first += *alpha_prime;
              ret->second = r;
              cache_.insert(*v_in, *alpha_prime);
              alpha_storage_[x_local] = ret->first;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              *beta_in = *beta_out;
              // backtracking line search
              while (pure_newton >= 2 || zeta_k > epsilon_ * beta_out->two_norm() / d_k->two_norm()) {
                thread_local auto beta_new = std::make_unique<BlockVectorType>();
                *beta_new = *d_k;
                *beta_new *= zeta_k;
                *beta_new += *beta_out;
                RangeFieldType f = calculate_scalar_integral(*beta_new, P_k) - *beta_new * *v_k;
                if (pure_newton >= 2 || f <= f_k + xi_ * zeta_k * (*g_k * *d_k)) {
                  *beta_in = *beta_new;
                  f_k = f;
                  pure_newton = 0;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              if (zeta_k <= epsilon_ * beta_out->two_norm() / d_k->two_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)
        mutex_.unlock();
        //        std::cerr << "Failed to converge for " << XT::Common::to_string(u_in, 15) << " with density "
        //                  << XT::Common::to_string(density, 15) << " at position "
        //                  << XT::Common::to_string(entity().geometry().center(), 15) << " due to too many iterations!"
        //                  << std::endl;
        DUNE_THROW(MathError, "Failed to converge");
      } // else ( value has not been calculated before )

    outside_all_loops:
      mutex_.unlock();
      return ret;
    }

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      ColRangeType col_ret;
      const auto alpharet = get_alpha(x_local, u, param, true);
      const auto& alpha = alpharet->first;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        evaluate_col_helper(dd, col_ret, alpha);
        for (size_t ii = 0; ii < dimRange; ++ii)
          helper<dimDomain>::get_ref(ret, ii, dd) = col_ret[ii];
      } // dd
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      const auto alpharet = get_alpha(x_local, u, param, true);
      const auto& alpha = alpharet->first;
      evaluate_col_helper(col, ret, alpha);
    }

    void evaluate_col_helper(const size_t col, ColRangeType& ret, const BlockVectorType& alpha) const
    {
      std::fill(ret.begin(), ret.end(), 0.);
      auto& work_vecs = working_storage();
      calculate_scalar_products(alpha, M_, work_vecs);
      apply_exponential(work_vecs);
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = block_size * jj;
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto factor = work_vecs[jj][ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][col];
          for (size_t ii = 0; ii < block_size; ++ii)
            ret[offset + ii] += M_[jj].get_entry(ll, ii) * factor;
        } // ll
      } // jj
    } // void evaluate_col_helper(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& /*u*/,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local auto H = XT::Common::make_unique<BlockMatrixType>();
      calculate_hessian(alpha, M_, *H);
      helper<dimDomain>::partial_u(M_, *H, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& /*u*/,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local auto H = XT::Common::make_unique<BlockMatrixType>();
      calculate_hessian(alpha, M_, *H);
      helper<dimDomain>::partial_u_col(col, M_, *H, ret, this);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

    template <size_t dim, class anything = void>
    struct helper;

    template <class anything>
    struct helper<1, anything>
    {
      static void partial_u(const BasisValuesMatrixType& M,
                            BlockMatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(M, ret, 0);
        entropy_flux->calculate_A_Binv(ret, H);
      } // void partial_u(...)

      static void partial_u_col(const size_t DXTC_DEBUG_ONLY(col),
                                const BasisValuesMatrixType& M,
                                BlockMatrixType& H,
                                ColPartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        assert(col == 0);
        partial_u(M, H, ret, entropy_flux);
      } // void partial_u(...)

      static RangeFieldType& get_ref(RangeType& ret, const size_t rr, const size_t DXTC_DEBUG_ONLY(cc))
      {
        assert(cc == 0);
        return ret[rr];
      }

      static void calculate_plane_coefficients(const BasisfunctionType& /*basis_functions*/) {}

      static bool is_realizable(const BlockVectorType& u, const BasisfunctionType& basis_functions)
      {
        for (size_t jj = 0; jj < num_blocks; ++jj) {
          const auto& u0 = u.block(jj)[0];
          const auto& u1 = u.block(jj)[1];
          const auto& v0 = basis_functions.triangulation()[jj];
          const auto& v1 = basis_functions.triangulation()[jj + 1];
          bool ret = (u0 >= 0) && (u1 <= v1 * u0) && (v0 * u0 <= u1);
          if (!ret)
            return false;
        } // jj
        return true;
      }
    }; // class helper<1, ...>

    template <class anything>
    struct helper<3, anything>
    {
      static void partial_u(const BasisValuesMatrixType& M,
                            BlockMatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        for (size_t dd = 0; dd < dimDomain; ++dd) {
          entropy_flux->calculate_J(M, ret[dd], dd);
          entropy_flux->calculate_A_Binv(ret[dd], H, dd > 0);
        }
      } // void partial_u(...)

      static void partial_u_col(const size_t col,
                                const BasisValuesMatrixType& M,
                                BlockMatrixType& H,
                                ColPartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(M, ret, col);
        entropy_flux->calculate_A_Binv(ret, H);
      } // void partial_u(...)

      static RangeFieldType& get_ref(RangeType& ret, const size_t rr, const size_t cc)
      {
        return ret[rr][cc];
      }

#if HAVE_QHULL
      static void calculate_plane_coefficients(const BasisfunctionType& basis_functions)
      {
        if (!basis_functions.plane_coefficients()[0].size())
          basis_functions.calculate_plane_coefficients();
      }

      static bool is_realizable(const BlockVectorType& u, const BasisfunctionType& basis_functions)
      {
        for (size_t jj = 0; jj < num_blocks; ++jj)
          for (const auto& coeff : basis_functions.plane_coefficients()[jj])
            if (!(u.block(jj) * coeff.first < coeff.second))
              return false;
        return true;
      }
#else
      static void calculate_plane_coefficients(const BasisfunctionType& /*basis_functions*/)
      {
        DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
      }

      static bool is_realizable(const BlockVectorType& /*u*/, const BasisfunctionType& /*basis_functions*/)
      {
        DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
        return false;
      }
#endif
    }; // class helper<3, ...>

    // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
    static void
    calculate_A_Binv(FieldMatrix<RangeFieldType, dimRange, dimRange>& A, BlockMatrixType& B, bool L_calculated = false)
    {
      // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
      // calculate B = LL^T first
      if (!L_calculated) {
        for (size_t jj = 0; jj < num_blocks; ++jj)
          XT::LA::cholesky(B.block(jj));
      }
      FieldVector<RangeFieldType, block_size> tmp_vec;
      FieldVector<RangeFieldType, block_size> tmp_A_row;
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        // calculate C = A (L^T)^{-1}
        const auto offset = block_size * jj;
        for (size_t ii = 0; ii < block_size; ++ii) {
          for (size_t kk = 0; kk < block_size; ++kk)
            tmp_A_row[kk] = A[offset + ii][offset + kk];
          XT::LA::solve_lower_triangular(B.block(jj), tmp_vec, tmp_A_row);
          // calculate ret = C L^{-1}
          XT::LA::solve_lower_triangular_transposed(B.block(jj), tmp_A_row, tmp_vec);
          for (size_t kk = 0; kk < block_size; ++kk)
            A[offset + ii][offset + kk] = tmp_A_row[kk];
        } // ii
      } // jj
    } // void calculate_A_Binv(...)

    void calculate_hessian(const BlockVectorType& alpha, const BasisValuesMatrixType& M, BlockMatrixType& H) const
    {
      auto& work_vec = working_storage();
      calculate_scalar_products(alpha, M, work_vec);
      apply_exponential(work_vec);
      // matrix is symmetric, we only use lower triangular part
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        std::fill(H.block(jj).begin(), H.block(jj).end(), 0.);
        const size_t num_quad_points = quad_weights_[jj].size();
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          auto factor_ll = work_vec[jj][ll] * quad_weights_[jj][ll];
          const auto* basis_ll = M[jj].get_ptr(ll);
          for (size_t ii = 0; ii < block_size; ++ii) {
            auto* H_row = &(H.block(jj)[ii][0]);
            const auto factor_ll_ii = basis_ll[ii] * factor_ll;
            for (size_t kk = 0; kk <= ii; ++kk)
              H_row[kk] += basis_ll[kk] * factor_ll_ii;
          } // ii
        } // ll
      } // jj
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) values
    void calculate_J(const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<RangeFieldType, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimDomain);
      std::fill(J_dd.begin(), J_dd.end(), 0.);
      const auto& work_vec = working_storage();
      // matrix is symmetric, we only use lower triangular part
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = jj * block_size;
        const size_t num_quad_points = quad_weights_[jj].size();
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          auto factor_ll = work_vec[jj][ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
          const auto* basis_ll = M[jj].get_ptr(ll);
          for (size_t ii = 0; ii < block_size; ++ii) {
            auto* J_row = &(J_dd[offset + ii][0]);
            const auto factor_ll_ii = basis_ll[ii] * factor_ll;
            for (size_t kk = 0; kk <= ii; ++kk)
              J_row[offset + kk] += basis_ll[kk] * factor_ll_ii;
          } // ii
        } // ll
      } // jj
      // symmetric update for upper triangular part of J
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = block_size * jj;
        for (size_t mm = 0; mm < block_size; ++mm)
          for (size_t nn = mm + 1; nn < block_size; ++nn)
            J_dd[offset + mm][offset + nn] = J_dd[offset + nn][offset + mm];
      }
    } // void calculate_J(...)

    void change_basis(const BlockVectorType& beta_in,
                      BlockVectorType& v_k,
                      BasisValuesMatrixType& P_k,
                      BlockMatrixType& T_k,
                      BlockVectorType& g_k,
                      BlockVectorType& beta_out,
                      BlockMatrixType& H) const
    {
      calculate_hessian(beta_in, P_k, H);
      FieldVector<RangeFieldType, block_size> tmp_vec;
      for (size_t jj = 0; jj < num_blocks; ++jj)
        XT::LA::cholesky(H.block(jj));
      const auto& L = H;
      T_k.rightmultiply(L);
      L.mtv(beta_in, beta_out);
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        XT::LA::solve_lower_triangular(L.block(jj), tmp_vec, v_k.block(jj));
        v_k.block(jj) = tmp_vec;
      } // jj
      apply_inverse_matrix(L, P_k);
      calculate_vector_integral(beta_out, P_k, P_k, g_k);
      g_k -= v_k;
    } // void change_basis(...)

    const BasisfunctionType& basis_functions_;
    const QuadraturePointsType& quad_points_;
    const QuadratureWeightsType& quad_weights_;
    const BasisValuesMatrixType& M_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const LocalMatrixType& T_minus_one_;
    const std::string name_;
    // constructor)
    LocalCacheType& cache_;
    AlphaStorageType& alpha_storage_;
    std::mutex& mutex_;
  }; // class Localfunction

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , M_(XT::LA::CommonDenseMatrix<RangeFieldType>())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , alpha_storage_(index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {
    XT::LA::eye_matrix(T_minus_one_);
    Localfunction::template helper<dimDomain>::calculate_plane_coefficients(basis_functions_);
    const auto& quadratures = basis_functions_.quadratures();
    assert(quadratures.size() == num_blocks);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (const auto& quad_point : quadratures[jj]) {
        quad_points_[jj].emplace_back(quad_point.position());
        quad_weights_[jj].emplace_back(quad_point.weight());
      }
    } // jj
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      while (quad_weights_[jj].size() % 8) { // align to 64 byte boundary
        quad_points_[jj].push_back(quad_points_[jj].back());
        quad_weights_[jj].push_back(0.);
      }
      M_[jj] = XT::LA::CommonDenseMatrix<RangeFieldType>(quad_points_[jj].size(), block_size, 0., 0);
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = basis_functions_.evaluate(quad_points_[jj][ll], jj);
        for (size_t ii = 0; ii < block_size; ++ii)
          M_[jj].set_entry(ll, ii, val[block_size * jj + ii]);
      } // ll
    } // jj
  }

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           quad_points_,
                                           quad_weights_,
                                           M_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           T_minus_one_,
                                           cache_[index],
                                           alpha_storage_[index],
                                           mutexes_[index]);
  }


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& /*u_i*/,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& /*param*/,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_stored_alpha(x_local_entity);
    const bool boundary = static_cast<bool>(param_neighbor.get("boundary")[0]);
    BlockVectorType alpha_j;
    if (boundary)
      alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, true)->first;
    else
      alpha_j = local_function_neighbor->get_stored_alpha(x_local_neighbor);
    thread_local FieldVector<TemporaryVectorsType, 2> work_vecs;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      work_vecs[0][jj].resize(quad_points_[jj].size());
      work_vecs[1][jj].resize(quad_points_[jj].size());
    }
    local_function_entity->calculate_scalar_products(alpha_i, M_, work_vecs[0]);
    local_function_entity->calculate_scalar_products(alpha_j, M_, work_vecs[1]);
    StateRangeType ret(0);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto position = quad_points_[jj][ll][dd];
        RangeFieldType factor =
            position * n_ij[dd] > 0. ? std::exp(work_vecs[0][jj][ll]) : std::exp(work_vecs[1][jj][ll]);
        factor *= quad_weights_[jj][ll] * position;
        for (size_t ii = 0; ii < block_size; ++ii)
          ret[offset + ii] += M_[jj].get_entry(ll, ii) * factor;
      } // ll
    } // jj
    ret *= n_ij[dd];
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

  const BasisfunctionType& basis_functions() const
  {
    return basis_functions_;
  }

private:
  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  LocalMatrixType T_minus_one_;
  const std::string name_;
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<AlphaStorageType> alpha_storage_;
  mutable std::vector<std::mutex> mutexes_;
};
#endif


#if 1
/**
 * Specialization of EntropyBasedLocalFlux for 3D Hatfunctions
 */
template <class GridLayerImp, class U, size_t refinements>
class EntropyBasedLocalFlux<
    HatFunctionMomentBasis<typename U::DomainFieldType, 3, typename U::RangeFieldType, refinements, 1>,
    GridLayerImp,
    U>
  : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                           typename U::DomainFieldType,
                                                           GridLayerImp::dimension,
                                                           U,
                                                           0,
                                                           typename U::RangeFieldType,
                                                           U::dimRange,
                                                           GridLayerImp::dimension>
{
  using BaseType =
      typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                               typename U::DomainFieldType,
                                                               GridLayerImp::dimension,
                                                               U,
                                                               0,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               GridLayerImp::dimension>;
  using ThisType = EntropyBasedLocalFlux;

public:
  using GridLayerType = GridLayerImp;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using typename BaseType::StateType;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, 3>;
  using LocalMatrixType = XT::Common::FieldMatrix<RangeFieldType, 3, 3>;
  using BasisValuesMatrixType = std::vector<std::vector<LocalVectorType>>;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, dimDomain>;
  using QuadraturePointsType = std::vector<std::vector<DomainType>>;
  using QuadratureWeightsType = std::vector<std::vector<RangeFieldType>>;
  using BasisfunctionType =
      HatFunctionMomentBasis<DomainFieldType, dimDomain, RangeFieldType, refinements, 1, dimDomain>;
#if HAVE_EIGEN
  using SparseMatrixType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::eigen_sparse>::MatrixType;
  using VectorType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::eigen_sparse>::VectorType;
#else
  using SparseMatrixType = typename XT::LA::Container<RangeFieldType, XT::LA::default_sparse_backend>::MatrixType;
  using VectorType = typename XT::LA::Container<RangeFieldType, XT::LA::default_sparse_backend>::VectorType;
#endif
  using AlphaReturnType = typename std::pair<VectorType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<VectorType, VectorType>;
  using AlphaStorageType = std::map<DomainType, VectorType, XT::Common::VectorFloatLess>;
  static constexpr size_t cache_size = 4 * dimDomain + 2;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const std::string name = "")
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , quad_points_(basis_functions_.triangulation().faces().size())
    , quad_weights_(basis_functions_.triangulation().faces().size())
    , M_(basis_functions_.triangulation().faces().size())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , alpha_storage_(index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {
    const auto& triangulation = basis_functions_.triangulation();
    const auto& vertices = triangulation.vertices();
    const auto& faces = triangulation.faces();
    assert(vertices.size() == dimRange);
    // create pattern
    XT::LA::SparsityPatternDefault pattern(dimRange);
    for (size_t vertex_index = 0; vertex_index < dimRange; ++vertex_index) {
      const auto& vertex = vertices[vertex_index];
      const auto& adjacent_faces = triangulation.get_face_indices(vertex->position());
      for (const auto& face_index : adjacent_faces) {
        const auto& face = faces[face_index];
        assert(face->vertices().size() == 3);
        for (size_t jj = 0; jj < 3; ++jj)
          pattern.insert(vertex_index, face->vertices()[jj]->index());
      }
    }
    pattern.sort();
    pattern_ = pattern;
    // store basis evaluations
    const auto& quadratures = basis_functions_.quadratures();
    assert(quadratures.size() == faces.size());
    for (size_t jj = 0; jj < faces.size(); ++jj) {
      for (const auto& quad_point : quadratures[jj]) {
        quad_points_[jj].emplace_back(quad_point.position());
        quad_weights_[jj].emplace_back(quad_point.weight());
      }
    } // jj
    for (size_t jj = 0; jj < faces.size(); ++jj) {
      M_[jj] = std::vector<LocalVectorType>(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        M_[jj][ll] = basis_functions_.evaluate_on_face(quad_points_[jj][ll], jj);
    } // jj

    // resize temporary vectors to store inner products and exponentials
    work_vecs_.resize(faces.size());
    for (size_t jj = 0; jj < faces.size(); ++jj)
      work_vecs_[jj].resize(quad_points_[jj].size());
  } // constructor

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using LocalfunctionType::dimRange;
    using typename LocalfunctionType::ColPartialURangeType;
    using typename LocalfunctionType::ColRangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const QuadraturePointsType& quad_points,
                  const QuadratureWeightsType& quad_weights,
                  const BasisValuesMatrixType& M,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  LocalCacheType& cache,
                  AlphaStorageType& alpha_storage,
                  std::mutex& mutex,
                  const XT::LA::SparsityPatternDefault& pattern,
                  std::vector<std::vector<RangeFieldType>>& work_vecs)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , quad_points_(quad_points)
      , quad_weights_(quad_weights)
      , M_(M)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , cache_(cache)
      , alpha_storage_(alpha_storage)
      , mutex_(mutex)
      , pattern_(pattern)
      , work_vecs_(work_vecs)
    {}

    static bool is_realizable(const VectorType& u)
    {
      for (const auto& u_i : u)
        if (!(u_i > 0.) || std::isinf(u_i))
          return false;
      return true;
    }

    void store_alpha(const DomainType& x_local, const VectorType& alpha)
    {
      alpha_storage_[x_local] = alpha;
    }

    VectorType get_stored_alpha(const DomainType& x_local) const
    {
      return alpha_storage_.at(x_local);
    }

    using LocalfunctionType::entity;

    template <class GridLayerType>
    void center_results_to_intersections(const GridLayerType& grid_layer)
    {
      const auto center = entity().geometry().local(entity().geometry().center());
      const auto center_alpha = get_stored_alpha(center);
      for (const auto& intersection : Dune::intersections(grid_layer, entity()))
        store_alpha(entity().geometry().local(intersection.geometry().center()), center_alpha);
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& x_local,
                                               const StateRangeType& u,
                                               const XT::Common::Parameter& param,
                                               const bool regularize) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      auto ret = std::make_unique<AlphaReturnType>();
      mutex_.lock();
      if (boundary)
        cache_.set_capacity(cache_size + dimDomain);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      if (!(density > 0.) || std::isinf(density)) {
        mutex_.unlock();
        DUNE_THROW(Dune::MathError, "Negative density!");
      }
      VectorType u_prime(dimRange, 0., 0);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_prime.set_entry(ii, u[ii] / density);
      VectorType alpha_iso(dimRange, 0., 0);
      basis_functions_.alpha_iso(alpha_iso);

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u_prime, 1e-14, 1e-14)) {
        const auto& alpha_prime = cache_iterator->second;
        ret->first = alpha_iso;
        ret->first *= std::log(density);
        ret->first += alpha_prime;
        ret->second = 0.;
        alpha_storage_[x_local] = ret->first;
        mutex_.unlock();
        return ret;
      } else {
        RangeFieldType tau_prime = std::min(
            tau_ / ((1 + std::sqrt(dimRange) * u_prime.l2_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
        thread_local SparseMatrixType H(dimRange, dimRange, pattern_, 0);
        thread_local auto solver = XT::LA::make_solver(H);

        // calculate moment vector for isotropic distribution
        VectorType u_iso(dimRange, 0., 0);
        basis_functions_.u_iso(u_iso);
        VectorType alpha_k = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;
        VectorType v(dimRange, 0., 0), g_k(dimRange, 0., 0), d_k(dimRange, 0., 0), tmp_vec(dimRange, 0., 0),
            alpha_prime(dimRange);
        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence_) {
          // regularize u
          v = u_prime;
          if (r > 0) {
            alpha_k = alpha_iso;
            tmp_vec = u_iso;
            tmp_vec *= r;
            v *= 1 - r;
            v += tmp_vec;
          }

          // calculate f_0
          RangeFieldType f_k = calculate_f(alpha_k, v);

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used
            if (kk > k_0_ && r < r_max)
              break;
            // calculate gradient g
            calculate_gradient(alpha_k, v, g_k);
            // calculate Hessian H
            calculate_hessian(alpha_k, M_, H, true);
            // calculate descent direction d_k;
            tmp_vec = g_k;
            tmp_vec *= -1;
            try {
              solver.apply(tmp_vec, d_k);
            } catch (const XT::LA::Exceptions::linear_solver_failed& error) {
              if (r < r_max) {
                break;
              } else {
                mutex_.unlock();
                DUNE_THROW(XT::LA::Exceptions::linear_solver_failed,
                           "Failure to converge, solver error was: " << error.what());
              }
            }

            const auto& alpha_tilde = alpha_k;
            auto& u_alpha_tilde = tmp_vec;
            u_alpha_tilde = g_k;
            u_alpha_tilde += v;
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            if (!(density_tilde > 0.) || std::isinf(density_tilde))
              break;
            alpha_prime = alpha_iso;
            alpha_prime *= -std::log(density_tilde);
            alpha_prime += alpha_tilde;
            auto& u_eps_diff = tmp_vec;
            calculate_u(alpha_prime, u_eps_diff); // store u_alpha_prime in u_eps_diff
            u_eps_diff *= -(1 - epsilon_gamma_);
            u_eps_diff += v;
            // checking realizability is cheap so we do not need the second stopping criterion
            if (g_k.l2_norm() < tau_prime && is_realizable(u_eps_diff)) {
              ret->first = alpha_iso;
              ret->first *= std::log(density);
              ret->first += alpha_prime;
              ret->second = r;
              cache_.insert(v, alpha_prime);
              alpha_storage_[x_local] = ret->first;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              // backtracking line search
              auto& alpha_new = tmp_vec;
              while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.l2_norm() / d_k.l2_norm()) {
                // calculate alpha_new = alpha_k + zeta_k d_k
                alpha_new = d_k;
                alpha_new *= zeta_k;
                alpha_new += alpha_k;
                // calculate f(alpha_new)
                RangeFieldType f_new = calculate_f(alpha_new, v);
                if (pure_newton >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  alpha_k = alpha_new;
                  f_k = f_new;
                  pure_newton = 0;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
              if (zeta_k <= epsilon_ * alpha_k.l2_norm() / d_k.l2_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)
        mutex_.unlock();
        DUNE_THROW(MathError, "Failed to converge");
      } // else ( value has not been calculated before )

    outside_all_loops:
      mutex_.unlock();
      return ret;
    } // ... get_alpha(...)

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      ColRangeType col_ret;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        evaluate_col(dd, x_local, u, col_ret, param);
        for (size_t ii = 0; ii < dimRange; ++ii)
          ret[ii][dd] = col_ret[ii];
      } // dd
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      std::fill(ret.begin(), ret.end(), 0.);
      const auto alpha = get_alpha(x_local, u, param, true)->first;
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      LocalVectorType local_alpha, local_ret;
      const auto& triangulation = basis_functions_.triangulation();
      const auto& faces = triangulation.faces();
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        local_ret *= 0.;
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        for (size_t ii = 0; ii < 3; ++ii)
          local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          auto factor_ll = std::exp(local_alpha * basis_ll) * quad_points_[jj][ll][col] * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 3; ++ii)
            local_ret[ii] += basis_ll[ii] * factor_ll;
        } // ll (quad points)
        for (size_t ii = 0; ii < 3; ++ii)
          ret[vertices[ii]->index()] += local_ret[ii];
      } // jj (faces)
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& /*u*/,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local SparseMatrixType H(dimRange, dimRange, pattern_, 0);
      thread_local SparseMatrixType J(dimRange, dimRange, pattern_, 0);
      calculate_hessian(alpha, M_, H);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        partial_u_col_helper(dd, M_, H, J, ret[dd]);
      }
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& /*u*/,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      thread_local SparseMatrixType H(dimRange, dimRange, pattern_, 0);
      thread_local SparseMatrixType J(dimRange, dimRange, pattern_, 0);
      calculate_hessian(alpha, M_, H);
      partial_u_col_helper(col, M_, H, J, ret);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    void partial_u_col_helper(const size_t col,
                              const BasisValuesMatrixType& M,
                              SparseMatrixType& H,
                              SparseMatrixType& J,
                              ColPartialURangeType& ret) const
    {
      assert(col < dimDomain);
      calculate_J(M, J, col);
      calculate_J_Hinv(J, H, ret);
    } // void partial_u_col(...)

    // calculates ret = J H^{-1}. H is assumed to be symmetric positive definite, which gives ret^T = H^{-T} J^T =
    // H^{-1} J^T, so we just have to solve y = H^{-1} x for each row x of J
    void calculate_J_Hinv(SparseMatrixType& J, const SparseMatrixType& H, ColPartialURangeType& ret) const
    {
      thread_local VectorType solution(dimRange, 0., 0), tmp_rhs(dimRange, 0., 0);
#if HAVE_EIGEN
      typedef ::Eigen::SparseMatrix<RangeFieldType, ::Eigen::ColMajor> ColMajorBackendType;
      ColMajorBackendType colmajor_copy(H.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
#else // HAVE_EIGEN
      auto solver = XT::LA::make_solver(H);
#endif // HAVE_EIGEN
      for (size_t ii = 0; ii < dimRange; ++ii) {
        // copy row to VectorType
        for (size_t kk = 0; kk < dimRange; ++kk)
          tmp_rhs.set_entry(kk, J.get_entry(ii, kk));
          // solve
#if HAVE_EIGEN
        solution.backend() = solver.solve(tmp_rhs.backend());
#else // HAVE_EIGEN
        solver.apply(tmp_rhs, solution);
#endif
        // copy result to C
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[ii][kk] = solution.get_entry(kk);
      }
    } // void calculate_J_Hinv(...)

    RangeFieldType calculate_f(const VectorType& alpha, const VectorType& v) const
    {
      RangeFieldType ret(0.);
      XT::Common::FieldVector<RangeFieldType, 3> local_alpha;
      const auto& triangulation = basis_functions_.triangulation();
      const auto& faces = triangulation.faces();
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        for (size_t ii = 0; ii < 3; ++ii)
          local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
          ret += std::exp(local_alpha * M_[jj][ll]) * quad_weights_[jj][ll];
      } // jj (faces)
      ret -= alpha * v;
      return ret;
    } // void calculate_u(...)

    void calculate_u(const VectorType& alpha, VectorType& u) const
    {
      u *= 0.;
      LocalVectorType local_alpha, local_u;
      const auto& triangulation = basis_functions_.triangulation();
      const auto& faces = triangulation.faces();
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        local_u *= 0.;
        for (size_t ii = 0; ii < 3; ++ii)
          local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          work_vecs_[jj][ll] = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 3; ++ii)
            local_u[ii] += basis_ll[ii] * work_vecs_[jj][ll];
        } // ll (quad points)
        for (size_t ii = 0; ii < 3; ++ii)
          u.add_to_entry(vertices[ii]->index(), local_u[ii]);
      } // jj (faces)
    } // void calculate_u(...)

    void calculate_gradient(const VectorType& alpha, const VectorType& v, VectorType& g_k) const
    {
      calculate_u(alpha, g_k);
      g_k -= v;
    }

    void calculate_hessian(const VectorType& alpha,
                           const BasisValuesMatrixType& M,
                           SparseMatrixType& H,
                           const bool use_work_vecs_results = false) const
    {
      H *= 0.;
      LocalVectorType local_alpha;
      LocalMatrixType H_local(0.);
      const auto& triangulation = basis_functions_.triangulation();
      const auto& faces = triangulation.faces();
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        H_local *= 0.;
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        for (size_t ii = 0; ii < 3; ++ii)
          local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M[jj][ll];
          if (!use_work_vecs_results)
            work_vecs_[jj][ll] = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 3; ++ii)
            for (size_t kk = 0; kk < 3; ++kk)
              H_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * work_vecs_[jj][ll];
        } // ll (quad points)
        for (size_t ii = 0; ii < 3; ++ii)
          for (size_t kk = 0; kk < 3; ++kk)
            H.add_to_entry(vertices[ii]->index(), vertices[kk]->index(), H_local[ii][kk]);
      } // jj (faces)
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) values
    void calculate_J(const BasisValuesMatrixType& M, SparseMatrixType& J_dd, const size_t dd) const
    {
      assert(dd < dimRangeCols);
      J_dd *= 0.;
      LocalMatrixType J_local(0.);
      const auto& triangulation = basis_functions_.triangulation();
      const auto& faces = triangulation.faces();
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        J_local *= 0.;
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M[jj][ll];
          for (size_t ii = 0; ii < 3; ++ii)
            for (size_t kk = 0; kk < 3; ++kk)
              J_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * work_vecs_[jj][ll] * quad_points_[jj][ll][dd];
        } // ll (quad points)
        for (size_t ii = 0; ii < 3; ++ii)
          for (size_t kk = 0; kk < 3; ++kk)
            J_dd.add_to_entry(vertices[ii]->index(), vertices[kk]->index(), J_local[ii][kk]);
      } // jj (faces)
    } // void calculate_J(...)

    const BasisfunctionType& basis_functions_;
    const QuadraturePointsType& quad_points_;
    const QuadratureWeightsType& quad_weights_;
    const BasisValuesMatrixType& M_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const std::string name_;
    // constructor)
    LocalCacheType& cache_;
    AlphaStorageType& alpha_storage_;
    std::mutex& mutex_;
    const XT::LA::SparsityPatternDefault& pattern_;
    std::vector<std::vector<RangeFieldType>>& work_vecs_;
  }; // class Localfunction

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           quad_points_,
                                           quad_weights_,
                                           M_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           cache_[index],
                                           alpha_storage_[index],
                                           mutexes_[index],
                                           pattern_,
                                           work_vecs_);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& /*u_i*/,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType& u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& /*param*/,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    const bool boundary = static_cast<bool>(param_neighbor.get("boundary")[0]);
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_stored_alpha(x_local_entity);
    VectorType alpha_j;
    if (boundary)
      alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, true)->first;
    else
      alpha_j = local_function_neighbor->get_stored_alpha(x_local_neighbor);
    StateRangeType ret(0);
    const auto& triangulation = basis_functions_.triangulation();
    const auto& faces = triangulation.faces();
    LocalVectorType local_alpha_i, local_alpha_j, local_ret;
    for (size_t jj = 0; jj < faces.size(); ++jj) {
      local_ret *= 0.;
      const auto& face = faces[jj];
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii) {
        local_alpha_i[ii] = alpha_i.get_entry(vertices[ii]->index());
        local_alpha_j[ii] = alpha_j.get_entry(vertices[ii]->index());
      }
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M_[jj][ll];
        const auto position = quad_points_[jj][ll][dd];
        RangeFieldType factor =
            position * n_ij[dd] > 0. ? std::exp(local_alpha_i * basis_ll) : std::exp(local_alpha_j * basis_ll);
        factor *= quad_weights_[jj][ll] * position;
        for (size_t ii = 0; ii < 3; ++ii)
          local_ret[ii] += basis_ll[ii] * factor;
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        ret[vertices[ii]->index()] += local_ret[ii];
    } // jj (faces)
    ret *= n_ij[dd];
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

  const BasisfunctionType& basis_functions() const
  {
    return basis_functions_;
  }

private:
  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const std::string name_;
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<AlphaStorageType> alpha_storage_;
  mutable std::vector<std::mutex> mutexes_;
  XT::LA::SparsityPatternDefault pattern_;
  static thread_local std::vector<std::vector<RangeFieldType>> work_vecs_;
};
#endif

template <class GridLayerImp, class U, size_t refinements>
thread_local std::vector<std::vector<typename EntropyBasedLocalFlux<
    HatFunctionMomentBasis<typename U::DomainFieldType, 3, typename U::RangeFieldType, refinements, 1>,
    GridLayerImp,
    U>::RangeFieldType>>
    EntropyBasedLocalFlux<
        HatFunctionMomentBasis<typename U::DomainFieldType, 3, typename U::RangeFieldType, refinements, 1>,
        GridLayerImp,
        U>::work_vecs_;

#if 0
/**
 * Specialization of EntropyBasedLocalFlux for 1D Hatfunctions (no change of basis, analytic integrals + Taylor)
 */
template <class GridLayerImp, class U>
class EntropyBasedLocalFlux<
    HatFunctionMomentBasis<typename U::DomainFieldType, 1, typename U::RangeFieldType, U::dimRange, 1, 1>,
    GridLayerImp,
    U>
  : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                           typename U::DomainFieldType,
                                                           GridLayerImp::dimension,
                                                           U,
                                                           0,
                                                           typename U::RangeFieldType,
                                                           U::dimRange,
                                                           1>
{
  using BaseType =
      typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                               typename U::DomainFieldType,
                                                               GridLayerImp::dimension,
                                                               U,
                                                               0,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               1>;
  using ThisType = EntropyBasedLocalFlux;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using typename BaseType::StateType;
  using BasisfunctionType = HatFunctionMomentBasis<DomainFieldType, 1, RangeFieldType, dimRange, 1, 1>;
  using GridLayerType = GridLayerImp;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, 1>;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using AlphaReturnType = typename std::pair<StateRangeType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, StateRangeType>;
  using AlphaStorageType = std::map<DomainType, StateRangeType, XT::Common::VectorFloatLess>;
  static const size_t cache_size = 4 * dimDomain + 2;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const RangeFieldType taylor_tol = 0.1,
      const size_t max_taylor_order = 200,
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , v_points_(basis_functions_.triangulation())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , taylor_tol_(taylor_tol)
    , max_taylor_order_(max_taylor_order)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , alpha_storage_(index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {}

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColPartialURangeType;
    using typename LocalfunctionType::ColRangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const std::vector<RangeFieldType>& v_points,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const RangeFieldType taylor_tol,
                  const size_t max_taylor_order,
                  LocalCacheType& cache,
                  AlphaStorageType& alpha_storage,
                  std::mutex& mutex)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , v_points_(v_points)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , taylor_tol_(taylor_tol)
      , max_taylor_order_(max_taylor_order)
      , cache_(cache)
      , alpha_storage_(alpha_storage)
      , mutex_(mutex)
    {}

    using LocalfunctionType::entity;

    static bool is_realizable(const RangeType& u)
    {
      for (const auto& u_i : u)
        if (!(u_i > 0.) || std::isinf(u_i))
          return false;
      return true;
    }

    void store_alpha(const DomainType& x_local, const StateRangeType& alpha)
    {
      alpha_storage_[x_local] = alpha;
    }

    StateRangeType get_stored_alpha(const DomainType& x_local) const
    {
      return alpha_storage_.at(x_local);
    }

    template <class GridLayerType>
    void center_results_to_intersections(const GridLayerType& grid_layer)
    {
      const auto center = entity().geometry().local(entity().geometry().center());
      const auto center_alpha = get_stored_alpha(center);
      for (const auto& intersection : Dune::intersections(grid_layer, entity()))
        store_alpha(entity().geometry().local(intersection.geometry().center()), center_alpha);
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& x_local,
                                               const StateRangeType& u,
                                               const XT::Common::Parameter& param,
                                               const bool regularize) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      auto ret = std::make_unique<AlphaReturnType>();
      mutex_.lock();
      if (boundary)
        cache_.set_capacity(cache_size + dimDomain);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      if (!(density > 0.) || std::isinf(density)) {
        mutex_.unlock();
        DUNE_THROW(Dune::MathError, "Negative density!");
      }
      RangeType u_prime = u / density;
      RangeType alpha_iso = basis_functions_.alpha_iso();

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u_prime, 1e-14, 1e-14)) {
        const auto alpha_prime = cache_iterator->second;
        ret->first = alpha_prime + alpha_iso * std::log(density);
        ret->second = 0.;
        alpha_storage_[x_local] = ret->first;
        mutex_.unlock();
        return ret;
      } else {
        RangeFieldType tau_prime = std::min(
            tau_ / ((1 + std::sqrt(dimRange) * u_prime.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
        // The hessian H is always symmetric and tridiagonal, so we only need to store the diagonal and subdiagonal
        // elements
        RangeType H_diag;
        FieldVector<RangeFieldType, dimRange - 1> H_subdiag;

        // calculate moment vector for isotropic distribution
        RangeType u_iso = basis_functions_.u_iso();
        RangeType v;
        RangeType alpha_k = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;
        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence_) {
          // regularize u
          v = u_prime;
          if (r > 0) {
            alpha_k = alpha_iso;
            RangeType r_times_u_iso(u_iso);
            r_times_u_iso *= r;
            v *= 1 - r;
            v += r_times_u_iso;
          }

          // calculate f_0
          RangeFieldType f_k = calculate_f(alpha_k, v);

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used
            if (kk > k_0_ && r < r_max)
              break;
            // calculate gradient g
            RangeType g_k = calculate_gradient(alpha_k, v);
            // calculate Hessian H
            calculate_hessian(alpha_k, H_diag, H_subdiag);
            // calculate descent direction d_k;
            RangeType d_k(0), minus_g_k(g_k);
            minus_g_k *= -1;
            try {
              d_k = minus_g_k;
              XT::LA::solve_sym_tridiag_posdef(H_diag, H_subdiag, d_k);
            } catch (const Dune::MathError&) {
              if (r < r_max) {
                break;
              } else {
                mutex_.unlock();
                //                std::cerr << "Failed to converge for " << XT::Common::to_string(u, 15) << " with
                //                density "
                //                          << XT::Common::to_string(density, 15) << " at position "
                //                          << XT::Common::to_string(entity().geometry().center(), 15)
                //                          << " due to errors in the Cholesky decomposition!" << std::endl;
                DUNE_THROW(Dune::MathError, "Failure to converge!");
              }
            }

            const auto& alpha_tilde = alpha_k;
            const auto u_alpha_tilde = g_k + v;
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            if (!(density_tilde > 0.) || std::isinf(density_tilde))
              break;
            const auto alpha_prime = alpha_tilde - alpha_iso * std::log(density_tilde);
            const auto u_alpha_prime = calculate_u(alpha_prime);
            auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
            // checking realizability is cheap so we do not need the second stopping criterion
            if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)) {
              ret->first = alpha_prime + alpha_iso * std::log(density);
              ret->second = r;
              cache_.insert(v, alpha_prime);
              alpha_storage_[x_local] = ret->first;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              // backtracking line search
              while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
                // while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.) {
                // calculate alpha_new = alpha_k + zeta_k d_k
                auto alpha_new = d_k;
                alpha_new *= zeta_k;
                alpha_new += alpha_k;
                // calculate f(alpha_new)
                RangeFieldType f_new = calculate_f(alpha_new, v);
                if (pure_newton >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  alpha_k = alpha_new;
                  f_k = f_new;
                  pure_newton = 0.;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
              if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)
        mutex_.unlock();
        //        std::cerr << "Failed to converge for " << XT::Common::to_string(u, 15) << " with density "
        //                  << XT::Common::to_string(density, 15) << " at position "
        //                  << XT::Common::to_string(entity().geometry().center(), 15) << " due to too many iterations!"
        //                  << std::endl;
        DUNE_THROW(MathError, "Failed to converge");
      } // else ( value has not been calculated before )

    outside_all_loops:
      mutex_.unlock();
      return ret;
    } // ... get_alpha(...)

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, true)->first;

      std::fill(ret.begin(), ret.end(), 0.);
      // calculate < \mu m G_\alpha(u) >
      for (size_t nn = 0; nn < dimRange; ++nn) {
        if (nn > 0) {
          if (std::abs(alpha[nn] - alpha[nn - 1]) > taylor_tol_) {
            ret[nn] += 2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                           * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
                       + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                             * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1])
                             * std::exp(alpha[nn]);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha[nn] - alpha[nn - 1];
            size_t ll = 0;
            auto pow_frac = 1. / 6.;
            while (ll <= max_taylor_order_ - 3 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * ((ll * ll + 3 * ll + 2) * v_points_[nn] + (ll + 1) * v_points_[nn - 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 3);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn - 1]);
          }
        }
        if (nn < dimRange - 1) {
          if (std::abs(alpha[nn + 1] - alpha[nn]) > taylor_tol_) {
            ret[nn] += -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                           * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
                       + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                             * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn])
                             * std::exp(alpha[nn]);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha[nn + 1] - alpha[nn];
            size_t ll = 0;
            auto pow_frac = 1. / 6.;
            while (ll < 3 || (ll <= max_taylor_order_ - 3 && XT::Common::FloatCmp::ne(update, 0.))) {
              update = pow_frac * (2 * v_points_[nn] + (ll + 1) * v_points_[nn + 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 3);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
          }
        } // if (nn < dimRange - 1)
      } // nn
    } // void evaluate(...)

    virtual void evaluate_col(const size_t DXTC_DEBUG_ONLY(col),
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      evaluate(x_local, u, ret, param);
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& /*u*/,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      RangeType H_diag, J_diag;
      FieldVector<RangeFieldType, dimRange - 1> H_subdiag, J_subdiag;
      calculate_hessian(alpha, H_diag, H_subdiag);
      calculate_J(alpha, J_diag, J_subdiag);
      calculate_J_Hinv(ret, J_diag, J_subdiag, H_diag, H_subdiag);
    }

    virtual void partial_u_col(const size_t DXTC_DEBUG_ONLY(col),
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      partial_u(x_local, u, ret, param);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    RangeFieldType calculate_f(const RangeType& alpha_k, const RangeType& v) const
    {
      RangeFieldType ret(0);
      for (size_t ii = 0; ii < dimRange - 1; ++ii) {
        if (std::abs(alpha_k[ii + 1] - alpha_k[ii]) > taylor_tol_) {
          ret += (v_points_[ii + 1] - v_points_[ii]) / (alpha_k[ii + 1] - alpha_k[ii])
                 * (std::exp(alpha_k[ii + 1]) - std::exp(alpha_k[ii]));
        } else {
          RangeFieldType update = 1.;
          RangeFieldType result = 0.;
          size_t ll = 1;
          RangeFieldType base = alpha_k[ii + 1] - alpha_k[ii];
          auto pow_frac = 1.;
          while (ll <= max_taylor_order_ && XT::Common::FloatCmp::ne(update, 0.)) {
            update = pow_frac;
            result += update;
            ++ll;
            pow_frac *= base / ll;
          }
          assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
          ret += result * (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_k[ii]);
        }
      } // ii
      ret -= alpha_k * v;
      return ret;
    } // .. calculate_f(...)

    RangeType calculate_u(const RangeType& alpha_k) const
    {
      RangeType u(0);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        if (nn > 0) {
          if (std::abs(alpha_k[nn] - alpha_k[nn - 1]) > taylor_tol_) {
            u[nn] += -(v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
                         * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                     + (v_points_[nn] - v_points_[nn - 1]) / (alpha_k[nn] - alpha_k[nn - 1]) * std::exp(alpha_k[nn]);
          } else {
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_k[nn - 1] - alpha_k[nn];
            size_t ll = 0;
            RangeFieldType update = 1;
            RangeFieldType pow_frac = 0.5;
            while (ll <= max_taylor_order_ - 2 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac;
              result += update;
              ++ll;
              pow_frac *= base / (ll + 2);
            }
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            u[nn] += result * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
          }
        } // if (nn > 0)
        if (nn < dimRange - 1) {
          if (std::abs(alpha_k[nn + 1] - alpha_k[nn]) > taylor_tol_) {
            u[nn] += (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                         * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                     - (v_points_[nn + 1] - v_points_[nn]) / (alpha_k[nn + 1] - alpha_k[nn]) * std::exp(alpha_k[nn]);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            size_t ll = 0;
            RangeFieldType base = alpha_k[nn + 1] - alpha_k[nn];
            auto pow_frac = 0.5;
            while (ll <= max_taylor_order_ - 2 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac;
              result += update;
              ++ll;
              pow_frac *= base / (ll + 2);
            }
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            u[nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
          }
        } // if (nn < dimRange-1)
      } // nn
      return u;
    } // RangeType calculate_u(...)

    RangeType calculate_gradient(const RangeType& alpha_k, const RangeType& v) const
    {
      return calculate_u(alpha_k) - v;
    }

    void calculate_hessian(const RangeType& alpha_k,
                           RangeType& diag,
                           FieldVector<RangeFieldType, dimRange - 1>& subdiag) const
    {
      std::fill(diag.begin(), diag.end(), 0.);
      std::fill(subdiag.begin(), subdiag.end(), 0.);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        if (nn > 0) {
          if (std::abs(alpha_k[nn] - alpha_k[nn - 1]) > taylor_tol_) {
            subdiag[nn - 1] =
                (v_points_[nn] - v_points_[nn - 1])
                * ((std::exp(alpha_k[nn]) + std::exp(alpha_k[nn - 1])) / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
                   - 2. * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                         / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3));
            diag[nn] = (v_points_[nn] - v_points_[nn - 1])
                       * ((-2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2) + 1. / (alpha_k[nn] - alpha_k[nn - 1]))
                              * std::exp(alpha_k[nn])
                          + 2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3)
                                * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1])));

          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_k[nn - 1] - alpha_k[nn];
            RangeFieldType factor = (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
            size_t ll = 2;
            auto pow_frac = 1. / 6.;
            while (ll <= max_taylor_order_ - 1 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * (ll - 1.);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 1);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            subdiag[nn - 1] += result * factor;

            result = 0.;
            update = 1;
            ll = 3;
            pow_frac = 2. / 6.;
            while (ll <= max_taylor_order_ && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac;
              result += update;
              ++ll;
              pow_frac *= base / ll;
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            diag[nn] += result * factor;
          }
        } // if (nn > 0)
        if (nn < dimRange - 1) {
          if (std::abs(alpha_k[nn + 1] - alpha_k[nn]) > taylor_tol_) {
            diag[nn] += (v_points_[nn + 1] - v_points_[nn])
                        * ((-2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2) - 1. / (alpha_k[nn + 1] - alpha_k[nn]))
                               * std::exp(alpha_k[nn])
                           + 2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 3)
                                 * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn])));
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_k[nn + 1] - alpha_k[nn];
            size_t ll = 3;
            auto pow_frac = 2. / 6.;
            while (ll <= max_taylor_order_ && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac;
              result += update;
              ++ll;
              pow_frac *= base / ll;
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            diag[nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
          }
        } // if (nn < dimRange - 1)
      } // nn
    } // void calculate_hessian(...)

    void
    calculate_J(const RangeType& alpha_k, RangeType& diag, FieldVector<RangeFieldType, dimRange - 1>& subdiag) const
    {
      std::fill(diag.begin(), diag.end(), 0.);
      std::fill(subdiag.begin(), subdiag.end(), 0.);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        if (nn > 0) {
          if (std::abs(alpha_k[nn] - alpha_k[nn - 1]) > taylor_tol_) {
            subdiag[nn - 1] =
                (v_points_[nn] - v_points_[nn - 1])
                    * ((v_points_[nn] * std::exp(alpha_k[nn]) + v_points_[nn - 1] * std::exp(alpha_k[nn - 1]))
                           / std::pow(alpha_k[nn - 1] - alpha_k[nn], 2)
                       + 2.
                             * ((2 * v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn])
                                - (2 * v_points_[nn - 1] - v_points_[nn]) * std::exp(alpha_k[nn - 1]))
                             / std::pow(alpha_k[nn - 1] - alpha_k[nn], 3))
                + 6. * std::pow(v_points_[nn] - v_points_[nn - 1], 2)
                      * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                      / std::pow(alpha_k[nn - 1] - alpha_k[nn], 4);
            diag[nn] = 6 * std::pow(v_points_[nn - 1] - v_points_[nn], 2)
                           * (std::exp(alpha_k[nn - 1]) - std::exp(alpha_k[nn]))
                           / std::pow(alpha_k[nn - 1] - alpha_k[nn], 4)
                       + 2. * (v_points_[nn] - v_points_[nn - 1])
                             * (v_points_[nn - 1] * std::exp(alpha_k[nn - 1])
                                - (3 * v_points_[nn] - 2 * v_points_[nn - 1]) * std::exp(alpha_k[nn]))
                             / std::pow(alpha_k[nn - 1] - alpha_k[nn], 3)
                       - v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn])
                             / (alpha_k[nn - 1] - alpha_k[nn])
                       - (std::pow(v_points_[nn - 1], 2) - 4 * v_points_[nn] * v_points_[nn - 1]
                          + 3. * std::pow(v_points_[nn], 2))
                             * std::exp(alpha_k[nn]) / std::pow(alpha_k[nn - 1] - alpha_k[nn], 2);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_k[nn - 1] - alpha_k[nn];
            RangeFieldType factor = (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
            size_t ll = 0;
            auto pow_frac = 1. / 24.;
            while (ll < 2 || (ll <= max_taylor_order_ - 4 && XT::Common::FloatCmp::ne(update, 0.))) {
              update = pow_frac * ((ll * ll + 3 * ll + 2) * v_points_[nn - 1] + (2 * ll + 2) * v_points_[nn]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 4);
            } // ll
            subdiag[nn - 1] += result * factor;
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));

            result = 0.;
            update = 1;
            ll = 0;
            pow_frac = 1. / 24.;
            while (ll < 4 || (ll <= max_taylor_order_ - 4 && XT::Common::FloatCmp::ne(update, 0.))) {
              update = pow_frac * (6 * v_points_[nn] + (2 * ll + 2) * v_points_[nn - 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 4);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            diag[nn] += result * factor;
          }
        } // if (nn > 0)
        if (nn < dimRange - 1) {
          if (std::abs(alpha_k[nn + 1] - alpha_k[nn]) > taylor_tol_) {
            diag[nn] += 6 * std::pow(v_points_[nn] - v_points_[nn + 1], 2)
                            * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn + 1]))
                            / std::pow(alpha_k[nn] - alpha_k[nn + 1], 4)
                        + 2. * (v_points_[nn] - v_points_[nn + 1])
                              * (v_points_[nn + 1] * std::exp(alpha_k[nn + 1])
                                 - (3 * v_points_[nn] - 2 * v_points_[nn + 1]) * std::exp(alpha_k[nn]))
                              / std::pow(alpha_k[nn] - alpha_k[nn + 1], 3)
                        - v_points_[nn] * (v_points_[nn] - v_points_[nn + 1]) * std::exp(alpha_k[nn])
                              / (alpha_k[nn] - alpha_k[nn + 1])
                        + (std::pow(v_points_[nn + 1], 2) - 4 * v_points_[nn] * v_points_[nn + 1]
                           + 3. * std::pow(v_points_[nn], 2))
                              * std::exp(alpha_k[nn]) / std::pow(alpha_k[nn] - alpha_k[nn + 1], 2);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_k[nn + 1] - alpha_k[nn];
            size_t ll = 0;
            auto pow_frac = 1. / 24.;
            while (ll < 4 || (ll <= max_taylor_order_ - 4 && XT::Common::FloatCmp::ne(update, 0.))) {
              update = pow_frac * (6 * v_points_[nn] + (2 * ll + 2) * v_points_[nn + 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 4);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            diag[nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
          }
        } // if (nn < dimRange - 1)
      } // nn
    } // void calculate_J(...)

    // calculates ret = J H^{-1}. Both J and H are symmetric tridiagonal, H is positive definite.
    static void calculate_J_Hinv(MatrixType& ret,
                                 const RangeType& J_diag,
                                 const FieldVector<RangeFieldType, dimRange - 1>& J_subdiag,
                                 RangeType& H_diag,
                                 FieldVector<RangeFieldType, dimRange - 1>& H_subdiag)
    {
      // factorize H = LDL^T, where L is unit lower bidiagonal and D is diagonal
      // H_diag is overwritten by the diagonal elements of D
      // H_subdiag is overwritten by the subdiagonal elements of L
      XT::LA::tridiagonal_ldlt(H_diag, H_subdiag);

      // copy J to dense matrix
      std::fill(ret.begin(), ret.end(), 0.);
      for (size_t ii = 0; ii < dimRange - 1; ++ii) {
        ret[ii][ii] = J_diag[ii];
        ret[ii + 1][ii] = J_subdiag[ii];
        ret[ii][ii + 1] = J_subdiag[ii];
      }
      ret[dimRange - 1][dimRange - 1] = J_diag[dimRange - 1];

      // Solve ret H = J which is equivalent to (as H and J are symmetric) to H ret^T = J;
      XT::LA::solve_tridiagonal_ldlt_factorized(H_diag, H_subdiag, ret);
      // transpose ret
      for (size_t ii = 0; ii < dimRange; ++ii)
        for (size_t jj = 0; jj < ii; ++jj)
          std::swap(ret[jj][ii], ret[ii][jj]);
    } // void calculate_J_Hinv(...)

    const BasisfunctionType& basis_functions_;
    const std::vector<RangeFieldType>& v_points_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const RangeFieldType taylor_tol_;
    const size_t max_taylor_order_;
    const std::string name_;
    LocalCacheType& cache_;
    AlphaStorageType& alpha_storage_;
    std::mutex& mutex_;
  }; // class Localfunction>

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           v_points_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           taylor_tol_,
                                           max_taylor_order_,
                                           cache_[index],
                                           alpha_storage_[index],
                                           mutexes_[index]);
  }

  // calculate \sum_{i=1}^d < v_i_+ m \psi >_+ n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& /*u_i*/,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t DXTC_DEBUG_ONLY(dd),
                                       const XT::Common::Parameter& /*param*/,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(dd == 0);
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_stored_alpha(x_local_entity);
    StateRangeType alpha_j;
    const bool boundary = bool(param_neighbor.get("boundary")[0]);
    if (boundary)
      alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, true)->first;
    else
      alpha_j = local_function_neighbor->get_stored_alpha(x_local_neighbor);
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (dimRange % 2 || nn != dimRange / 2) {
          const auto& alpha = (n_ij[0] * (v_points_[nn - 1] + v_points_[nn]) / 2. > 0.) ? alpha_i : alpha_j;
          if (std::abs(alpha[nn] - alpha[nn - 1]) > taylor_tol_) {
            ret[nn] += 2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                           * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
                       + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                             * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1])
                             * std::exp(alpha[nn]);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha[nn - 1] - alpha[nn];
            size_t ll = 0;
            auto pow_frac = 1. / 6.;
            while (ll < 3 || (ll <= max_taylor_order_ - 3 && XT::Common::FloatCmp::ne(update, 0.))) {
              update = pow_frac * (2 * v_points_[nn] + (ll + 1) * v_points_[nn - 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 3);
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn]);
          }
        } else { //  if (dimRange % 2 || nn != dimRange/2)
          const auto& alpha_pos = n_ij[0] > 0. ? alpha_i : alpha_j;
          const auto& alpha_neg = n_ij[0] > 0. ? alpha_j : alpha_i;
          if (std::abs(alpha_neg[nn] - alpha_neg[nn - 1]) > taylor_tol_) {
            ret[nn] += -2. * std::pow(v_points_[nn], 2)
                       * (4. / std::pow(alpha_neg[nn - 1] - alpha_neg[nn], 3)
                              * (std::exp((alpha_neg[nn] + alpha_neg[nn - 1]) / 2.) - std::exp(alpha_neg[nn - 1]))
                          + 1. / std::pow(alpha_neg[nn - 1] - alpha_neg[nn], 2)
                                * (std::exp((alpha_neg[nn] + alpha_neg[nn - 1]) / 2.) + std::exp(alpha_neg[nn - 1])));

          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_neg[nn] - alpha_neg[nn - 1];
            size_t ll = 2;
            auto pow_frac = 1. / 24.;
            while (ll <= max_taylor_order_ - 1 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * (ll - 1.);
              result += update;
              ++ll;
              pow_frac *= base / (2. * (ll + 1));
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * -2. * std::pow(v_points_[nn], 2) * std::exp(alpha_neg[nn - 1]);
          }
          if (std::abs(alpha_pos[nn] - alpha_pos[nn - 1]) > taylor_tol_) {
            ret[nn] += 2. * std::pow(v_points_[nn], 2)
                       * (4. / std::pow(alpha_pos[nn - 1] - alpha_pos[nn], 3)
                              * (std::exp((alpha_pos[nn] + alpha_pos[nn - 1]) / 2.) - std::exp(alpha_pos[nn]))
                          + 1. / std::pow(alpha_pos[nn - 1] - alpha_pos[nn], 2)
                                * (std::exp((alpha_pos[nn] + alpha_pos[nn - 1]) / 2.) - 3. * std::exp(alpha_pos[nn]))
                          - 1. / (alpha_pos[nn - 1] - alpha_pos[nn]) * std::exp(alpha_pos[nn]));
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_pos[nn - 1] - alpha_pos[nn];
            auto pow_frac = 1. / 24.;
            size_t ll = 2;
            while (ll <= max_taylor_order_ - 1 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * (ll + 3);
              result += update;
              ++ll;
              pow_frac *= base / (2. * (ll + 1));
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * 2. * std::pow(v_points_[nn], 2) * std::exp(alpha_pos[nn]);
          } // else (alpha_n - alpha_{n-1} != 0)
        } // else (dimRange % 2 || nn != dimRange/2)
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        if (dimRange % 2 || nn != dimRange / 2 - 1) {
          const auto& alpha = (n_ij[0] * (v_points_[nn] + v_points_[nn + 1]) / 2. > 0.) ? alpha_i : alpha_j;
          if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], 0., taylor_tol_)) {
            ret[nn] += -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                           * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
                       + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                             * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn])
                             * std::exp(alpha[nn]);
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha[nn + 1] - alpha[nn];
            size_t ll = 0;
            auto pow_frac = 1. / 6.;
            while (ll < 3
                   || (ll <= max_taylor_order_ - 3 && XT::Common::FloatCmp::ne(result, result + update, 1e-16, 0.))) {
              update = pow_frac * (2 * v_points_[nn] + (ll + 1) * v_points_[nn + 1]);
              result += update;
              ++ll;
              pow_frac *= base / (ll + 3);
            } // ll
            ret[nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
          }
        } else { // if (dimRange % 2 || nn != dimRange / 2 - 1)
          const auto& alpha_pos = n_ij[0] > 0. ? alpha_i : alpha_j;
          const auto& alpha_neg = n_ij[0] > 0. ? alpha_j : alpha_i;
          if (std::abs(alpha_neg[nn + 1] - alpha_neg[nn]) > taylor_tol_) {
            ret[nn] += -2. * std::pow(v_points_[nn + 1], 2)
                       * (-4. / std::pow(alpha_neg[nn + 1] - alpha_neg[nn], 3)
                              * (std::exp(alpha_neg[nn]) - std::exp((alpha_neg[nn + 1] + alpha_neg[nn]) / 2.))
                          - 1. / std::pow(alpha_neg[nn + 1] - alpha_neg[nn], 2)
                                * (3 * std::exp(alpha_neg[nn]) - std::exp((alpha_neg[nn + 1] + alpha_neg[nn]) / 2.))
                          - 1. / (alpha_neg[nn + 1] - alpha_neg[nn]) * std::exp(alpha_neg[nn]));
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_neg[nn + 1] - alpha_neg[nn];
            auto pow_frac = 1. / 24.;
            size_t ll = 2;
            while (ll <= max_taylor_order_ - 1 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * (ll + 3);
              result += update;
              ++ll;
              pow_frac *= base / (2. * (ll + 1));
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * -2. * std::pow(v_points_[nn + 1], 2) * std::exp(alpha_neg[nn]);
          }
          if (std::abs(alpha_pos[nn + 1] - alpha_pos[nn]) > taylor_tol_) {
            ret[nn] += 2. * std::pow(v_points_[nn + 1], 2)
                       * (4. / std::pow(alpha_pos[nn + 1] - alpha_pos[nn], 3)
                              * (std::exp((alpha_pos[nn + 1] + alpha_pos[nn]) / 2.) - std::exp(alpha_pos[nn + 1]))
                          + 1. / std::pow(alpha_pos[nn + 1] - alpha_pos[nn], 2)
                                * (std::exp((alpha_pos[nn + 1] + alpha_pos[nn]) / 2.) + std::exp(alpha_pos[nn + 1])));
          } else {
            RangeFieldType update = 1.;
            RangeFieldType result = 0.;
            RangeFieldType base = alpha_pos[nn] - alpha_pos[nn + 1];
            auto pow_frac = 1. / 24.;
            size_t ll = 2;
            while (ll <= max_taylor_order_ - 1 && XT::Common::FloatCmp::ne(update, 0.)) {
              update = pow_frac * (ll - 1.);
              result += update;
              ++ll;
              pow_frac *= base / (2. * (ll + 1));
            } // ll
            assert(!(std::isinf(pow_frac) || std::isnan(pow_frac)));
            ret[nn] += result * 2. * std::pow(v_points_[nn + 1], 2) * std::exp(alpha_pos[nn + 1]);
          } // else (alpha_n - alpha_{n-1} != 0)
        } // else (dimRange % 2 || nn != dimRange / 2 - 1)
      } // if (nn < dimRange - 1)
    } // nn
    ret *= n_ij[0];
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

  const BasisfunctionType& basis_functions() const
  {
    return basis_functions_;
  }

private:
  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  const std::vector<RangeFieldType>& v_points_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const RangeFieldType taylor_tol_;
  const size_t max_taylor_order_;
  const std::string name_;
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see constructor)
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<AlphaStorageType> alpha_storage_;
  mutable std::vector<std::mutex> mutexes_;
};
#endif


#if 1
/**
 * Specialization of EntropyBasedLocalFlux for 1D Hatfunctions (no change of basis)
 */
template <class GridLayerImp, class U>
class EntropyBasedLocalFlux<
    HatFunctionMomentBasis<typename U::DomainFieldType, 1, typename U::RangeFieldType, U::dimRange, 1, 1>,
    GridLayerImp,
    U>
  : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                           typename U::DomainFieldType,
                                                           GridLayerImp::dimension,
                                                           U,
                                                           0,
                                                           typename U::RangeFieldType,
                                                           U::dimRange,
                                                           1>
{
  using BaseType =
      typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                               typename U::DomainFieldType,
                                                               GridLayerImp::dimension,
                                                               U,
                                                               0,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               1>;
  using ThisType = EntropyBasedLocalFlux;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using typename BaseType::StateType;
  using BasisfunctionType = HatFunctionMomentBasis<DomainFieldType, 1, RangeFieldType, dimRange, 1, 1>;
  using GridLayerType = GridLayerImp;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, 1>;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using AlphaReturnType = typename std::pair<StateRangeType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, StateRangeType>;
  using AlphaStorageType = std::map<DomainType, StateRangeType, XT::Common::VectorFloatLess>;
  static const size_t cache_size = 4 * dimDomain + 2;
  static const size_t num_intervals = dimRange - 1;
  static const size_t block_size = 2;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;
  using BasisValuesMatrixType = FieldVector<std::vector<LocalVectorType>, num_intervals>;
  using QuadraturePointsType = FieldVector<std::vector<RangeFieldType>, num_intervals>;
  using QuadratureWeightsType = FieldVector<std::vector<RangeFieldType>, num_intervals>;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , grid_points_(basis_functions_.triangulation())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , alpha_storage_(index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {
    const auto& quadratures = basis_functions_.quadratures();
    assert(quadratures.size() == grid_points_.size() - 1);
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (const auto& quad_point : quadratures[jj]) {
        quad_points_[jj].emplace_back(quad_point.position()[0]);
        quad_weights_[jj].emplace_back(quad_point.weight());
      }
    } // jj
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      M_[jj].resize(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        M_[jj][ll] = basis_functions_.evaluate_on_interval(quad_points_[jj][ll], jj);
    } // jj
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColPartialURangeType;
    using typename LocalfunctionType::ColRangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const QuadraturePointsType& quad_points,
                  const QuadratureWeightsType& quad_weights,
                  const std::vector<RangeFieldType>& grid_points,
                  const BasisValuesMatrixType& M,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  LocalCacheType& cache,
                  AlphaStorageType& alpha_storage,
                  std::mutex& mutex)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , quad_points_(quad_points)
      , quad_weights_(quad_weights)
      , grid_points_(grid_points)
      , M_(M)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , cache_(cache)
      , alpha_storage_(alpha_storage)
      , mutex_(mutex)
    {}

    using LocalfunctionType::entity;

    static bool is_realizable(const RangeType& u)
    {
      for (const auto& u_i : u)
        if (!(u_i > 0.) || std::isinf(u_i))
          return false;
      return true;
    }

    void store_alpha(const DomainType& x_local, const StateRangeType& alpha)
    {
      alpha_storage_[x_local] = alpha;
    }

    StateRangeType get_stored_alpha(const DomainType& x_local) const
    {
      return alpha_storage_.at(x_local);
    }

    // temporary vectors to store inner products and exponentials
    FieldVector<std::vector<RangeFieldType>, num_intervals>& working_storage() const
    {
      thread_local FieldVector<std::vector<RangeFieldType>, num_intervals> work_vec;
      for (size_t jj = 0; jj < num_intervals; ++jj)
        work_vec[jj].resize(quad_points_[jj].size());
      return work_vec;
    }

    template <class GridLayerType>
    void center_results_to_intersections(const GridLayerType& grid_layer)
    {
      const auto center = entity().geometry().local(entity().geometry().center());
      const auto center_alpha = get_stored_alpha(center);
      for (const auto& intersection : Dune::intersections(grid_layer, entity()))
        store_alpha(entity().geometry().local(intersection.geometry().center()), center_alpha);
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& x_local,
                                               const StateRangeType& u,
                                               const XT::Common::Parameter& param,
                                               const bool regularize) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      auto ret = std::make_unique<AlphaReturnType>();
      mutex_.lock();
      if (boundary)
        cache_.set_capacity(cache_size + dimDomain);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      if (!(density > 0.) || std::isinf(density)) {
        mutex_.unlock();
        DUNE_THROW(Dune::MathError, "Negative density!");
      }
      RangeType u_prime = u / density;
      RangeType alpha_iso = basis_functions_.alpha_iso();

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u_prime, 1e-14, 1e-14)) {
        const auto alpha_prime = cache_iterator->second;
        ret->first = alpha_prime + alpha_iso * std::log(density);
        ret->second = 0.;
        alpha_storage_[x_local] = ret->first;
        mutex_.unlock();
        return ret;
      } else {
        RangeFieldType tau_prime = std::min(
            tau_ / ((1 + std::sqrt(dimRange) * u_prime.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
        // The hessian H is always symmetric and tridiagonal, so we only need to store the diagonal and subdiagonal
        // elements
        RangeType H_diag;
        FieldVector<RangeFieldType, dimRange - 1> H_subdiag;

        // calculate moment vector for isotropic distribution
        RangeType u_iso = basis_functions_.u_iso();
        RangeType v;
        RangeType alpha_k = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;
        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence_) {
          // regularize u
          v = u_prime;
          if (r > 0) {
            alpha_k = alpha_iso;
            RangeType r_times_u_iso(u_iso);
            r_times_u_iso *= r;
            v *= 1 - r;
            v += r_times_u_iso;
          }

          // calculate f_0
          RangeFieldType f_k = calculate_f(alpha_k, v);

          int pure_newton = 0;
          StateRangeType g_k, d_k, minus_g_k, u_alpha_prime;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used
            if (kk > k_0_ && r < r_max)
              break;
            // calculate gradient g
            calculate_gradient(alpha_k, v, g_k);
            // calculate Hessian H
            calculate_hessian(alpha_k, M_, H_diag, H_subdiag);
            // calculate descent direction d_k;
            minus_g_k = g_k;
            minus_g_k *= -1;
            try {
              d_k = minus_g_k;
              XT::LA::solve_sym_tridiag_posdef(H_diag, H_subdiag, d_k);
            } catch (const Dune::MathError&) {
              if (r < r_max) {
                break;
              } else {
                mutex_.unlock();
                DUNE_THROW(Dune::MathError, "Failure to converge!");
              }
            }

            const auto& alpha_tilde = alpha_k;
            const auto u_alpha_tilde = g_k + v;
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            if (!(density_tilde > 0.) || std::isinf(density_tilde))
              break;
            const auto alpha_prime = alpha_tilde - alpha_iso * std::log(density_tilde);
            calculate_u(alpha_prime, u_alpha_prime);
            auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
            // checking realizability is cheap so we do not need the second stopping criterion
            if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)) {
              ret->first = alpha_prime + alpha_iso * std::log(density);
              ret->second = r;
              cache_.insert(v, alpha_prime);
              alpha_storage_[x_local] = ret->first;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              // backtracking line search
              while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
                // while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.) {
                // calculate alpha_new = alpha_k + zeta_k d_k
                auto alpha_new = d_k;
                alpha_new *= zeta_k;
                alpha_new += alpha_k;
                // calculate f(alpha_new)
                RangeFieldType f_new = calculate_f(alpha_new, v);
                if (pure_newton >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  alpha_k = alpha_new;
                  f_k = f_new;
                  pure_newton = 0.;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
              if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)
        mutex_.unlock();
        DUNE_THROW(MathError, "Failed to converge");
      } // else ( value has not been calculated before )

    outside_all_loops:
      mutex_.unlock();
      return ret;
    } // ... get_alpha(...)

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      std::fill(ret.begin(), ret.end(), 0.);
      const auto alpha = get_alpha(x_local, u, param, true)->first;
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      LocalVectorType local_alpha;
      for (size_t jj = 0; jj < num_intervals; ++jj) {
        for (size_t ii = 0; ii < 2; ++ii)
          local_alpha[ii] = alpha[jj + ii];
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          auto factor_ll = std::exp(local_alpha * basis_ll) * quad_points_[jj][ll] * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 2; ++ii)
            ret[jj + ii] += basis_ll[ii] * factor_ll;
        } // ll (quad points)
      } // jj (intervals)
    } // void evaluate(...)

    virtual void evaluate_col(const size_t DXTC_DEBUG_ONLY(col),
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      evaluate(x_local, u, ret, param);
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& /*u*/,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& /*param*/) const override
    {
      const auto alpha = get_stored_alpha(x_local);
      RangeType H_diag, J_diag;
      FieldVector<RangeFieldType, dimRange - 1> H_subdiag, J_subdiag;
      calculate_hessian(alpha, M_, H_diag, H_subdiag);
      calculate_J(M_, J_diag, J_subdiag);
      calculate_J_Hinv(ret, J_diag, J_subdiag, H_diag, H_subdiag);
    }

    virtual void partial_u_col(const size_t DXTC_DEBUG_ONLY(col),
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      partial_u(x_local, u, ret, param);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    RangeFieldType calculate_f(const StateRangeType& alpha, const StateRangeType& v) const
    {
      RangeFieldType ret(0.);
      XT::Common::FieldVector<RangeFieldType, block_size> local_alpha;
      for (size_t jj = 0; jj < num_intervals; ++jj) {
        for (size_t ii = 0; ii < 2; ++ii)
          local_alpha[ii] = alpha[jj + ii];
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
          ret += std::exp(local_alpha * M_[jj][ll]) * quad_weights_[jj][ll];
      } // jj (intervals)
      ret -= alpha * v;
      return ret;
    } // void calculate_u(...)

    void calculate_u(const StateRangeType& alpha, StateRangeType& u) const
    {
      std::fill(u.begin(), u.end(), 0.);
      LocalVectorType local_alpha;
      for (size_t jj = 0; jj < num_intervals; ++jj) {
        for (size_t ii = 0; ii < 2; ++ii)
          local_alpha[ii] = alpha[jj + ii];
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          auto factor_ll = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 2; ++ii)
            u[jj + ii] += basis_ll[ii] * factor_ll;
        } // ll (quad points)
      } // jj (intervals)
    } // void calculate_u(...)

    void calculate_gradient(const StateRangeType& alpha, const StateRangeType& v, StateRangeType& g_k) const
    {
      calculate_u(alpha, g_k);
      g_k -= v;
    }

    void calculate_hessian(const StateRangeType& alpha,
                           const BasisValuesMatrixType& M,
                           StateRangeType& H_diag,
                           FieldVector<RangeFieldType, dimRange - 1>& H_subdiag) const
    {
      std::fill(H_diag.begin(), H_diag.end(), 0.);
      std::fill(H_subdiag.begin(), H_subdiag.end(), 0.);
      LocalVectorType local_alpha;
      auto& work_vecs = working_storage();
      for (size_t jj = 0; jj < num_intervals; ++jj) {
        for (size_t ii = 0; ii < 2; ++ii)
          local_alpha[ii] = alpha[jj + ii];
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M[jj][ll];
          work_vecs[jj][ll] = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 2; ++ii)
            H_diag[jj + ii] += std::pow(basis_ll[ii], 2) * work_vecs[jj][ll];
          H_subdiag[jj] += basis_ll[0] * basis_ll[1] * work_vecs[jj][ll];
        } // ll (quad points)
      } // jj (intervals)
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) values
    void calculate_J(const BasisValuesMatrixType& M,
                     StateRangeType& J_diag,
                     FieldVector<RangeFieldType, dimRange - 1>& J_subdiag) const
    {
      std::fill(J_diag.begin(), J_diag.end(), 0.);
      std::fill(J_subdiag.begin(), J_subdiag.end(), 0.);
      const auto& work_vecs = working_storage();
      for (size_t jj = 0; jj < num_intervals; ++jj) {
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M[jj][ll];
          for (size_t ii = 0; ii < 2; ++ii)
            J_diag[jj + ii] += std::pow(basis_ll[ii], 2) * work_vecs[jj][ll] * quad_points_[jj][ll];
          J_subdiag[jj] += basis_ll[0] * basis_ll[1] * work_vecs[jj][ll] * quad_points_[jj][ll];
        } // ll (quad points)
      } // jj (intervals)
    } // void calculate_J(...)

    // calculates ret = J H^{-1}. Both J and H are symmetric tridiagonal, H is positive definite.
    static void calculate_J_Hinv(MatrixType& ret,
                                 const StateRangeType& J_diag,
                                 const FieldVector<RangeFieldType, dimRange - 1>& J_subdiag,
                                 StateRangeType& H_diag,
                                 FieldVector<RangeFieldType, dimRange - 1>& H_subdiag)
    {
      // factorize H = LDL^T, where L is unit lower bidiagonal and D is diagonal
      // H_diag is overwritten by the diagonal elements of D
      // H_subdiag is overwritten by the subdiagonal elements of L
      XT::LA::tridiagonal_ldlt(H_diag, H_subdiag);

      // copy J to dense matrix
      std::fill(ret.begin(), ret.end(), 0.);
      for (size_t ii = 0; ii < dimRange - 1; ++ii) {
        ret[ii][ii] = J_diag[ii];
        ret[ii + 1][ii] = J_subdiag[ii];
        ret[ii][ii + 1] = J_subdiag[ii];
      }
      ret[dimRange - 1][dimRange - 1] = J_diag[dimRange - 1];

      // Solve ret H = J which is equivalent to (as H and J are symmetric) to H ret^T = J;
      XT::LA::solve_tridiagonal_ldlt_factorized(H_diag, H_subdiag, ret);
      // transpose ret
      for (size_t ii = 0; ii < dimRange; ++ii)
        for (size_t jj = 0; jj < ii; ++jj)
          std::swap(ret[jj][ii], ret[ii][jj]);
    } // void calculate_J_Hinv(...)

    const BasisfunctionType& basis_functions_;
    const QuadraturePointsType& quad_points_;
    const QuadratureWeightsType& quad_weights_;
    const std::vector<RangeFieldType>& grid_points_;
    const BasisValuesMatrixType& M_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const std::string name_;
    LocalCacheType& cache_;
    AlphaStorageType& alpha_storage_;
    std::mutex& mutex_;
  }; // class Localfunction>

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           quad_points_,
                                           quad_weights_,
                                           grid_points_,
                                           M_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           cache_[index],
                                           alpha_storage_[index],
                                           mutexes_[index]);
  }

  // calculate \sum_{i=1}^d < v_i_+ m \psi >_+ n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& /*u_i*/,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t DXTC_DEBUG_ONLY(dd),
                                       const XT::Common::Parameter& /*param*/,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(dd == 0);
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    const bool boundary = static_cast<bool>(param_neighbor.get("boundary")[0]);
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_stored_alpha(x_local_entity);
    StateRangeType alpha_j;
    if (boundary)
      alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, true)->first;
    else
      alpha_j = local_function_neighbor->get_stored_alpha(x_local_neighbor);
    StateRangeType ret(0);
    LocalVectorType local_alpha_i, local_alpha_j;
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ii = 0; ii < 2; ++ii) {
        local_alpha_i[ii] = alpha_i[jj + ii];
        local_alpha_j[ii] = alpha_j[jj + ii];
      }
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M_[jj][ll];
        const auto position = quad_points_[jj][ll];
        RangeFieldType factor =
            position * n_ij[0] > 0. ? std::exp(local_alpha_i * basis_ll) : std::exp(local_alpha_j * basis_ll);
        factor *= quad_weights_[jj][ll] * position;
        for (size_t ii = 0; ii < 2; ++ii)
          ret[jj + ii] += basis_ll[ii] * factor;
      } // ll (quad points)
    } // jj (intervals)
    ret *= n_ij[0];
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

  const BasisfunctionType& basis_functions() const
  {
    return basis_functions_;
  }

private:
  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  const std::vector<RangeFieldType>& grid_points_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const std::string name_;
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see constructor)
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<AlphaStorageType> alpha_storage_;
  mutable std::vector<std::mutex> mutexes_;
};
#endif


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
