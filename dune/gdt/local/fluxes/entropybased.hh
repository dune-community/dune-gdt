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
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/memory.hh>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/container/pattern.hh>

#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/hatfunctions.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/piecewise_monomials.hh>

#if HAVE_CLP
#include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

namespace Dune {
namespace GDT {


template <class StateRangeType, class VectorType>
class EntropyLocalCache
{
public:
  typedef typename std::map<StateRangeType, VectorType, XT::Common::FieldVectorLess> MapType;
  typedef typename MapType::iterator IteratorType;
  typedef typename MapType::const_iterator ConstIteratorType;
  using RangeFieldType = typename StateRangeType::value_type;

  EntropyLocalCache(const size_t capacity)
    : capacity_(capacity)
  {
  }

  void insert(const StateRangeType& u, const VectorType& alpha)
  {
    cache_.insert(std::make_pair(u, alpha));
    keys_.push_back(u);
    if (cache_.size() > capacity_) {
      cache_.erase(keys_.front());
      keys_.pop_front();
    }
  }

  void keep(const StateRangeType& u)
  {
    keys_.remove(u);
    keys_.push_back(u);
  }

  ConstIteratorType find_closest(const StateRangeType& u) const
  {
    ConstIteratorType ret = cache_.begin();
    if (ret == end())
      return ret;
    RangeFieldType distance = (u - ret->first).two_norm2();
    RangeFieldType new_distance = distance;
    auto it = ret;
    while (++it != end()) {
      if ((new_distance = (u - it->first).two_norm2()) < distance) {
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

  void increase_capacity(const size_t new_capacity)
  {
    if (new_capacity > capacity_)
      capacity_ = new_capacity;
  }

private:
  std::size_t capacity_;
  MapType cache_;
  std::list<StateRangeType> keys_;
};

template <class BasisfunctionImp, class GridLayerImp, class U>
class EntropyBasedLocalFlux;

#if HAVE_CLP
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
  typedef typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                                   typename BasisfunctionImp::DomainFieldType,
                                                                   BasisfunctionImp::dimFlux,
                                                                   U,
                                                                   0,
                                                                   typename BasisfunctionImp::RangeFieldType,
                                                                   BasisfunctionImp::dimRange,
                                                                   BasisfunctionImp::dimFlux>
      BaseType;
  typedef EntropyBasedLocalFlux ThisType;

public:
  typedef BasisfunctionImp BasisfunctionType;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using VectorType = FieldVector<RangeFieldType, dimRange>;
  using BasisValuesMatrixType = XT::LA::CommonDenseMatrix<RangeFieldType>;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureRuleType;
  typedef std::pair<VectorType, RangeFieldType> AlphaReturnType;
  typedef EntropyLocalCache<StateRangeType, VectorType> LocalCacheType;
  static const size_t cache_size = 2 * dimDomain + 2;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& quadrature,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = XT::LA::eye_matrix<MatrixType>(dimRange,
                                                                     XT::LA::dense_pattern(dimRange, dimRange)),
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , quad_points_(quadrature.size())
    , quad_weights_(quadrature.size())
    , M_(dimRange, quadrature.size(), 0., 0)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , T_minus_one_(T_minus_one)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , mutexes_(index_set_.size(0))
  {
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      quad_points_[ll] = quadrature[ll].position();
      quad_weights_[ll] = quadrature[ll].weight();
      const auto val = basis_functions_.evaluate(quad_points_[ll]);
      for (size_t ii = 0; ii < dimRange; ++ii)
        M_.set_entry(ii, ll, val[ii]);
    }
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using LocalfunctionType::dimRange;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

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
                  std::mutex& mutex,
                  XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp)
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
      , mutex_(mutex)
      , lp_(lp)
    {
    }

    void setup_linear_program() const
    {
      if (!*lp_) {
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
          row_indices[ii] = ii;

        // set columns for quadrature points
        for (int ii = 0; ii < num_cols; ++ii) {
          const auto v_i = basis_functions_.evaluate(quad_points_[ii]);
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

    //  RangeFieldType solve_linear_program(const RangeType& u_bar, const RangeType& u_l, const size_t index)
    bool is_realizable(const StateRangeType& u) const
    {
      const auto u_prime = u / basis_functions_.density(u);
      setup_linear_program();
      auto& lp = **lp_;
      constexpr int num_rows = static_cast<int>(dimRange);
      // set rhs (equality constraints, so set both bounds equal
      for (int ii = 0; ii < num_rows; ++ii) {
        lp.setRowLower(ii, u_prime[ii]);
        lp.setRowUpper(ii, u_prime[ii]);
      }
      // Now check solvability
      lp.primal();
      return lp.primalFeasible();
    }

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
      const size_t num_quad_points = quad_points_.size();
      const auto& M_backend = M.backend();
      std::fill(scalar_products.begin(), scalar_products.end(), 0.);
      for (size_t ii = 0; ii < dimRange; ++ii) {
        const auto factor = beta_in[ii];
        for (size_t ll = 0; ll < num_quad_points; ++ll)
          scalar_products[ll] += M_backend.get_entry_ref(ii, ll) * factor;
      }
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
                                   bool only_first_component = false) const
    {
      auto& work_vec = working_storage();
      calculate_scalar_products(beta_in, M2, work_vec);
      apply_exponential(work_vec);
      std::fill(ret.begin(), ret.end(), 0.);
      const size_t num_quad_points = quad_weights_.size();
      const auto& M1_backend = M1.backend();
      for (size_t ii = 0; ii < (only_first_component ? 1 : dimRange); ++ii) {
        RangeFieldType result(0.);
        for (size_t ll = 0; ll < num_quad_points; ++ll)
          result += M1_backend.get_entry_ref(ii, ll) * work_vec[ll] * quad_weights_[ll];
        ret[ii] = result;
      } // ii
    }

    void apply_inverse_matrix(const MatrixType& T_k, BasisValuesMatrixType& M) const
    {
      assert(quad_points_.size() < std::numeric_limits<int>::max());
      XT::Common::Blas::dtrsm(XT::Common::Blas::row_major(),
                              XT::Common::Blas::left(),
                              XT::Common::Blas::lower(),
                              XT::Common::Blas::no_trans(),
                              XT::Common::Blas::non_unit(),
                              dimRange,
                              static_cast<int>(quad_points_.size()),
                              1.,
                              &(T_k[0][0]),
                              dimRange,
                              M.data(),
                              static_cast<int>(quad_points_.size()));
    }

    AlphaReturnType get_alpha(const DomainType& /*x_local*/,
                              const StateRangeType& u,
                              const XT::Common::Parameter& param,
                              const bool regularize,
                              const bool only_cache) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      AlphaReturnType ret;
      mutex_.lock();
      if (boundary)
        cache_.increase_capacity(2 * cache_size);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      RangeType u_prime = u / density;
      RangeType alpha_iso = basis_functions_.alpha_iso();

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u_prime)) {
        const auto alpha_prime = cache_iterator->second;
        ret.first = alpha_prime + alpha_iso * std::log(density);
        ret.second = 0.;
        cache_.keep(cache_iterator->first);
        mutex_.unlock();
        return ret;
      } else if (only_cache) {
        DUNE_THROW(Dune::MathError, "Cache was not used!");
      } else {
        RangeType u_iso = basis_functions_.integrated() * 0.5;

        // define further variables
        VectorType g_k, beta_in, beta_out, v;
        beta_in = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;
        static thread_local auto T_k = XT::Common::make_unique<MatrixType>();

        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence) {
          // regularize u
          v = u_prime;
          if (r > 0) {
            beta_in = alpha_iso;
            VectorType r_times_u_iso = u_iso;
            r_times_u_iso *= r;
            v *= 1 - r;
            v += r_times_u_iso;
          }
          *T_k = T_minus_one_;
          // calculate T_k u
          VectorType v_k = v;
          // calculate values of basis p = S_k m
          static thread_local BasisValuesMatrixType P_k = M_;
          P_k = M_;
          // calculate f_0
          RangeFieldType f_k = calculate_scalar_integral(beta_in, P_k);
          f_k -= beta_in * v_k;

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
            if (kk > k_0_ && r < r_max)
              break;
            try {
              change_basis(beta_in, v_k, P_k, *T_k, g_k, beta_out);
            } catch (const Dune::MathError&) {
              if (r < r_max)
                break;
              mutex_.unlock();
              DUNE_THROW(Dune::MathError, "Failure to converge!");
            }
            // calculate current error
            VectorType error(0);
            T_k->mv(g_k, error);
            // calculate descent direction d_k;
            VectorType d_k = g_k;
            d_k *= -1;
            // Calculate stopping criteria (in original basis). Variables with _k are in current basis, without k in
            // original basis.
            RangeType alpha_tilde;
            XT::LA::solve_lower_triangular_transposed(*T_k, alpha_tilde, beta_out);
            auto u_alpha_tilde_k = g_k + v_k;
            RangeType u_alpha_tilde;
            T_k->mv(u_alpha_tilde_k, u_alpha_tilde);
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            const auto alpha_prime = alpha_tilde - alpha_iso * std::log(density_tilde);
            RangeType u_alpha_prime;
            calculate_vector_integral(alpha_prime, M_, M_, u_alpha_prime);
            auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
            VectorType d_alpha_tilde;
            XT::LA::solve_lower_triangular_transposed(*T_k, d_alpha_tilde, d_k);
            if (error.two_norm() < tau_
                && 1 - epsilon_gamma_ < std::exp(d_alpha_tilde.one_norm() + std::abs(std::log(density_tilde)))
                && is_realizable(u_eps_diff)) {
              ret.first = alpha_prime + alpha_iso * std::log(density);
              ret.second = r;
              cache_.insert(v, alpha_prime);
              mutex_.unlock();
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
              // if (zeta_k <= epsilon_ * beta_out.two_norm() / d_k.two_norm() * 100.)
              if (zeta_k <= epsilon_ * beta_out.two_norm() / d_k.two_norm())
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)

        mutex_.unlock();
        DUNE_THROW(MathError, "Failed to converge");
      } // else ( value has not been calculated before )

    outside_all_loops:
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
      const auto alpha = get_alpha(x_local, u, param, true, false).first;
      auto& work_vecs = working_storage();
      calculate_scalar_products(alpha, M_, work_vecs);
      apply_exponential(work_vecs);
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      for (size_t ll = 0; ll < quad_weights_.size(); ++ll) {
        const auto factor = work_vecs[ll] * quad_weights_[ll] * quad_points_[ll][col];
        for (size_t ii = 0; ii < dimRange; ++ii)
          ret[ii] += M_.get_entry(ii, ll) * factor;
      } // ll
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, false, true).first;
      thread_local auto H = XT::Common::make_unique<MatrixType>();
      calculate_hessian(alpha, M_, *H);
      helper<dimDomain>::partial_u(M_, *H, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, false, true).first;
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
    static void calculate_A_Binv(MatrixType& A, MatrixType& B, bool L_calculated = false)
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
      for (size_t ll = 0; ll < num_quad_points; ++ll)
        work_vec[ll] *= quad_weights_[ll];
      // matrix is symmetric, we only use lower triangular part
      const auto& M_backend = M.backend();
      for (size_t ii = 0; ii < dimRange; ++ii) {
        for (size_t kk = 0; kk <= ii; ++kk) {
          RangeFieldType result(0);
          for (size_t ll = 0; ll < num_quad_points; ++ll)
            result += M_backend.get_entry_ref(ii, ll) * M_backend.get_entry_ref(kk, ll) * work_vec[ll];
          H[ii][kk] = result;
        } // kk
      } // ii
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) * quad_weights values
    void calculate_J(const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<RangeFieldType, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimRangeCols);
      const auto& work_vecs = working_storage();
      std::fill(J_dd.begin(), J_dd.end(), 0);
      const auto& M_backend = M.backend();
      const size_t num_quad_points = quad_points_.size();
      for (size_t ii = 0; ii < dimRange; ++ii) {
        for (size_t kk = 0; kk <= ii; ++kk) {
          RangeFieldType result(0.);
          for (size_t ll = 0; ll < num_quad_points; ++ll)
            result += M_backend.get_entry_ref(ii, ll) * M_backend.get_entry_ref(kk, ll) * work_vecs[ll]
                      * quad_points_[ll][dd];
          J_dd[ii][kk] = result;
        } // kk
      } // ii
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
                      VectorType& beta_out) const
    {
      thread_local auto H = XT::Common::make_unique<MatrixType>(0.);
      calculate_hessian(beta_in, P_k, *H);
      XT::LA::cholesky(*H);
      const auto& L = *H;
      T_k.rightmultiply(L);
      L.mtv(beta_in, beta_out);
      StateRangeType tmp_vec;
      XT::LA::solve_lower_triangular(L, tmp_vec, v_k);
      v_k = tmp_vec;
      apply_inverse_matrix(L, P_k);
      calculate_vector_integral(beta_out, P_k, P_k, g_k);
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
    std::mutex& mutex_;
    XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp_;
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
                                           T_minus_one_,
                                           cache_[index],
                                           mutexes_[index],
                                           lp_);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType& u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    const bool boundary = static_cast<bool>(param_neighbor.get("boundary")[0]);
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param, false, true).first;
    const auto alpha_j =
        local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, boundary, !boundary).first;
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
      for (size_t ii = 0; ii < dimRange; ++ii)
        ret[ii] += M_.get_entry(ii, ll) * factor;
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
  const MatrixType T_minus_one_;
  const std::string name_;
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see
  // constructor)
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<std::mutex> mutexes_;
  mutable XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>> lp_;
};
#endif

#if 1

// explicit specialization for 4, 4 because dune-common's operator+/-= do not compile if the two sizes are equal
template <class FieldType>
FieldVector<FieldVector<FieldType, 4>, 4>& operator-=(FieldVector<FieldVector<FieldType, 4>, 4>& vec1,
                                                      const FieldVector<FieldVector<FieldType, 4>, 4>& vec2)
{
  for (size_t ii = 0; ii < vec1.size(); ++ii)
    vec1[ii] -= vec2[ii];
  return vec1;
}

template <class FieldType>
FieldVector<FieldVector<FieldType, 4>, 4>& operator+=(FieldVector<FieldVector<FieldType, 4>, 4>& vec1,
                                                      const FieldVector<FieldVector<FieldType, 4>, 4>& vec2)
{
  for (size_t ii = 0; ii < vec1.size(); ++ii)
    vec1[ii] += vec2[ii];
  return vec1;
}

// explicit specialization for 2, 2 because dune-common's operator+/-= do not compile if the two sizes are equal
template <class FieldType>
FieldVector<FieldVector<FieldType, 2>, 2>& operator-=(FieldVector<FieldVector<FieldType, 2>, 2>& vec1,
                                                      const FieldVector<FieldVector<FieldType, 2>, 2>& vec2)
{
  for (size_t ii = 0; ii < vec1.size(); ++ii)
    vec1[ii] -= vec2[ii];
  return vec1;
}

template <class FieldType>
FieldVector<FieldVector<FieldType, 2>, 2>& operator+=(FieldVector<FieldVector<FieldType, 2>, 2>& vec1,
                                                      const FieldVector<FieldVector<FieldType, 2>, 2>& vec2)
{
  for (size_t ii = 0; ii < vec1.size(); ++ii)
    vec1[ii] += vec2[ii];
  return vec1;
}

/**
 * Specialization for DG basis
 */
template <class GridLayerImp, class U, size_t domainDim>
class EntropyBasedLocalFlux<Hyperbolic::Problems::PiecewiseMonomials<typename U::DomainFieldType,
                                                                     domainDim,
                                                                     typename U::RangeFieldType,
                                                                     U::dimRange,
                                                                     1>,
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
  typedef typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                                   typename U::DomainFieldType,
                                                                   GridLayerImp::dimension,
                                                                   U,
                                                                   0,
                                                                   typename U::RangeFieldType,
                                                                   U::dimRange,
                                                                   GridLayerImp::dimension>
      BaseType;
  typedef EntropyBasedLocalFlux ThisType;

public:
  typedef GridLayerImp GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  static const size_t block_size = (dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = dimRange / block_size;
  using BlockMatrixType = FieldVector<FieldMatrix<RangeFieldType, block_size, block_size>, num_blocks>;
  typedef FieldVector<FieldVector<RangeFieldType, block_size>, num_blocks> BlockVectorType;
  using BasisValuesMatrixType = FieldVector<XT::LA::CommonDenseMatrix<RangeFieldType>, num_blocks>;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureRuleType;
  using QuadraturePointsType =
      FieldVector<std::vector<DomainType, boost::alignment::aligned_allocator<DomainType, 64>>, num_blocks>;
  using QuadratureWeightsType =
      FieldVector<std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>, num_blocks>;
  typedef Hyperbolic::Problems::PiecewiseMonomials<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimDomain>
      BasisfunctionType;
  typedef std::pair<BlockVectorType, RangeFieldType> AlphaReturnType;
  typedef EntropyLocalCache<StateRangeType, BlockVectorType> LocalCacheType;
  using TemporaryVectorType =
      FieldVector<std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>, num_blocks>;
  static const size_t cache_size = 2 * dimDomain + 2;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& quadrature,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52),
      const FieldMatrix<RangeFieldType, block_size, block_size>& T_minus_one =
          XT::LA::eye_matrix<FieldMatrix<RangeFieldType, block_size, block_size>>(block_size,
                                                                                  XT::LA::dense_pattern(block_size,
                                                                                                        block_size)),
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
    , T_minus_one_(T_minus_one)
    , name_(name)
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , mutexes_(index_set_.size(0))
  {
    for (size_t ii = 0; ii < quadrature.size(); ++ii) {
      const auto face_indices = basis_functions_.get_face_indices(quadrature[ii].position());
      const size_t num_adjacent_faces = face_indices.size();
      for (const auto& kk : face_indices) {
        quad_points_[kk].emplace_back(quadrature[ii].position());
        quad_weights_[kk].emplace_back(quadrature[ii].weight() / num_adjacent_faces);
      }
    } // ii
    size_t num_faces;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      while (quad_weights_[jj].size() % 8) { // align to 64 byte boundary
        quad_points_[jj].push_back(quad_points_[jj].back());
        quad_weights_[jj].push_back(0.);
      }
      M_[jj] = XT::LA::CommonDenseMatrix<RangeFieldType>(block_size, quad_points_[jj].size(), 0., 0);
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = basis_functions_.evaluate(quad_points_[jj][ll], false, num_faces);
        for (size_t ii = 0; ii < block_size; ++ii)
          M_[jj].set_entry(ii, ll, val[block_size * jj + ii]);
      } // ii
    } // jj
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using LocalfunctionType::dimRange;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

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
                  const BlockMatrixType& T_minus_one,
                  LocalCacheType& cache,
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
      , mutex_(mutex)
    {
    }

    using LocalfunctionType::entity;

    // temporary vectors to store inner products and exponentials
    TemporaryVectorType& working_storage() const
    {
      thread_local TemporaryVectorType work_vecs;
      for (size_t jj = 0; jj < num_blocks; ++jj)
        work_vecs[jj].resize(quad_points_[jj].size());
      return work_vecs;
    }

    RangeFieldType one_norm(const BlockVectorType& vec) const
    {
      RangeFieldType ret(0);
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ii = 0; ii < block_size; ++ii)
          ret += std::abs(vec[jj][ii]);
      return ret;
    }

    RangeFieldType two_norm(const BlockVectorType& vec) const
    {
      RangeFieldType ret(0);
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ii = 0; ii < block_size; ++ii)
          ret += std::pow(vec[jj][ii], 2);
      ret = std::sqrt(ret);
      return ret;
    }

    void vector_to_block_vector(const StateRangeType& vec, BlockVectorType& block_vec) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ii = 0; ii < block_size; ++ii)
          block_vec[jj][ii] = vec[jj * block_size + ii];
    }

    void block_vector_to_vector(const BlockVectorType& block_vec, StateRangeType& vec) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ii = 0; ii < block_size; ++ii)
          vec[jj * block_size + ii] = block_vec[jj][ii];
    }

    void copy_basis_matrix(const BasisValuesMatrixType& source_mat, BasisValuesMatrixType& range_mat) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        range_mat[jj].backend() = source_mat[jj].backend();
    }

    void calculate_scalar_products(const BlockVectorType& beta_in,
                                   const BasisValuesMatrixType& M,
                                   TemporaryVectorType& scalar_products) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t num_quad_points = quad_points_[jj].size();
        const auto& M_backend = M[jj].backend();
        std::fill(scalar_products[jj].begin(), scalar_products[jj].end(), 0.);
        for (size_t ii = 0; ii < block_size; ++ii) {
          const auto factor = beta_in[jj][ii];
          for (size_t ll = 0; ll < num_quad_points; ++ll)
            scalar_products[jj][ll] += M_backend.get_entry_ref(ii, ll) * factor;
        }
      }
    }

    void apply_exponential(TemporaryVectorType& values) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        assert(values[jj].size() < std::numeric_limits<int>::max());
        XT::Common::Mkl::exp(static_cast<int>(values[jj].size()), values[jj].data(), values[jj].data());
      }
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
    void calculate_vector_integral(const BlockVectorType& beta_in,
                                   const BasisValuesMatrixType& M1,
                                   const BasisValuesMatrixType& M2,
                                   BlockVectorType& ret,
                                   bool only_first_component = false) const
    {
      auto& work_vecs = working_storage();
      calculate_scalar_products(beta_in, M2, work_vecs);
      apply_exponential(work_vecs);
      std::fill(ret.begin(), ret.end(), 0.);
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t num_quad_points = quad_weights_[jj].size();
        const auto& M1_backend = M1[jj].backend();
        for (size_t ii = 0; ii < (only_first_component ? 1 : block_size); ++ii) {
          RangeFieldType result(0.);
          for (size_t ll = 0; ll < num_quad_points; ++ll)
            result += M1_backend.get_entry_ref(ii, ll) * work_vecs[jj][ll] * quad_weights_[jj][ll];
          ret[jj][ii] = result;
        } // ii
      } // jj
    }

    void apply_inverse_matrix(const BlockMatrixType& T_k, BasisValuesMatrixType& M) const
    {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        assert(quad_points_[jj].size() < std::numeric_limits<int>::max());
        XT::Common::Blas::dtrsm(XT::Common::Blas::row_major(),
                                XT::Common::Blas::left(),
                                XT::Common::Blas::lower(),
                                XT::Common::Blas::no_trans(),
                                XT::Common::Blas::non_unit(),
                                block_size,
                                static_cast<int>(quad_points_[jj].size()),
                                1.,
                                &(T_k[jj][0][0]),
                                block_size,
                                M[jj].data(),
                                static_cast<int>(quad_points_[jj].size()));
      } // jj
    }

    AlphaReturnType get_alpha(const DomainType& /*x_local*/,
                              const StateRangeType& u_in,
                              const XT::Common::Parameter& param,
                              const bool regularize,
                              const bool only_cache = false) const
    {
      const bool boundary = static_cast<bool>(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      AlphaReturnType ret;
      StateRangeType v_in;
      mutex_.lock();
      if (boundary)
        cache_.increase_capacity(2 * cache_size);
      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_in);
      if (cache_iterator != cache_.end() && cache_iterator->first == u_in) {
        ret.first = cache_iterator->second;
        ret.second = 0.;
        cache_.keep(cache_iterator->first);
        mutex_.unlock();
        return ret;
      } else if (only_cache) {
        DUNE_THROW(Dune::MathError, "Cache was not used!");
      } else {
        StateRangeType u_iso_in, alpha_iso_in;
        std::tie(u_iso_in, alpha_iso_in) = basis_functions_.calculate_isotropic_distribution(u_in);
        BlockVectorType alpha_iso;
        vector_to_block_vector(alpha_iso_in, alpha_iso);

        // define further variables
        BlockVectorType g_k, beta_in, beta_out, v;
        thread_local auto T_k = XT::Common::make_unique<BlockMatrixType>();
        beta_in = cache_iterator != cache_.end() ? cache_iterator->second : alpha_iso;

        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence) {
          // normalize u
          v_in = u_in;
          if (r > 0.) {
            beta_in = alpha_iso;
            StateRangeType r_times_u_iso = u_iso_in;
            r_times_u_iso *= r;
            v_in *= 1 - r;
            v_in += r_times_u_iso;
          }
          vector_to_block_vector(v_in, v);
          *T_k = T_minus_one_;
          // calculate T_k u
          BlockVectorType v_k = v;
          // calculate values of basis p = S_k m
          thread_local BasisValuesMatrixType P_k(XT::LA::CommonDenseMatrix<RangeFieldType>(0, 0, 0., 0));
          copy_basis_matrix(M_, P_k);
          // calculate f_0
          RangeFieldType f_k = calculate_scalar_integral(beta_in, P_k);
          f_k -= beta_in * v_k;

          int pure_newton = 0;
          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used or cholesky decomposition fails
            if (kk > k_0_ && r < r_max && r_max > 0)
              break;
            try {
              change_basis(beta_in, v_k, P_k, *T_k, g_k, beta_out);
            } catch (const Dune::MathError&) {
              if (r < r_max)
                break;
              mutex_.unlock();
              DUNE_THROW(Dune::MathError, "Failure to converge!");
            }
            // calculate vector of errors e = \int { m exp(beta_out * T_k^{-1}m) } - v
            BlockVectorType error_vec(FieldVector<RangeFieldType, block_size>(0.));
            thread_local BasisValuesMatrixType tmp_mat(XT::LA::CommonDenseMatrix<RangeFieldType>(0, 0, 0., 0));
            copy_basis_matrix(M_, tmp_mat); // calculate T_k^{-1} M_ and store in tmp_mat
            apply_inverse_matrix(*T_k, tmp_mat);
            calculate_vector_integral(beta_out, M_, tmp_mat, error_vec);
            error_vec -= v;

            // calculate descent direction d_k;
            BlockVectorType d_k = g_k;
            d_k *= -1.;
            BlockVectorType T_k_inv_transp_d_k;
            for (size_t jj = 0; jj < num_blocks; ++jj)
              XT::LA::solve_lower_triangular_transposed((*T_k)[jj], T_k_inv_transp_d_k[jj], d_k[jj]);
            if (two_norm(error_vec) < tau_ && std::exp(5 * one_norm(T_k_inv_transp_d_k)) < 1 + epsilon_gamma_) {
              for (size_t jj = 0; jj < num_blocks; ++jj)
                XT::LA::solve_lower_triangular_transposed((*T_k)[jj], ret.first[jj], beta_out[jj]);
              ret.second = r;
              if (kk > 40)
                std::cout << XT::Common::to_string(kk) << std::endl;
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              beta_in = beta_out;
              // backtracking line search
              while (pure_newton >= 2 || zeta_k > epsilon_ * two_norm(beta_out) / two_norm(d_k) * 100.) {
                BlockVectorType beta_new = d_k;
                beta_new *= zeta_k;
                beta_new += beta_out;
                RangeFieldType f = calculate_scalar_integral(beta_new, P_k);
                f -= beta_new * v_k;
                if (pure_newton >= 2 || f <= f_k + xi_ * zeta_k * (g_k * d_k)) {
                  beta_in = beta_new;
                  f_k = f;
                  pure_newton = 0;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
              if (zeta_k <= epsilon_ * two_norm(beta_out) / two_norm(d_k) * 100)
                ++pure_newton;
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)

        mutex_.unlock();
        DUNE_THROW(MathError, "Failed to converge");

      outside_all_loops:
        // store values as initial conditions for next time step on this entity
        cache_.insert(v_in, ret.first);
        mutex_.unlock();
      } // else ( value has not been calculated before )
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
      const auto alpha = get_alpha(x_local, u, param, true).first;
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
      const auto alpha = get_alpha(x_local, u, param, true).first;
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
            ret[offset + ii] += M_[jj].get_entry(ii, ll) * factor;
        } // ll
      } // jj
    } // void evaluate_col_helper(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, false, true).first;
      thread_local auto H = XT::Common::make_unique<BlockMatrixType>();
      calculate_hessian(alpha, M_, *H);
      helper<dimDomain>::partial_u(M_, *H, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, false, true).first;
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
    }; // class helper<3, ...>

    // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
    static void
    calculate_A_Binv(FieldMatrix<RangeFieldType, dimRange, dimRange>& A, BlockMatrixType& B, bool L_calculated = false)
    {
      // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
      // calculate B = LL^T first
      if (!L_calculated) {
        for (size_t jj = 0; jj < num_blocks; ++jj)
          XT::LA::cholesky(B[jj]);
      }
      FieldVector<RangeFieldType, block_size> tmp_vec;
      FieldVector<RangeFieldType, block_size> tmp_A_row;
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        // calculate C = A (L^T)^{-1}
        const auto offset = block_size * jj;
        for (size_t ii = 0; ii < block_size; ++ii) {
          for (size_t kk = 0; kk < block_size; ++kk)
            tmp_A_row[kk] = A[offset + ii][offset + kk];
          XT::LA::solve_lower_triangular(B[jj], tmp_vec, tmp_A_row);
          // calculate ret = C L^{-1}
          XT::LA::solve_lower_triangular_transposed(B[jj], tmp_A_row, tmp_vec);
          for (size_t kk = 0; kk < block_size; ++kk)
            A[offset + ii][offset + kk] = tmp_A_row[kk];
        } // ii
      } // jj
    } // void calculate_A_Binv(...)

    void calculate_hessian(const BlockVectorType& alpha, const BasisValuesMatrixType& M, BlockMatrixType& H) const
    {
      auto& work_vecs = working_storage();
      calculate_scalar_products(alpha, M, work_vecs);
      apply_exponential(work_vecs);
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
          work_vecs[jj][ll] *= quad_weights_[jj][ll];
      // matrix is symmetric, we only use lower triangular part
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const size_t num_quad_points = quad_weights_[jj].size();
        const auto& M_backend = M[jj].backend();
        for (size_t ii = 0; ii < block_size; ++ii) {
          for (size_t kk = 0; kk <= ii; ++kk) {
            RangeFieldType result(0);
            for (size_t ll = 0; ll < num_quad_points; ++ll)
              result += M_backend.get_entry_ref(ii, ll) * M_backend.get_entry_ref(kk, ll) * work_vecs[jj][ll];
            H[jj][ii][kk] = result;
          } // kk
        } // ii
      } // jj
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    // assumes work_vecs already contains the needed exp(alpha * m) * quad_weights values
    void calculate_J(const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<RangeFieldType, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimRangeCols);
      const auto& work_vecs = working_storage();
      std::fill(J_dd.begin(), J_dd.end(), 0.);
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = jj * block_size;
        const auto& M_backend = M[jj].backend();
        const size_t num_quad_points = quad_points_[jj].size();
        for (size_t ii = 0; ii < block_size; ++ii) {
          for (size_t kk = 0; kk <= ii; ++kk) {
            RangeFieldType result(0.);
            for (size_t ll = 0; ll < num_quad_points; ++ll)
              result += M_backend.get_entry_ref(ii, ll) * M_backend.get_entry_ref(kk, ll) * work_vecs[jj][ll]
                        * quad_points_[jj][ll][dd];
            J_dd[offset + ii][offset + kk] = result;
          } // kk
        } // ii
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
                      BlockVectorType& beta_out) const
    {
      thread_local auto H = XT::Common::make_unique<BlockMatrixType>(0.);
      calculate_hessian(beta_in, P_k, *H);
      for (size_t jj = 0; jj < num_blocks; ++jj)
        XT::LA::cholesky((*H)[jj]);
      const auto& L = *H;
      FieldVector<RangeFieldType, block_size> tmp_vec;
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        T_k[jj].rightmultiply(L[jj]);
        L[jj].mtv(beta_in[jj], beta_out[jj]);
        XT::LA::solve_lower_triangular(L[jj], tmp_vec, v_k[jj]);
        v_k[jj] = tmp_vec;
      } // jj
      apply_inverse_matrix(L, P_k);
      calculate_vector_integral(beta_out, P_k, P_k, g_k, true);
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
    const BlockMatrixType& T_minus_one_;
    const std::string name_;
    // constructor)
    LocalCacheType& cache_;
    std::mutex& mutex_;
  }; // class Localfunction

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
                                           mutexes_[index]);
  }


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(XT::Common::FloatCmp::ne(n_ij[dd], 0.));
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param, false, true).first;
    const auto alpha_j =
        local_function_neighbor
            ->get_alpha(
                x_local_neighbor, u_j, param_neighbor, false, !static_cast<bool>(param_neighbor.get("boundary")[0]))
            .first;
    thread_local FieldVector<TemporaryVectorType, 2> work_vecs;
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
          ret[offset + ii] += M_[jj].get_entry(ii, ll) * factor;
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
  const BlockMatrixType T_minus_one_;
  const std::string name_;
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<std::mutex> mutexes_;
};
#endif

template <class FieldType, int size>
FieldVector<FieldType, size> pow(const FieldVector<FieldType, size>& base, const size_t exponent)
{
  FieldVector<FieldType, size> ret;
  for (int ii = 0; ii < size; ++ii)
    ret[ii] = std::pow(base[ii], exponent);
  return ret;
}

template <class FieldType, int size>
FieldVector<FieldType, size> multiply_componentwise(const FieldVector<FieldType, size>& first,
                                                    const FieldVector<FieldType, size>& second)
{
  auto ret = first;
  for (int ii = 0; ii < size; ++ii)
    ret[ii] *= second[ii];
  return ret;
}

template <class FieldType, int size>
FieldVector<DynamicVector<FieldType>, size>
multiply_componentwise(const FieldVector<DynamicVector<FieldType>, size>& first,
                       const FieldVector<DynamicVector<FieldType>, size>& second)
{
  auto ret = first;
  for (int ii = 0; ii < 3; ++ii)
    for (size_t kk = 0; kk < second[ii].size(); ++kk)
      ret[ii][kk] *= second[ii][kk];
  return ret;
}

template <class FieldType>
DynamicVector<FieldType>& operator*=(DynamicVector<FieldType>& first, const DynamicVector<FieldType>& second)
{
  assert(first.size() == second.size());
  for (size_t ii = 0; ii < second.size(); ++ii)
    first[ii] *= second[ii];
  return first;
}

template <class FieldType>
DynamicVector<FieldVector<FieldType, 3>>& operator*=(DynamicVector<FieldVector<FieldType, 3>>& first,
                                                     const DynamicVector<FieldVector<FieldType, 3>>& second)
{
  assert(first.size() == second.size());
  for (size_t ii = 0; ii < second.size(); ++ii)
    for (size_t jj = 0; jj < 3; ++jj)
      first[ii][jj] *= second[ii][jj];
  return first;
}

template <class FieldType>
DynamicVector<FieldMatrix<FieldType, 3, 3>>& operator*=(DynamicVector<FieldMatrix<FieldType, 3, 3>>& first,
                                                        const DynamicVector<FieldMatrix<FieldType, 3, 3>>& second)
{
  assert(first.size() == second.size());
  for (size_t ii = 0; ii < second.size(); ++ii)
    for (size_t jj = 0; jj < 3; ++jj)
      for (size_t kk = 0; kk < 3; ++kk)
        first[ii][jj][kk] *= second[ii][jj][kk];
  return first;
}

#if 1
/**
 * Specialization of EntropyBasedLocalFlux for 3D Hatfunctions
 */
template <class GridLayerImp, class U>
class EntropyBasedLocalFlux<Hyperbolic::Problems::HatFunctions<typename U::DomainFieldType,
                                                               3,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               1>,
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
  typedef GridLayerImp GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using BasisValuesMatrixType = std::vector<std::pair<StateRangeType, RangeFieldType>>;
  typedef FieldVector<RangeFieldType, 3> ThreeVectorType;
  typedef FieldMatrix<RangeFieldType, 3, 3> ThreeMatrixType;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureRuleType;
  typedef Hyperbolic::Problems::HatFunctions<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimDomain>
      BasisfunctionType;
  // index1: order l, index2: l1, index3: k, index4: qq
  typedef std::vector<std::vector<DynamicVector<FieldVector<RangeFieldType, 3>>>> P1Type;
  // index1: order l, index2: l1, index3: j \in \{k_1, k_2, k_3}, index4: k, index5: qq
  typedef std::vector<std::vector<FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>>> P2Type;
  // index1: dimension d, index2 :order l, index3: l1, index4: j \in \{k_1, k_2, k_3}, index5: k, index6: qq
  typedef FieldVector<std::vector<std::vector<FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>>>, 3>
      P3Type;
  // index1: order l, index2: l1, index3: j \in \{k_1, k_2, k_3}, index4: m \in \{k_1, k_2, k_3}, index5: k, index6: qq
  typedef std::vector<std::vector<std::array<std::array<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>, 3>>> P4Type;
  // index1: dimension d, index2 :order l, index3: l1, index4: j \in \{k_1, k_2, k_3}, index5: m \in \{k_1, k_2, k_3},
  // index6: k, index7: qq
  typedef FieldVector<std::vector<std::vector<std::array<std::array<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>,
                                                         3>>>,
                      3>
      P5Type;
  using AlphaReturnType = typename std::pair<StateRangeType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, AlphaReturnType>;
  static constexpr size_t cache_size = 2 * dimDomain + 2;

private:
  class PowCache
  {
  public:
    void update(const StateRangeType& alpha)
    {
      if (alpha != last_alpha_) {
        pow_cache_.resize(2);
        pow_cache_[0] = StateRangeType(1.);
        pow_cache_[1] = alpha;
        last_alpha_ = alpha;
      }
    }

    void add(const size_t ll)
    {
      if (ll == pow_cache_.size())
        pow_cache_.push_back(multiply_componentwise(pow_cache_[ll / 2], pow_cache_[ll - ll / 2]));
    }

    const StateRangeType& get(const size_t ll)
    {
      return pow_cache_[ll];
    }

  private:
    StateRangeType last_alpha_;
    std::vector<StateRangeType> pow_cache_;
  };

  class BlockedPowCache
  {
  public:
    const std::vector<size_t>&
    update(const StateRangeType& alpha, const size_t num_blocks, const FieldVector<std::vector<size_t>, 3>& k)
    {
      if (alpha != last_alpha_) {
        qqs_.resize(num_blocks);
        pow_cache_.resize(2,
                          FieldVector<DynamicVector<RangeFieldType>, 2>(DynamicVector<RangeFieldType>(num_blocks, 1.)));
        for (size_t kk = 0; kk < num_blocks; ++kk) {
          qqs_[kk] = 0;
          for (size_t qq = 1; qq < 3; ++qq)
            if (alpha[k[qq][kk]] < alpha[k[qqs_[kk]][kk]])
              qqs_[kk] = qq;
          pow_cache_[1][0][kk] = alpha[k[(qqs_[kk] + 1) % 3][kk]] - alpha[k[qqs_[kk]][kk]];
          pow_cache_[1][1][kk] = alpha[k[(qqs_[kk] + 2) % 3][kk]] - alpha[k[qqs_[kk]][kk]];
        } // kk
        last_alpha_ = alpha;
      } // if (alpha != last_alpha_)
      return qqs_;
    }

    void add(const size_t ll)
    {
      if (ll == pow_cache_.size())
        pow_cache_.push_back(multiply_componentwise(pow_cache_[ll / 2], pow_cache_[ll - ll / 2]));
    }

    const DynamicVector<RangeFieldType>& get_first(const size_t ll) const
    {
      return pow_cache_[ll][0];
    }

    const DynamicVector<RangeFieldType>& get_second(const size_t ll) const
    {
      return pow_cache_[ll][1];
    }

  private:
    StateRangeType last_alpha_;
    std::vector<FieldVector<DynamicVector<RangeFieldType>, 2>> pow_cache_;
    std::vector<size_t> qqs_;
  };

  static bool not_equal(const FieldVector<DynamicVector<RangeFieldType>, 3>& first,
                        const FieldVector<DynamicVector<RangeFieldType>, 3>& second,
                        const RangeFieldType tol)
  {
    for (size_t ii = 0; ii < 3; ++ii)
      if (XT::Common::FloatCmp::ne(first[ii], second[ii], tol))
        return true;
    return false;
  } // bool not_equal(..)

  static bool not_equal(const std::array<std::array<DynamicVector<RangeFieldType>, 3>, 3>& first,
                        const std::array<std::array<DynamicVector<RangeFieldType>, 3>, 3>& second,
                        const RangeFieldType tol)
  {
    for (size_t ii = 0; ii < 3; ++ii)
      for (size_t jj = 0; jj < 3; ++jj)
        if (XT::Common::FloatCmp::ne(first[ii][jj], second[ii][jj], tol))
          return true;
    return false;
  } // bool not_equal(..)

public:
  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& quadrature,
      const size_t max_order = 150,
      const RangeFieldType tol = 1e-10,
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
    , quadrature_(quadrature)
    , max_order_(max_order)
    , tol_(tol)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , name_(name)
    , num_blocks_(basis_functions_.triangulation().faces().size())
    , p1_(max_order_ + 1)
    , p2_(max_order_ + 1)
    , p3_minus_(std::vector<std::vector<FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>>>(max_order_ + 1))
    , p3_plus_(p3_minus_)
    , p4_(max_order_ + 1)
    , p5_(std::vector<std::vector<std::array<std::array<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>, 3>>>(
          max_order_ + 1))
    , k_(std::vector<size_t>(num_blocks_, 0))
    , cache_(index_set_.size(0), LocalCacheType(cache_size))
    , mutexes_(index_set_.size(0))
  {
    std::vector<BasisValuesMatrixType> M(num_blocks_);
    std::vector<std::vector<size_t>> quadrature_indices(num_blocks_);
    for (size_t ii = 0; ii < quadrature_.size(); ++ii) {
      const auto face_indices = basis_functions_.triangulation().get_face_indices(quadrature_[ii].position());
      const size_t num_adjacent_faces = face_indices.size();
      for (const auto& kk : face_indices) {
        M[kk].push_back(std::make_pair(basis_functions_.evaluate(quadrature_[ii].position()),
                                       quadrature_[ii].weight() / num_adjacent_faces));
        quadrature_indices[kk].push_back(ii);
      } // face_indices
    } // ii

    // precalculate factors for the adaptive quadrature. We need to calculate the integrals
    // < exp(alpha^T b) >, < b exp(alpha^T b) >, < v b exp(alpha^T b) >, < b b^T exp(alpha^T b) >, < v b b^T exp(alpha^T
    // b) >,
    // which will be approximated as
    // \sum_{l=0}^{order} \sum_{l_1 + l_2 + l_3 = l} \sum_{k=1}^{N} \prod_{i=1}^{3} \frac{alpha_{k_i}^{l_i}}{l_i!} p(k,
    // l_1, l_2, l_3),
    // where N is the number of faces of the triangulation, order is the taylor order, S_k is the k-th face of the
    // triangulation,
    // k_1, k_2, k_3 are the indices of the three basis_functions that belong to the vertices of S_k and the factor p(k,
    // l_1, l_2, l_3)
    // can be precomputed and is
    // p1(k, l_1, l_2, l_3) = \int_{S_k} \prod_{i=1}^3 b_{k_i}^{l_i} for the integral < exp(alpha^T b) >
    // p2(k, l_1, l_2, l_3) = \int_{S_k} \prod_{i=1}^3 b_{k_i}^{l_i + \delta_{k_i,j}} for the integral < b_j exp(alpha^T
    // b) >
    // p3plus(k, l_1, l_2, l_3) = \int_{S_k} (v_d)_+ \prod_{i=1}^3 b_{k_i}^{l_i + \delta_{k_i,j}} for the integral < v_d
    // b_j exp(alpha^T b) >_+exp(-11.7413)
    // p3minus(k, l_1, l_2, l_3) = \int_{S_k} (v_d)_- \prod_{i=1}^3 b_{k_i}^{l_i + \delta_{k_i,j}} for the integral <
    // v_d b_j exp(alpha^T b) >_-
    // p4(k, l_1, l_2, l_3) = \int_{S_k} \prod_{i=1}^3 b_{k_i}^{l_i + \delta_{k_i,j} + \delta_{k_i,m}} for the integral
    // < b_j b_m exp(alpha^T b) >
    // p5(k, l_1, l_2, l_3) = \int_{S_k} \prod_{i=1}^3 v_d b_{k_i}^{l_i + \delta_{k_i,j} + \delta_{k_i,m}} for the
    // integral < v_d b_j b_m exp(alpha^T b) >
    // where b_i is the i-th basis function.
    // now walk over faces of the triangulation and precompute factors
    const auto& faces = basis_functions_.triangulation().faces();
    PowCache pow_cache;
    std::vector<RangeFieldType> factorial_cache(2, 1.);
    factorial_cache.reserve(max_order_ + 1);
    bool resized = false;
    for (const auto& face : faces) {
      const size_t kk = face->index();
      const auto& vertices = face->vertices();
      assert(vertices.size() == 3);
      for (size_t ii = 0; ii < 3; ++ii)
        k_[ii][kk] = vertices[ii]->index();
      for (size_t qq = 0; qq < 3; ++qq) {
        for (size_t ii = 0; ii < M[kk].size(); ++ii) {
          const auto& h = M[kk][ii].first;
          pow_cache.update(h);
          const auto& v = quadrature_[quadrature_indices[kk][ii]].position();

          for (size_t ll = 0; ll <= max_order_; ++ll) {
            pow_cache.add(ll);
            if (factorial_cache.size() == ll)
              factorial_cache.push_back(factorial_cache[ll - 1] * ll);
            if (!resized) {
              const size_t l1_size = ll + 1;
              p1_[ll].resize(
                  l1_size,
                  DynamicVector<FieldVector<RangeFieldType, 3>>(num_blocks_, FieldVector<RangeFieldType, 3>(0.)));
              p2_[ll].resize(
                  l1_size,
                  FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>(
                      DynamicVector<FieldVector<RangeFieldType, 3>>(num_blocks_, FieldVector<RangeFieldType, 3>(0.))));
              p4_[ll].resize(l1_size);
              for (size_t l1 = 0; l1 <= ll; ++l1)
                std::fill_n(
                    p4_[ll][l1][0].data(),
                    9,
                    DynamicVector<FieldVector<RangeFieldType, 3>>(num_blocks_, FieldVector<RangeFieldType, 3>(0.)));
              for (size_t dd = 0; dd < dimDomain; ++dd) {
                p3_minus_[dd][ll].resize(l1_size,
                                         FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>(
                                             DynamicVector<FieldVector<RangeFieldType, 3>>(
                                                 num_blocks_, FieldVector<RangeFieldType, 3>(0.))));
                p3_plus_[dd][ll].resize(l1_size,
                                        FieldVector<DynamicVector<FieldVector<RangeFieldType, 3>>, 3>(
                                            DynamicVector<FieldVector<RangeFieldType, 3>>(
                                                num_blocks_, FieldVector<RangeFieldType, 3>(0.))));
                p5_[dd][ll].resize(l1_size);
                for (size_t l1 = 0; l1 <= ll; ++l1)
                  std::fill_n(
                      p5_[dd][ll][l1][0].data(),
                      9,
                      DynamicVector<FieldVector<RangeFieldType, 3>>(num_blocks_, FieldVector<RangeFieldType, 3>(0.)));
              }
            } // if (!resized)
            for (size_t l1 = 0; l1 <= ll; ++l1) {
              const size_t l2 = ll - l1;
              const auto product = pow_cache.get(l1)[k_[(qq + 1) % 3][kk]] / factorial_cache[l1]
                                   * pow_cache.get(l2)[k_[(qq + 2) % 3][kk]] / factorial_cache[l2] * M[kk][ii].second;
              p1_[ll][l1][kk][qq] += product;
              for (size_t jj = 0; jj < 3; ++jj) {
                auto product_h_jj = h[k_[jj][kk]] * product;
                p2_[ll][l1][jj][kk][qq] += product_h_jj;
                for (size_t mm = 0; mm < 3; ++mm)
                  p4_[ll][l1][jj][mm][kk][qq] += h[k_[mm][kk]] * product_h_jj;
                for (size_t dd = 0; dd < dimDomain; ++dd) {
                  auto product_h_jj_v_dd = v[dd] * product_h_jj;
                  if (XT::Common::FloatCmp::eq(v[dd], 0.)) {
                    p3_minus_[dd][ll][l1][jj][kk][qq] += 0.5 * product_h_jj_v_dd;
                    p3_plus_[dd][ll][l1][jj][kk][qq] += 0.5 * product_h_jj_v_dd;
                  } else if (v[dd] > 0.) {
                    p3_plus_[dd][ll][l1][jj][kk][qq] += product_h_jj_v_dd;
                  } else {
                    p3_minus_[dd][ll][l1][jj][kk][qq] += product_h_jj_v_dd;
                  }
                  for (size_t mm = 0; mm < 3; ++mm)
                    p5_[dd][ll][l1][jj][mm][kk][qq] += h[k_[mm][kk]] * product_h_jj_v_dd;
                } // dd
              } // jj
            } // l1
          } // ll
          resized = true;
        } // ii
      } // qq
    } // faces)
  } // constructor

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const QuadratureRuleType& quadrature,
                  const size_t max_order,
                  const RangeFieldType tol,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const size_t num_blocks,
                  const P1Type& p1,
                  const P2Type& p2,
                  const P3Type& p3_minus,
                  const P3Type& p3_plus,
                  const P4Type& p4,
                  const P5Type& p5,
                  const FieldVector<std::vector<size_t>, 3>& k,
                  LocalCacheType& cache,
                  BlockedPowCache& blocked_pow_cache,
                  std::mutex& mutex)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
      , quadrature_(quadrature)
      , max_order_(max_order)
      , tol_(tol)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , num_blocks_(num_blocks)
      , p1_(p1)
      , p2_(p2)
      , p3_minus_(p3_minus)
      , p3_plus_(p3_plus)
      , p4_(p4)
      , p5_(p5)
      , k_(k)
      , cache_(cache)
      , blocked_pow_cache_(blocked_pow_cache)
      , mutex_(mutex)
    {
    }

    using LocalfunctionType::entity;

    AlphaReturnType get_alpha(const DomainType& /*x_local*/,
                              const StateRangeType& u,
                              const XT::Common::Parameter& param,
                              const bool regularize = true) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      AlphaReturnType ret;
      mutex_.lock();
      if (boundary)
        cache_.increase_capacity(2 * cache_size);
      const auto cache_iterator = cache_.find_closest(u);
      if (cache_iterator != cache_.end() && XT::Common::FloatCmp::eq(cache_iterator->first, u)) {
        ret = cache_iterator->second;
        mutex_.unlock();
        return ret;
      } else {
        // define further variables
        thread_local std::unique_ptr<MatrixType> H_k = XT::Common::make_unique<MatrixType>();

        // calculate moment vector for isotropic distribution
        StateRangeType u_iso, alpha_iso, v;
        std::tie(u_iso, alpha_iso) = basis_functions_.calculate_isotropic_distribution(u);
        StateRangeType alpha_k = cache_iterator != cache_.end() ? cache_iterator->second.first : alpha_iso;

        const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
        const auto r_max = r_sequence.back();
        for (const auto& r : r_sequence) {
          // get initial alpha
          v = u;
          if (r > 0) {
            alpha_k = alpha_iso;
            StateRangeType r_times_u_iso = u_iso;
            r_times_u_iso *= r;
            v *= 1 - r;
            v += r_times_u_iso;
          }

          // calculate f_0
          RangeFieldType f_k = calculate_f(alpha_k, v);

          for (size_t kk = 0; kk < k_max_; ++kk) {
            // exit inner for loop to increase r if too many iterations are used
            if (kk > k_0_ && r < r_max && r_max > 0)
              break;

            // calculate gradient g
            StateRangeType g_k = calculate_gradient(alpha_k, v);
            // calculate Hessian H
            calculate_hessian(alpha_k, *H_k);
            // calculate descent direction d_k;
            StateRangeType d_k(0);
            StateRangeType minus_g_k(g_k);
            minus_g_k *= -1;
            try {
              XT::LA::cholesky(*H_k);
              const auto& L = *H_k;
              StateRangeType tmp_vec;
              XT::LA::solve_lower_triangular(L, tmp_vec, minus_g_k);
              XT::LA::solve_lower_triangular_transposed(L, d_k, tmp_vec);
            } catch (const Dune::MathError&) {
              if (r < r_max) {
                break;
              } else {
                mutex_.unlock();
                DUNE_THROW(Dune::MathError, "Failure to converge!");
              }
            }
            if (g_k.two_norm() < tau_ && std::exp(5 * d_k.one_norm()) < 1 + epsilon_gamma_) {
              ret = std::make_pair(alpha_k, r);
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              // backtracking line search
              while (zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
                // calculate alpha_new = alpha_k + zeta_k d_k
                auto alpha_new = d_k;
                alpha_new *= zeta_k;
                alpha_new += alpha_k;
                RangeFieldType f_new = calculate_f(alpha_new, v);
                if (XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  alpha_k = alpha_new;
                  f_k = f_new;
                  break;
                }
                zeta_k = chi_ * zeta_k;
              } // backtracking linesearch while
            } // else (stopping conditions)
          } // k loop (Newton iterations)
        } // r loop (Regularization parameter)

        mutex_.unlock();
        DUNE_THROW(MathError, "Failed to converge");

      outside_all_loops:
        cache_.insert(v, ret);
        mutex_.unlock();
      } // else ( value has not been calculated before )

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
      ColRangeType ret_col(0.);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        evaluate_col(dd, x_local, u, ret_col, param);
        for (size_t ii = 0; ii < dimRange; ++ii)
          ret[ii][dd] = ret_col[ii];
      } // dd
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      std::fill(ret.begin(), ret.end(), 0.);
      const auto alpha = get_alpha(x_local, u, param).first;
      const auto& qqs = blocked_pow_cache_.update(alpha, num_blocks_, k_);
      DynamicVector<RangeFieldType> exp_alpha3(num_blocks_);
      for (size_t kk = 0; kk < num_blocks_; ++kk)
        exp_alpha3[kk] = std::exp(alpha[k_[qqs[kk]][kk]]);
      StateRangeType update(1.);
      FieldVector<DynamicVector<RangeFieldType>, 3> update_vec((DynamicVector<RangeFieldType>(num_blocks_)));
      DynamicVector<RangeFieldType> tmp_vec(num_blocks_);
      StateRangeType zero_vector(0.);
      size_t ll = 0;
      while (ll <= max_order_ && XT::Common::FloatCmp::ne(update, zero_vector, tol_)) {
        blocked_pow_cache_.add(ll);
        // calculate_integral
        std::fill(update_vec.begin(), update_vec.end(), 0.);
        std::fill(update.begin(), update.end(), 0.);
        for (size_t l1 = 0; l1 <= ll; ++l1) {
          const size_t l2 = ll - l1;
          tmp_vec = blocked_pow_cache_.get_first(l1);
          const auto& alpha2m3_pow_l2 = blocked_pow_cache_.get_second(l2);
          tmp_vec *= alpha2m3_pow_l2;
          for (size_t jj = 0; jj < 3; ++jj)
            for (size_t kk = 0; kk < num_blocks_; ++kk)
              update_vec[jj][kk] +=
                  tmp_vec[kk] * (p3_plus_[col][ll][l1][jj][kk][qqs[kk]] + p3_minus_[col][ll][l1][jj][kk][qqs[kk]]);
        }
        for (size_t jj = 0; jj < 3; ++jj)
          for (size_t kk = 0; kk < num_blocks_; ++kk)
            update[k_[jj][kk]] += update_vec[jj][kk] * exp_alpha3[kk];
        ret += update;
        ++ll;
      } // ll
    } // void evaluate_col(...)

    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param).first;
      thread_local std::unique_ptr<MatrixType> H = XT::Common::make_unique<MatrixType>(0.);
      calculate_hessian(alpha, *H);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        calculate_J(alpha, ret[dd], dd);
        calculate_A_Binv(ret[dd], *H, dd > 0);
      }
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param).first;
      thread_local std::unique_ptr<MatrixType> H = XT::Common::make_unique<MatrixType>(0.);
      calculate_hessian(alpha, *H);
      calculate_J(alpha, ret, col);
      calculate_A_Binv(ret, *H);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    RangeFieldType calculate_f(const StateRangeType& alpha, const StateRangeType& v) const
    {
      RangeFieldType ret(0);
      const auto& qqs = blocked_pow_cache_.update(alpha, num_blocks_, k_);
      DynamicVector<RangeFieldType> tmp_vec(num_blocks_, 1.);
      DynamicVector<RangeFieldType> tmp_update_vec(num_blocks_, 1.);
      DynamicVector<RangeFieldType> exp_alpha3(num_blocks_);
      for (size_t kk = 0; kk < num_blocks_; ++kk)
        exp_alpha3[kk] = std::exp(alpha[k_[qqs[kk]][kk]]);
      DynamicVector<RangeFieldType> update_vec(num_blocks_, 0.);
      static const auto zero_vec = update_vec;
      size_t ll = 0;
      while (ll <= max_order_ && XT::Common::FloatCmp::ne(tmp_update_vec, zero_vec, tol_)) {
        blocked_pow_cache_.add(ll);
        // calculate_integral
        std::fill(tmp_update_vec.begin(), tmp_update_vec.end(), 0.);
        for (size_t l1 = 0; l1 <= ll; ++l1) {
          const size_t l2 = ll - l1;
          tmp_vec = blocked_pow_cache_.get_first(l1);
          const auto& alpha2m3_pow_l2 = blocked_pow_cache_.get_second(l2);
          tmp_vec *= alpha2m3_pow_l2;
          for (size_t kk = 0; kk < num_blocks_; ++kk)
            tmp_update_vec[kk] += tmp_vec[kk] * p1_[ll][l1][kk][qqs[kk]];
        } // l1
        update_vec += tmp_update_vec;
        ++ll;
      } // ll
      update_vec *= exp_alpha3;
      ret = std::accumulate(update_vec.begin(), update_vec.end(), 0.);
      ret -= alpha * v;
      return ret;
    } // .. calculate_f(...)

    StateRangeType calculate_gradient(const StateRangeType& alpha, const StateRangeType& v) const
    {
      StateRangeType ret(0);
      const auto& qqs = blocked_pow_cache_.update(alpha, num_blocks_, k_);
      DynamicVector<RangeFieldType> exp_alpha3(num_blocks_);
      for (size_t kk = 0; kk < num_blocks_; ++kk)
        exp_alpha3[kk] = std::exp(alpha[k_[qqs[kk]][kk]]);
      FieldVector<DynamicVector<RangeFieldType>, 3> update_vec((DynamicVector<RangeFieldType>(num_blocks_, 0.)));
      FieldVector<DynamicVector<RangeFieldType>, 3> tmp_update_vec((DynamicVector<RangeFieldType>(num_blocks_, 1.)));
      static const auto zero_update_vec = update_vec;
      DynamicVector<RangeFieldType> tmp_vec(num_blocks_);
      size_t ll = 0;
      while (ll <= max_order_ && not_equal(tmp_update_vec, zero_update_vec, tol_)) {
        blocked_pow_cache_.add(ll);
        // calculate_integral
        std::fill(tmp_update_vec.begin(), tmp_update_vec.end(), 0.);
        for (size_t l1 = 0; l1 <= ll; ++l1) {
          const size_t l2 = ll - l1;
          tmp_vec = blocked_pow_cache_.get_first(l1);
          const auto& alpha2m3_pow_l2 = blocked_pow_cache_.get_second(l2);
          tmp_vec *= alpha2m3_pow_l2;
          for (size_t jj = 0; jj < 3; ++jj)
            for (size_t kk = 0; kk < num_blocks_; ++kk)
              tmp_update_vec[jj][kk] += tmp_vec[kk] * p2_[ll][l1][jj][kk][qqs[kk]];
        } // l1
        for (size_t jj = 0; jj < 3; ++jj)
          update_vec[jj] += tmp_update_vec[jj];
        ++ll;
      } // ll
      for (size_t jj = 0; jj < 3; ++jj)
        for (size_t kk = 0; kk < num_blocks_; ++kk)
          ret[k_[jj][kk]] += update_vec[jj][kk] * exp_alpha3[kk];
      ret -= v;
      return ret;
    }

    void calculate_hessian(const StateRangeType& alpha, MatrixType& H) const
    {
      std::fill(H.begin(), H.end(), 0.);
      const auto& qqs = blocked_pow_cache_.update(alpha, num_blocks_, k_);
      DynamicVector<RangeFieldType> exp_alpha3(num_blocks_);
      for (size_t kk = 0; kk < num_blocks_; ++kk)
        exp_alpha3[kk] = std::exp(alpha[k_[qqs[kk]][kk]]);
      std::array<std::array<DynamicVector<RangeFieldType>, 3>, 3> update_vecs;
      std::fill_n(update_vecs[0].data(), 9, DynamicVector<RangeFieldType>(num_blocks_, 0.));
      auto tmp_update_vecs = update_vecs;
      DynamicVector<RangeFieldType> tmp_vec(num_blocks_);
      static const auto zero_update_vecs = update_vecs;
      size_t ll = 0;
      tmp_update_vecs[0][0][0] = 1.;
      while (ll <= max_order_ && not_equal(tmp_update_vecs, zero_update_vecs, tol_)) {
        blocked_pow_cache_.add(ll);
        for (auto& row : tmp_update_vecs)
          for (auto& vec : row)
            std::fill(vec.begin(), vec.end(), 0.);
        for (size_t l1 = 0; l1 <= ll; ++l1) {
          const size_t l2 = ll - l1;
          tmp_vec = blocked_pow_cache_.get_first(l1);
          const auto& alpha2m3_pow_l2 = blocked_pow_cache_.get_second(l2);
          tmp_vec *= alpha2m3_pow_l2;
          for (size_t jj = 0; jj < 3; ++jj)
            for (size_t mm = 0; mm < 3; ++mm)
              for (size_t kk = 0; kk < num_blocks_; ++kk)
                tmp_update_vecs[jj][mm][kk] += tmp_vec[kk] * p4_[ll][l1][jj][mm][kk][qqs[kk]];
        } // l1
        for (size_t jj = 0; jj < 3; ++jj)
          for (size_t mm = 0; mm < 3; ++mm)
            update_vecs[jj][mm] += tmp_update_vecs[jj][mm];
        ++ll;
      } // ll
      for (size_t jj = 0; jj < 3; ++jj)
        for (size_t mm = 0; mm < 3; ++mm)
          for (size_t kk = 0; kk < num_blocks_; ++kk)
            H[k_[jj][kk]][k_[mm][kk]] += update_vecs[jj][mm][kk] * exp_alpha3[kk];
    } // void calculate_hessian(...)

    void calculate_J(const StateRangeType& alpha, MatrixType& J, const size_t dd) const
    {
      std::fill(J.begin(), J.end(), 0.);
      const auto& qqs = blocked_pow_cache_.update(alpha, num_blocks_, k_);
      DynamicVector<RangeFieldType> exp_alpha3(num_blocks_);
      for (size_t kk = 0; kk < num_blocks_; ++kk)
        exp_alpha3[kk] = std::exp(alpha[k_[qqs[kk]][kk]]);
      std::array<std::array<DynamicVector<RangeFieldType>, 3>, 3> update_vecs;
      std::fill_n(update_vecs[0].data(), 9, DynamicVector<RangeFieldType>(num_blocks_, 0.));
      auto tmp_update_vecs = update_vecs;
      DynamicVector<RangeFieldType> tmp_vec(num_blocks_);
      static const auto zero_update_vecs = update_vecs;
      size_t ll = 0;
      tmp_update_vecs[0][0][0] = 1.;
      while (ll <= max_order_ && not_equal(tmp_update_vecs, zero_update_vecs, tol_)) {
        blocked_pow_cache_.add(ll);
        for (auto& row : tmp_update_vecs)
          for (auto& vec : row)
            std::fill(vec.begin(), vec.end(), 0.);
        for (size_t l1 = 0; l1 <= ll; ++l1) {
          const size_t l2 = ll - l1;
          tmp_vec = blocked_pow_cache_.get_first(l1);
          const auto& alpha2m3_pow_l2 = blocked_pow_cache_.get_second(l2);
          tmp_vec *= alpha2m3_pow_l2;
          for (size_t jj = 0; jj < 3; ++jj)
            for (size_t mm = 0; mm < 3; ++mm)
              for (size_t kk = 0; kk < num_blocks_; ++kk)
                tmp_update_vecs[jj][mm][kk] += tmp_vec[kk] * p5_[dd][ll][l1][jj][mm][kk][qqs[kk]];
        } // l1
        for (size_t jj = 0; jj < 3; ++jj)
          for (size_t mm = 0; mm < 3; ++mm)
            update_vecs[jj][mm] += tmp_update_vecs[jj][mm];
        ++ll;
      } // ll
      for (size_t jj = 0; jj < 3; ++jj)
        for (size_t mm = 0; mm < 3; ++mm)
          for (size_t kk = 0; kk < num_blocks_; ++kk)
            J[k_[jj][kk]][k_[mm][kk]] += update_vecs[jj][mm][kk] * exp_alpha3[kk];
    } // void calculate_J(...)

    // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
    static void calculate_A_Binv(MatrixType& A, MatrixType& B, bool L_calculated = false)
    {
      // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
      // calculate B = LL^T
      if (!L_calculated)
        XT::LA::cholesky(B);
      const auto& L = B;
      StateRangeType tmp_vec;
      for (size_t ii = 0; ii < dimRange; ++ii) {
        // calculate C = A (L^T)^{-1} and store in B
        XT::LA::solve_lower_triangular(L, tmp_vec, A[ii]);
        // calculate ret = C L^{-1}
        XT::LA::solve_lower_triangular_transposed(L, A[ii], tmp_vec);
      } // ii
    } // void calculate_A_Binv(...)

    const BasisfunctionType& basis_functions_;
    const QuadratureRuleType& quadrature_;
    const size_t max_order_;
    const RangeFieldType tol_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const size_t num_blocks_;
    const P1Type& p1_;
    const P2Type& p2_;
    const P3Type& p3_minus_;
    const P3Type& p3_plus_;
    const P4Type& p4_;
    const P5Type& p5_;
    const FieldVector<std::vector<size_t>, 3>& k_;
    const std::string name_;
    LocalCacheType& cache_;
    BlockedPowCache& blocked_pow_cache_;
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
                                           quadrature_,
                                           max_order_,
                                           tol_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           num_blocks_,
                                           p1_,
                                           p2_,
                                           p3_minus_,
                                           p3_plus_,
                                           p4_,
                                           p5_,
                                           k_,
                                           cache_[index],
                                           blocked_pow_cache_,
                                           mutexes_[index]);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType& u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param, false).first;
    thread_local BlockedPowCache pow_cache_i;
    pow_cache_i = blocked_pow_cache_;
    const auto alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor, false).first;
    auto& pow_cache_j = blocked_pow_cache_;
    StateRangeType ret(0);
    auto& pow_cache_pos = n_ij[dd] > 0 ? pow_cache_i : pow_cache_j;
    auto& pow_cache_neg = n_ij[dd] > 0 ? pow_cache_j : pow_cache_i;
    auto& alpha_pos = n_ij[dd] > 0 ? alpha_i : alpha_j;
    auto& alpha_neg = n_ij[dd] > 0 ? alpha_j : alpha_i;
    const auto& qqs_pos = pow_cache_pos.update(alpha_pos, num_blocks_, k_);
    const auto& qqs_neg = pow_cache_neg.update(alpha_neg, num_blocks_, k_);
    DynamicVector<RangeFieldType> exp_alpha3_pos(num_blocks_);
    auto exp_alpha3_neg = exp_alpha3_pos;
    for (size_t kk = 0; kk < num_blocks_; ++kk) {
      exp_alpha3_pos[kk] = std::exp(alpha_pos[k_[qqs_pos[kk]][kk]]);
      exp_alpha3_neg[kk] = std::exp(alpha_neg[k_[qqs_neg[kk]][kk]]);
    }
    FieldVector<DynamicVector<RangeFieldType>, 3> update_vec_pos((DynamicVector<RangeFieldType>(num_blocks_, 0.)));
    FieldVector<DynamicVector<RangeFieldType>, 3> update_vec_neg((DynamicVector<RangeFieldType>(num_blocks_, 0.)));
    FieldVector<DynamicVector<RangeFieldType>, 3> tmp_update_vec_pos((DynamicVector<RangeFieldType>(num_blocks_, 1.)));
    FieldVector<DynamicVector<RangeFieldType>, 3> tmp_update_vec_neg((DynamicVector<RangeFieldType>(num_blocks_, 1.)));
    static const auto zero_update_vec = update_vec_pos;
    DynamicVector<RangeFieldType> tmp_vec_pos(num_blocks_);
    DynamicVector<RangeFieldType> tmp_vec_neg(num_blocks_);
    size_t ll = 0;
    while (ll <= max_order_ && (not_equal(tmp_update_vec_pos, zero_update_vec, tol_)
                                || not_equal(tmp_update_vec_neg, zero_update_vec, tol_))) {
      pow_cache_pos.add(ll);
      pow_cache_neg.add(ll);
      // calculate_integral
      std::fill(tmp_update_vec_pos.begin(), tmp_update_vec_pos.end(), 0.);
      std::fill(tmp_update_vec_neg.begin(), tmp_update_vec_neg.end(), 0.);
      for (size_t l1 = 0; l1 <= ll; ++l1) {
        const size_t l2 = ll - l1;
        tmp_vec_pos = pow_cache_pos.get_first(l1);
        tmp_vec_neg = pow_cache_neg.get_first(l1);
        const auto& alpha2m3_pow_l2_pos = pow_cache_pos.get_second(l2);
        const auto& alpha2m3_pow_l2_neg = pow_cache_neg.get_second(l2);
        tmp_vec_pos *= alpha2m3_pow_l2_pos;
        tmp_vec_neg *= alpha2m3_pow_l2_neg;
        for (size_t jj = 0; jj < 3; ++jj)
          for (size_t kk = 0; kk < num_blocks_; ++kk) {
            tmp_update_vec_pos[jj][kk] += tmp_vec_pos[kk] * p3_plus_[dd][ll][l1][jj][kk][qqs_pos[kk]];
            tmp_update_vec_neg[jj][kk] += tmp_vec_neg[kk] * p3_minus_[dd][ll][l1][jj][kk][qqs_neg[kk]];
          }
      } // l1
      for (size_t jj = 0; jj < 3; ++jj) {
        update_vec_pos[jj] += tmp_update_vec_pos[jj];
        update_vec_neg[jj] += tmp_update_vec_neg[jj];
      }
      ++ll;
    } // ll
    for (size_t jj = 0; jj < 3; ++jj) {
      for (size_t kk = 0; kk < num_blocks_; ++kk) {
        ret[k_[jj][kk]] += update_vec_pos[jj][kk] * exp_alpha3_pos[kk];
        ret[k_[jj][kk]] += update_vec_neg[jj][kk] * exp_alpha3_neg[kk];
      } // kk
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
  const QuadratureRuleType quadrature_;
  const size_t max_order_;
  const RangeFieldType tol_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const std::string name_;
  const size_t num_blocks_;
  P1Type p1_;
  P2Type p2_;
  P3Type p3_minus_;
  P3Type p3_plus_;
  P4Type p4_;
  P5Type p5_;
  // contains the three basis function indices k_1, k_2, k_3 for each face k
  FieldVector<std::vector<size_t>, 3> k_;
  mutable std::vector<LocalCacheType> cache_;
  mutable std::vector<std::mutex> mutexes_;
  static thread_local BlockedPowCache blocked_pow_cache_;
};

template <class GridLayerImp, class U>
thread_local typename EntropyBasedLocalFlux<Hyperbolic::Problems::HatFunctions<typename U::DomainFieldType,
                                                                               3,
                                                                               typename U::RangeFieldType,
                                                                               U::dimRange,
                                                                               1>,
                                            GridLayerImp,
                                            U>::BlockedPowCache
    EntropyBasedLocalFlux<Hyperbolic::Problems::
                              HatFunctions<typename U::DomainFieldType, 3, typename U::RangeFieldType, U::dimRange, 1>,
                          GridLayerImp,
                          U>::blocked_pow_cache_;
#endif

#if 1
/**
 * Specialization of EntropyBasedLocalFlux for 1D Hatfunctions
 */
template <class GridLayerImp, class U>
class EntropyBasedLocalFlux<Hyperbolic::Problems::HatFunctions<typename U::DomainFieldType,
                                                               1,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               1,
                                                               1>,
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
  typedef typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                                   typename U::DomainFieldType,
                                                                   GridLayerImp::dimension,
                                                                   U,
                                                                   0,
                                                                   typename U::RangeFieldType,
                                                                   U::dimRange,
                                                                   1>
      BaseType;
  typedef EntropyBasedLocalFlux ThisType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using BasisfunctionType = Hyperbolic::Problems::HatFunctions<DomainFieldType, 1, RangeFieldType, dimRange, 1, 1>;
  using GridLayerType = GridLayerImp;
  using QuadratureRuleType = Dune::QuadratureRule<DomainFieldType, 1>;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using AlphaReturnType = typename std::pair<StateRangeType, RangeFieldType>;
  using LocalCacheType = EntropyLocalCache<StateRangeType, StateRangeType>;
  static const size_t cache_size = 2 * dimDomain + 2;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& /*quadrature*/,
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
    , mutexes_(index_set_.size(0))
  {
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const RangeType& v_points,
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
      , mutex_(mutex)
    {
    }

    using LocalfunctionType::entity;

    static bool is_realizable(const RangeType& u, const RangeFieldType eps = 0.)
    {
      for (const auto& u_i : u)
        if (u_i < eps)
          return false;
      return true;
    }

    AlphaReturnType get_alpha(const DomainType& /*x_local*/,
                              const StateRangeType& u,
                              const XT::Common::Parameter& param,
                              const bool regularize,
                              const bool only_cache) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      AlphaReturnType ret;
      mutex_.lock();
      if (boundary)
        cache_.increase_capacity(2 * cache_size);

      // rescale u such that the density <psi> is 1
      RangeFieldType density = basis_functions_.density(u);
      RangeType u_prime = u / density;
      RangeType alpha_iso = basis_functions_.alpha_iso();

      // if value has already been calculated for these values, skip computation
      const auto cache_iterator = cache_.find_closest(u_prime);
      if (cache_iterator != cache_.end() && cache_iterator->first == u_prime) {
        const auto alpha_prime = cache_iterator->second;
        ret.first = alpha_prime + alpha_iso * std::log(density);
        ret.second = 0.;
        cache_.keep(cache_iterator->first);
        mutex_.unlock();
        return ret;
      } else if (only_cache) {
        DUNE_THROW(Dune::MathError, "Cache was not used!");
      } else {
        // The hessian H is always symmetric and tridiagonal, so we only need to store the diagonal and subdiagonal
        // elements
        RangeType H_diag;
        FieldVector<RangeFieldType, dimRange - 1> H_subdiag;

        // calculate moment vector for isotropic distribution
        RangeType u_iso = basis_functions_.integrated() * 0.5;
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
                DUNE_THROW(Dune::MathError, "Failure to converge!");
              }
            }


            const auto& alpha_tilde = alpha_k;
            const auto u_alpha_tilde = g_k + v;
            auto density_tilde = basis_functions_.density(u_alpha_tilde);
            const auto alpha_prime = alpha_tilde - alpha_iso * std::log(density_tilde);
            const auto u_alpha_prime = calculate_u(alpha_prime);
            auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
            // checking realizability is cheap so we do not need the second stopping criterion
            if (g_k.two_norm() < tau_ && is_realizable(u_eps_diff)) {
              ret.first = alpha_prime + alpha_iso * std::log(density);
              ret.second = r;
              cache_.insert(v, alpha_prime);
              mutex_.unlock();
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
      const auto alpha = get_alpha(x_local, u, param, true, false).first;

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
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param, false, true).first;
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
                       + 2. * ((2 * v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn])
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
    const RangeType& v_points_;
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
                                           mutexes_[index]);
  }

  // calculate \sum_{i=1}^d < v_i_+ m \psi >_+ n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t DXTC_DEBUG_ONLY(dd),
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(dd == 0);
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param, false, true).first;
    const auto alpha_j =
        local_function_neighbor
            ->get_alpha(
                x_local_neighbor, u_j, param_neighbor, false, !static_cast<bool>(param_neighbor.get("boundary")[0]))
            .first;
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
  const RangeType& v_points_;
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
  mutable std::vector<std::mutex> mutexes_;
};
#endif


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
