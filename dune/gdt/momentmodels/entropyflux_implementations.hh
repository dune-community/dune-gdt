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

#ifndef DUNE_GDT_MOMENTMODELS_ENTROPYFLUX_IMPLEMENTATIONS_HH
#define DUNE_GDT_MOMENTMODELS_ENTROPYFLUX_IMPLEMENTATIONS_HH

#include <algorithm>
#include <cmath>
#include <list>
#include <memory>

#include <boost/align/aligned_allocator.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/container/pattern.hh>

#include <dune/xt/functions/interfaces/function.hh>

#include <dune/gdt/momentmodels/basisfunctions.hh>
#include <dune/gdt/type_traits.hh>

#if HAVE_CLP
#  include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

namespace Dune {
namespace GDT {


// choose specializations
#ifndef ENTROPY_FLUX_UNSPECIALIZED_USE_ADAPTIVE_CHANGE_OF_BASIS
#  define ENTROPY_FLUX_UNSPECIALIZED_USE_ADAPTIVE_CHANGE_OF_BASIS 1
#endif

#ifndef ENTROPY_FLUX_USE_PARTIAL_MOMENTS_SPECIALIZATION
#  define ENTROPY_FLUX_USE_PARTIAL_MOMENTS_SPECIALIZATION 1
#endif

#ifndef ENTROPY_FLUX_USE_3D_HATFUNCTIONS_SPECIALIZATION
#  define ENTROPY_FLUX_USE_3D_HATFUNCTIONS_SPECIALIZATION 1
#endif

#ifndef ENTROPY_FLUX_USE_1D_HATFUNCTIONS_SPECIALIZATION
#  define ENTROPY_FLUX_USE_1D_HATFUNCTIONS_SPECIALIZATION 1
#endif

#ifndef ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS
#  define ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS 0
#endif


template <class MomentBasisImp>
class EntropyBasedFluxImplementationUnspecializedBase
  : public XT::Functions::FunctionInterface<MomentBasisImp::dimRange,
                                            MomentBasisImp::dimDomain,
                                            MomentBasisImp::dimRange,
                                            typename MomentBasisImp::R>
{
  using BaseType = typename XT::Functions::FunctionInterface<MomentBasisImp::dimRange,
                                                             MomentBasisImp::dimDomain,
                                                             MomentBasisImp::dimRange,
                                                             typename MomentBasisImp::R>;
  using ThisType = EntropyBasedFluxImplementationUnspecializedBase;

public:
  using MomentBasis = MomentBasisImp;
  using BaseType::d;
  using BaseType::r;
  static const size_t basis_dimDomain = MomentBasis::dimDomain;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  // make matrices a little larger to align to 64 byte boundary
  static constexpr size_t matrix_num_cols =
      basis_dimRange % 8 ? basis_dimRange : basis_dimRange + (8 - basis_dimRange % 8);
  using MatrixType = XT::Common::FieldMatrix<RangeFieldType, basis_dimRange, basis_dimRange>;
  using VectorType = XT::Common::FieldVector<RangeFieldType, basis_dimRange>;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using BasisValuesMatrixType = XT::LA::CommonDenseMatrix<RangeFieldType>;
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;

  explicit EntropyBasedFluxImplementationUnspecializedBase(const MomentBasis& basis_functions,
                                                           const RangeFieldType tau,
                                                           const RangeFieldType epsilon_gamma,
                                                           const RangeFieldType chi,
                                                           const RangeFieldType xi,
                                                           const std::vector<RangeFieldType> r_sequence,
                                                           const size_t k_0,
                                                           const size_t k_max,
                                                           const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , quad_points_(XT::Data::merged_quadrature(basis_functions_.quadratures()).size())
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
    , realizability_helper_(basis_functions_,
                            quad_points_
#if HAVE_CLP
                            ,
                            lp_
#endif
      )
  {
    size_t ll = 0;
    for (const auto& quad_point : XT::Data::merged_quadrature(basis_functions_.quadratures())) {
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
      for (size_t ii = 0; ii < basis_dimRange; ++ii)
        M_.set_entry(ll, ii, val[ii]);
    }
  }

  virtual int order(const XT::Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  VectorType get_isotropic_alpha(const DomainType& u) const
  {
    static const auto alpha_iso = basis_functions_.alpha_iso();
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    return alpha_iso + alpha_iso_prime * std::log(basis_functions_.density(u));
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    auto& work_vecs = working_storage();
    calculate_scalar_products(alpha, M_, work_vecs);
    apply_exponential(work_vecs);
    for (size_t dd = 0; dd < basis_dimDomain; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t ll = 0; ll < quad_weights_.size(); ++ll) {
        const auto factor = work_vecs[ll] * quad_weights_[ll] * quad_points_[ll][dd];
        for (size_t ii = 0; ii < basis_dimRange; ++ii)
          ret[dd][ii] += M_.get_entry(ll, ii) * factor;
      } // ll
    } // dd
    return ret;
  } // void evaluate(...)

  virtual void jacobian(const DomainType& u,
                        DynamicDerivativeRangeType& result,
                        const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  virtual void jacobian_with_alpha(const VectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    thread_local auto H = XT::Common::make_unique<MatrixType>();
    calculate_hessian(alpha, M_, *H);
    for (size_t dd = 0; dd < basis_dimDomain; ++dd)
      row_jacobian(dd, M_, *H, result[dd], dd > 0);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType evaluate_kinetic_flux(const DomainType& u_i,
                                   const DomainType& u_j,
                                   const BasisDomainType& n_ij,
                                   const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const BasisDomainType& n_ij,
                                               const size_t dd) const

  {
    thread_local FieldVector<std::vector<RangeFieldType>, 2> work_vecs;
    work_vecs[0].resize(quad_points_.size());
    work_vecs[1].resize(quad_points_.size());
    calculate_scalar_products(alpha_i, M_, work_vecs[0]);
    calculate_scalar_products(alpha_j, M_, work_vecs[1]);
    DomainType ret(0);
    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      const auto position = quad_points_[ll][dd];
      RangeFieldType factor = position * n_ij[dd] > 0. ? std::exp(work_vecs[0][ll]) : std::exp(work_vecs[1][ll]);
      factor *= quad_weights_[ll] * position;
      const auto* basis_ll = M_.get_ptr(ll);
      for (size_t ii = 0; ii < basis_dimRange; ++ii)
        ret[ii] += basis_ll[ii] * factor;
    } // ll
    ret *= n_ij[dd];
    return ret;
  } // DomainType evaluate_kinetic_flux_with_alphas(...)

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const = 0;

protected:
  // get permutation instead of sorting directly to be able to sort two vectors the same way
  // see
  // https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  template <typename T, typename Compare>
  static std::vector<std::size_t> get_sort_permutation(const std::vector<T>& vec, const Compare& compare)
  {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
    return p;
  }

  template <typename T>
  static void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p)
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
  static void join_duplicate_quadpoints(std::vector<BasisDomainType>& quad_points,
                                        std::vector<RangeFieldType>& quad_weights)
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

#if HAVE_CLP
  template <class BasisFuncImp = MomentBasis, bool anything = true>
  struct RealizabilityHelper
  {
    static_assert(std::is_same<BasisFuncImp, MomentBasis>::value, "BasisFuncImp has to be MomentBasis!");

    RealizabilityHelper(const MomentBasis& basis_functions,
                        const std::vector<BasisDomainType>& quad_points,
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
        // We start with creating a model with basis_dimRange rows and num_quad_points columns */
        constexpr int num_rows = static_cast<int>(basis_dimRange);
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

    bool is_realizable(const DomainType& u, const bool reinitialize) const
    {
      const auto density = basis_functions_.density(u);
      if (!(density > 0.) || std::isinf(density))
        return false;
      const auto u_prime = u / density;
      setup_linear_program(reinitialize);
      auto& lp = **lp_;
      constexpr int num_rows = static_cast<int>(basis_dimRange);
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
    const MomentBasis& basis_functions_;
    const std::vector<BasisDomainType>& quad_points_;
    XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp_;
  }; // struct RealizabilityHelper<...>
#else // HAVE_CLP
  template <class BasisFuncImp = MomentBasis, bool anything = true>
  struct RealizabilityHelper
  {
    RealizabilityHelper(const MomentBasis& /*basis_functions*/, const std::vector<BasisDomainType>& /*quad_points*/)
    {
      DUNE_THROW(Dune::NotImplemented, "You are missing Clp!");
    }

    bool is_realizable(const DomainType& /*u*/, const bool /*reinitialize*/) const
    {
      DUNE_THROW(Dune::NotImplemented, "You are missing Clp!");
      return false;
    }
  }; // struct RealizabilityHelper<...>
#endif // HAVE_CLP

  // specialization for hatfunctions
  template <size_t dimRange_or_refinements, bool anything>
  struct RealizabilityHelper<HatFunctionMomentBasis<DomainFieldType,
                                                    basis_dimDomain,
                                                    RangeFieldType,
                                                    dimRange_or_refinements,
                                                    1,
                                                    basis_dimDomain>,
                             anything>
  {
    RealizabilityHelper(const MomentBasis& /*basis_functions*/,
                        const std::vector<BasisDomainType>& /*quad_points*/
#if HAVE_CLP
                        ,
                        XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& /*lp*/)
#else
    )
#endif
    {}

    static bool is_realizable(const DomainType& u, const bool /*reinitialize*/)
    {
      for (const auto& u_i : u)
        if (!(u_i > 0.) || std::isinf(u_i))
          return false;
      return true;
    }
  }; // struct RealizabilityHelper<Hatfunctions, ...>

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
                            basis_dimRange,
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
      for (size_t ii = 0; ii < (only_first_component ? 1 : basis_dimRange); ++ii)
        ret[ii] += basis_ll[ii] * factor_ll;
    } // ll
  }

  void copy_transposed(const MatrixType& T_k, MatrixType& T_k_trans) const
  {
    for (size_t ii = 0; ii < basis_dimRange; ++ii)
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
                            basis_dimRange,
                            1.,
                            &((*T_k_trans)[0][0]),
                            basis_dimRange,
                            M.data(),
                            matrix_num_cols);
#else
    assert(quad_points_.size() == M.rows());
    VectorType tmp_vec, tmp_vec2;
    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      std::copy_n(M.get_ptr(ll), basis_dimRange, tmp_vec.begin());
      XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
      std::copy_n(tmp_vec2.begin(), basis_dimRange, M.get_ptr(ll));
    }
#endif
  }

  void row_jacobian(const size_t row,
                    const BasisValuesMatrixType& M,
                    MatrixType& H,
                    DynamicRowDerivativeRangeType& ret,
                    bool L_calculated = false) const
  {
    assert(row < basis_dimDomain);
    calculate_J(M, ret, row);
    calculate_A_Binv(ret, H, L_calculated);
  } // void partial_u_col(...)

  // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
  static void calculate_A_Binv(DynamicRowDerivativeRangeType& A, MatrixType& B, bool L_calculated = false)
  {
    // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
    // calculate B = LL^T first
    if (!L_calculated)
      XT::LA::cholesky(B);
    VectorType tmp_vec;
    for (size_t ii = 0; ii < basis_dimRange; ++ii) {
      // calculate C = A (L^T)^{-1} and store in B
      auto&& row_view = A[ii];
      XT::LA::solve_lower_triangular(B, tmp_vec, row_view);
      // calculate ret = C L^{-1}
      XT::LA::solve_lower_triangular_transposed(B, row_view, tmp_vec);
    } // ii
  } // void calculate_A_Binv(...)

  void calculate_hessian(const VectorType& alpha,
                         const BasisValuesMatrixType& M,
                         MatrixType& H,
                         const bool use_work_vec_data = false) const
  {
    std::fill(H.begin(), H.end(), 0.);
    auto& work_vec = working_storage();
    if (!use_work_vec_data) {
      calculate_scalar_products(alpha, M, work_vec);
      apply_exponential(work_vec);
    }
    const size_t num_quad_points = quad_weights_.size();
    // matrix is symmetric, we only use lower triangular part
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      auto factor_ll = work_vec[ll] * quad_weights_[ll];
      const auto* basis_ll = M.get_ptr(ll);
      for (size_t ii = 0; ii < basis_dimRange; ++ii) {
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
  void calculate_J(const BasisValuesMatrixType& M, DynamicRowDerivativeRangeType& J_dd, const size_t dd) const
  {
    assert(dd < basis_dimDomain);
    const auto& work_vecs = working_storage();
    J_dd.set_all_entries(0.);
    const size_t num_quad_points = quad_points_.size();
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto factor_ll = work_vecs[ll] * quad_points_[ll][dd] * quad_weights_[ll];
      const auto* basis_ll = M.get_ptr(ll);
      for (size_t ii = 0; ii < basis_dimRange; ++ii) {
        const auto factor_ll_ii = factor_ll * basis_ll[ii];
        if (!XT::Common::is_zero(factor_ll_ii)) {
          for (size_t kk = 0; kk <= ii; ++kk)
            J_dd.add_to_entry(ii, kk, basis_ll[kk] * factor_ll_ii);
        }
      } // ii
    } // ll
    // symmetric update for upper triangular part of J
    for (size_t mm = 0; mm < basis_dimRange; ++mm)
      for (size_t nn = mm + 1; nn < basis_dimRange; ++nn)
        J_dd.set_entry(mm, nn, J_dd.get_entry(nn, mm));
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
    VectorType tmp_vec;
    XT::LA::solve_lower_triangular(L, tmp_vec, v_k);
    v_k = tmp_vec;
    apply_inverse_matrix(L, P_k);
    calculate_vector_integral(beta_out, P_k, P_k, g_k, false);
    g_k -= v_k;
  } // void change_basis(...)

  const MomentBasis& basis_functions_;
  std::vector<BasisDomainType> quad_points_;
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
  const RealizabilityHelper<> realizability_helper_;
#if HAVE_CLP
  mutable XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>> lp_;
#endif
};

#if ENTROPY_FLUX_UNSPECIALIZED_USE_ADAPTIVE_CHANGE_OF_BASIS
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 */
template <class MomentBasisImp>
class EntropyBasedFluxImplementation : public EntropyBasedFluxImplementationUnspecializedBase<MomentBasisImp>
{
  using BaseType = EntropyBasedFluxImplementationUnspecializedBase<MomentBasisImp>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using BaseType::basis_dimDomain;
  using typename BaseType::AlphaReturnType;
  using typename BaseType::BasisValuesMatrixType;
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::VectorType;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : BaseType(basis_functions, tau, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon)
    , T_minus_one_(std::make_unique<MatrixType>())
  {
    XT::LA::eye_matrix(*T_minus_one_);
  }

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const override final
  {
    auto ret = std::make_unique<AlphaReturnType>();

    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");

    VectorType u_prime = u / density;
    VectorType alpha_initial = alpha_in - alpha_iso_prime * std::log(density);
    VectorType beta_in = alpha_initial;
    VectorType v, u_eps_diff, g_k, beta_out;
    RangeFieldType first_error_cond, second_error_cond, tau_prime;

    auto u_iso = basis_functions_.u_iso();
    const RangeFieldType dim_factor = is_full_moment_basis<MomentBasis>::value ? 1. : std::sqrt(basis_dimDomain);
    tau_prime = std::min(tau_ / ((1 + dim_factor * u_prime.two_norm()) * density + dim_factor * tau_), tau_);

    thread_local auto T_k = XT::Common::make_unique<MatrixType>();

    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence) {
      // regularize u
      v = u_prime;
      if (r > 0) {
        beta_in = get_isotropic_alpha(u);
        VectorType r_times_u_iso = u_iso;
        r_times_u_iso *= r;
        v *= 1 - r;
        v += r_times_u_iso;
      }
      *T_k = *T_minus_one_;
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
          const std::string err_msg =
              "Failed to converge for " + XT::Common::to_string(u) + " with density " + XT::Common::to_string(density)
              + " and multiplier " + XT::Common::to_string(beta_in)
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
        const auto alpha_prime = alpha_tilde - alpha_iso_prime * std::log(density_tilde);
        VectorType u_alpha_prime;
        calculate_vector_integral(alpha_prime, M_, M_, u_alpha_prime);
        u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        VectorType d_alpha_tilde;
        XT::LA::solve_lower_triangular_transposed(*T_k, d_alpha_tilde, d_k);
        first_error_cond = g_alpha_tilde.two_norm();
        second_error_cond = std::exp(d_alpha_tilde.one_norm() + std::abs(std::log(density_tilde)));
        if (first_error_cond < tau_prime && 1 - epsilon_gamma_ < second_error_cond
            && realizability_helper_.is_realizable(u_eps_diff, kk == static_cast<size_t>(0.8 * k_0_))) {
          ret->first = alpha_prime + alpha_iso_prime * std::log(density);
          ret->second = std::make_pair(v * density, r);
          return ret;
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
    const std::string err_msg = "Failed to converge for " + XT::Common::to_string(u) + " with density "
                                + XT::Common::to_string(density) + " and multiplier " + XT::Common::to_string(beta_in)
                                + " due to too many iterations! Last u_eps_diff = " + XT::Common::to_string(u_eps_diff)
                                + ", first_error_cond = " + XT::Common::to_string(first_error_cond)
                                + ", second_error_cond = " + XT::Common::to_string(second_error_cond)
                                + ", tau_prime = " + XT::Common::to_string(tau_prime);
    DUNE_THROW(MathError, err_msg);

    return ret;
  }

private:
  using BaseType::apply_inverse_matrix;
  using BaseType::calculate_hessian;
  using BaseType::calculate_scalar_integral;
  using BaseType::calculate_vector_integral;
  using BaseType::get_isotropic_alpha;

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
    VectorType tmp_vec;
    XT::LA::solve_lower_triangular(L, tmp_vec, v_k);
    v_k = tmp_vec;
    apply_inverse_matrix(L, P_k);
    calculate_vector_integral(beta_out, P_k, P_k, g_k, false);
    g_k -= v_k;
  } // void change_basis(...)

  using BaseType::basis_functions_;
  using BaseType::chi_;
  using BaseType::epsilon_;
  using BaseType::epsilon_gamma_;
  using BaseType::k_0_;
  using BaseType::k_max_;
  using BaseType::M_;
  using BaseType::quad_points_;
  using BaseType::quad_weights_;
  using BaseType::r_sequence_;
  using BaseType::realizability_helper_;
  using BaseType::tau_;
  using BaseType::xi_;
  const std::unique_ptr<MatrixType> T_minus_one_;
};

#else // ENTROPY_FLUX_UNSPECIALIZED_USE_ADAPTIVE_CHANGE_OF_BASIS

/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * Simple backtracking Newton without change of basis
 */
template <class MomentBasisImp>
class EntropyBasedFluxImplementation : public EntropyBasedFluxImplementationUnspecializedBase<MomentBasisImp>
{
  using BaseType = EntropyBasedFluxImplementationUnspecializedBase<MomentBasisImp>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using BaseType::basis_dimDomain;
  using typename BaseType::AlphaReturnType;
  using typename BaseType::BasisValuesMatrixType;
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::VectorType;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : BaseType(basis_functions, tau, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon)
  {}

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const override final
  {
    auto ret = std::make_unique<AlphaReturnType>();
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    VectorType u_prime = u / density;
    VectorType alpha_initial = alpha_in - alpha_iso_prime * std::log(density);
    VectorType v, g_k, d_k, tmp_vec, alpha_prime;
    RangeFieldType first_error_cond, second_error_cond, tau_prime;
    auto u_iso = basis_functions_.u_iso();
    const RangeFieldType dim_factor = is_full_moment_basis<MomentBasis>::value ? 1. : std::sqrt(basis_dimDomain);
    tau_prime = std::min(tau_ / ((1 + dim_factor * u_prime.two_norm()) * density + dim_factor * tau_), tau_);
    VectorType alpha_k = alpha_initial;
    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence) {
      // regularize u
      v = u_prime;
      if (r > 0) {
        alpha_k = get_isotropic_alpha(u);
        VectorType r_times_u_iso = u_iso;
        r_times_u_iso *= r;
        v *= 1 - r;
        v += r_times_u_iso;
      }
      // calculate T_k u
      VectorType v_k = v;
      // calculate f_0
      RangeFieldType f_k = calculate_scalar_integral(alpha_k, M_);
      f_k -= alpha_k * v_k;

      thread_local auto H = XT::Common::make_unique<MatrixType>(0.);

      int pure_newton = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase r if too many iterations are used
        if (kk > k_0_ && r < r_max)
          break;
        // calculate gradient g
        calculate_vector_integral(alpha_k, M_, M_, g_k);
        g_k -= v_k;
        // calculate Hessian H
        calculate_hessian(alpha_k, M_, *H, true);
        // calculate descent direction d_k;
        d_k = g_k;
        d_k *= -1;
        try {
          // if H = LL^T, then we have to calculate d_k = - L^{-T} L^{-1} g_k
          // calculate H = LL^T first
          XT::LA::cholesky(*H);
          // calculate d_tmp = -L^{-1} g_k and store in B
          XT::LA::solve_lower_triangular(*H, tmp_vec, d_k);
          // calculate d_k = L^{-T} d_tmp
          XT::LA::solve_lower_triangular_transposed(*H, d_k, tmp_vec);
        } catch (const Dune::MathError&) {
          if (r < r_max)
            break;
          const std::string err_msg =
              "Failed to converge for " + XT::Common::to_string(u) + " with density " + XT::Common::to_string(density);
          DUNE_THROW(MathError, err_msg);
        }

        const auto& alpha_tilde = alpha_k;
        auto& u_alpha_tilde = tmp_vec;
        u_alpha_tilde = g_k + v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;
        alpha_prime = alpha_tilde - alpha_iso_prime * std::log(density_tilde);
        auto& u_eps_diff = tmp_vec;
        calculate_vector_integral(alpha_prime, M_, M_, u_eps_diff);
        u_eps_diff *= -(1 - epsilon_gamma_);
        u_eps_diff += v;

        first_error_cond = g_k.two_norm();
        second_error_cond = std::exp(d_k.one_norm() + std::abs(std::log(density_tilde)));
        if (first_error_cond < tau_prime && 1 - epsilon_gamma_ < second_error_cond
            && realizability_helper_.is_realizable(u_eps_diff, kk == static_cast<size_t>(0.8 * k_0_))) {
          ret->first = alpha_prime + alpha_iso_prime * std::log(density);
          ret->second = std::make_pair(v * density, r);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          auto& alpha_new = tmp_vec;
          while (pure_newton >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
            // calculate alpha_new = alpha_k + zeta_k d_k
            alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;
            // calculate f(alpha_new)
            RangeFieldType f_new = calculate_scalar_integral(alpha_new, M_);
            f_new -= alpha_new * v_k;
            if (pure_newton >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              pure_newton = 0;
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
    const std::string err_msg =
        "Failed to converge for " + XT::Common::to_string(u) + " with density " + XT::Common::to_string(density);
    DUNE_THROW(MathError, err_msg);
    return ret;
  }

private:
  using BaseType::calculate_hessian;
  using BaseType::calculate_scalar_integral;
  using BaseType::calculate_vector_integral;
  using BaseType::get_isotropic_alpha;

  using BaseType::basis_functions_;
  using BaseType::chi_;
  using BaseType::epsilon_;
  using BaseType::epsilon_gamma_;
  using BaseType::k_0_;
  using BaseType::k_max_;
  using BaseType::M_;
  using BaseType::quad_points_;
  using BaseType::quad_weights_;
  using BaseType::r_sequence_;
  using BaseType::realizability_helper_;
  using BaseType::tau_;
  using BaseType::xi_;
};
#endif

#if ENTROPY_FLUX_USE_PARTIAL_MOMENTS_SPECIALIZATION
/**
 * Specialization for DG basis
 */
template <class D, size_t d, class R, size_t dimRange_or_refinements>
class EntropyBasedFluxImplementation<PartialMomentBasis<D, d, R, dimRange_or_refinements, 1>>
  : public XT::Functions::FunctionInterface<PartialMomentBasis<D, d, R, dimRange_or_refinements, 1>::dimRange,
                                            d,
                                            PartialMomentBasis<D, d, R, dimRange_or_refinements, 1>::dimRange,
                                            R>
{
public:
  using MomentBasis = PartialMomentBasis<D, d, R, dimRange_or_refinements, 1>;
  using BaseType = typename XT::Functions::
      FunctionInterface<MomentBasis::dimRange, MomentBasis::dimDomain, MomentBasis::dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;
  using BaseType::d;
  using BaseType::r;
  static const size_t basis_dimDomain = MomentBasis::dimDomain;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  static const size_t block_size = (basis_dimDomain == 1) ? 2 : 4;
  static const size_t num_blocks = basis_dimRange / block_size;
  using BlockMatrixType = XT::Common::BlockedFieldMatrix<RangeFieldType, num_blocks, block_size>;
  using LocalMatrixType = typename BlockMatrixType::BlockType;
  using BlockVectorType = XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>;
  using VectorType = BlockVectorType;
  using LocalVectorType = typename BlockVectorType::BlockType;
  using BasisValuesMatrixType = FieldVector<XT::LA::CommonDenseMatrix<RangeFieldType>, num_blocks>;
  using QuadraturePointsType =
      FieldVector<std::vector<BasisDomainType, boost::alignment::aligned_allocator<BasisDomainType, 64>>, num_blocks>;
  using QuadratureWeightsType =
      FieldVector<std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>, num_blocks>;
  using TemporaryVectorType = std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>;
  using TemporaryVectorsType = FieldVector<TemporaryVectorType, num_blocks>;
  using AlphaReturnType = std::pair<BlockVectorType, std::pair<DomainType, RangeFieldType>>;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , M_(XT::LA::CommonDenseMatrix<RangeFieldType>())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
  {
    XT::LA::eye_matrix(T_minus_one_);
    helper<basis_dimDomain>::calculate_plane_coefficients(basis_functions_);
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

  virtual int order(const XT::Common::Parameter& /*param*/) const override
  {
    return 1;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = std::make_unique<BlockVectorType>(get_alpha(u, *get_isotropic_alpha(u), true)->first);
    return evaluate_with_alpha(*alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const BlockVectorType& alpha) const
  {
    RangeReturnType ret(0.);
    auto& work_vecs = working_storage();
    calculate_scalar_products(alpha, M_, work_vecs);
    apply_exponential(work_vecs);
    for (size_t dd = 0; dd < basis_dimDomain; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = block_size * jj;
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto factor = work_vecs[jj][ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
          for (size_t ii = 0; ii < block_size; ++ii)
            ret[dd][offset + ii] += M_[jj].get_entry(ll, ii) * factor;
        } // ll
      } // jj
    } // dd
    return ret;
  } // void evaluate(...)

  virtual void jacobian(const DomainType& u,
                        DynamicDerivativeRangeType& result,
                        const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = std::make_unique<BlockVectorType>(get_alpha(u, *get_isotropic_alpha(u), true)->first);
    jacobian_with_alpha(*alpha, result);
  }

  virtual void jacobian_with_alpha(const BlockVectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    thread_local auto H = XT::Common::make_unique<BlockMatrixType>();
    calculate_hessian(alpha, M_, *H);
    for (size_t dd = 0; dd < basis_dimDomain; ++dd)
      row_jacobian(dd, M_, *H, result[dd], dd > 0);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType evaluate_kinetic_flux(const DomainType& u_i,
                                   const DomainType& u_j,
                                   const BasisDomainType& n_ij,
                                   const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = std::make_unique<BlockVectorType>(get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first);
    const auto alpha_j = std::make_unique<BlockVectorType>(get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first);
    evaluate_kinetic_flux_with_alphas(*alpha_i, *alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const BlockVectorType& alpha_i,
                                               const BlockVectorType& alpha_j,
                                               const BasisDomainType& n_ij,
                                               const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    thread_local FieldVector<TemporaryVectorsType, 2> work_vecs;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      work_vecs[0][jj].resize(quad_points_[jj].size());
      work_vecs[1][jj].resize(quad_points_[jj].size());
    }
    calculate_scalar_products(alpha_i, M_, work_vecs[0]);
    calculate_scalar_products(alpha_j, M_, work_vecs[1]);
    DomainType ret(0);
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
  } // DomainType evaluate_kinetic_flux(...)

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();

    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_iso_prime = std::make_unique<BlockVectorType>(basis_functions_.alpha_iso_prime());
    auto alpha_initial = std::make_unique<BlockVectorType>(*alpha_iso_prime);
    *alpha_initial *= -std::log(density);
    *alpha_initial += alpha_in;
    if (!(density > 0. || !(basis_functions_.min_density(u) > 0.)) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    auto u_prime = std::make_unique<const BlockVectorType>(u / density);

    // if value has already been calculated for these values, skip computation
    RangeFieldType tau_prime = std::min(
        tau_ / ((1 + std::sqrt(basis_dimRange) * u_prime->two_norm()) * density + std::sqrt(basis_dimRange) * tau_),
        tau_);

    // calculate moment vector for isotropic distribution
    auto u_iso = std::make_unique<const BlockVectorType>(basis_functions_.u_iso());

    // define further variables
    auto g_k = std::make_unique<BlockVectorType>();
    auto beta_out = std::make_unique<BlockVectorType>();
    auto v = std::make_unique<BlockVectorType>();
    thread_local auto T_k = XT::Common::make_unique<BlockMatrixType>();
    auto beta_in = std::make_unique<BlockVectorType>(*alpha_initial);

    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence) {
      // regularize u
      *v = *u_prime;
      if (r > 0.) {
        *beta_in = *get_isotropic_alpha(u);
        // calculate v = (1-r) u + r u_iso
        // use beta_out as storage for u_iso_in * r
        *v *= (1 - r);
        *beta_out = *u_iso;
        *beta_out *= r;
        *v += *beta_out;
      }
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
        if (!(density_tilde > 0.) || !(basis_functions_.min_density(*u_alpha_tilde) > 0.) || std::isinf(density_tilde))
          break;
        *alpha_prime = *alpha_iso_prime;
        *alpha_prime *= -std::log(density_tilde);
        *alpha_prime += *alpha_tilde;
        calculate_vector_integral(*alpha_prime, M_, M_, *u_alpha_prime);
        *u_eps_diff = *u_alpha_prime;
        *u_eps_diff *= -(1 - epsilon_gamma_);
        *u_eps_diff += *v;
        if (g_alpha_tilde->two_norm() < tau_prime
            && 1 - epsilon_gamma_ < std::exp(d_alpha_tilde->one_norm() + std::abs(std::log(density_tilde)))
            && helper<basis_dimDomain>::is_realizable(*u_eps_diff, basis_functions_)) {
          ret->first = *alpha_iso_prime;
          ret->first *= std::log(density);
          ret->first += *alpha_prime;
          ret->second.first = *v;
          ret->second.first *= density;
          ret->second.second = r;
          return ret;
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
    DUNE_THROW(MathError, "Failed to converge");

    return ret;
  }

private:
  // temporary vectors to store inner products and exponentials
  TemporaryVectorsType& working_storage() const
  {
    thread_local TemporaryVectorsType work_vecs;
    for (size_t jj = 0; jj < num_blocks; ++jj)
      work_vecs[jj].resize(quad_points_[jj].size());
    return work_vecs;
  }

  std::unique_ptr<BlockVectorType> get_isotropic_alpha(const DomainType& u) const
  {
    static const auto alpha_iso = basis_functions_.alpha_iso();
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    return std::make_unique<BlockVectorType>(alpha_iso + alpha_iso_prime * std::log(basis_functions_.density(u)));
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
#  if HAVE_MKL || HAVE_CBLAS
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
#  else
      LocalVectorType tmp_vec, tmp_vec2;
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        std::copy_n(M.get_ptr(ll), block_size, tmp_vec.begin());
        XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
        std::copy_n(tmp_vec2.begin(), block_size, M.get_ptr(ll));
      }
#  endif
    }
  }

  void apply_inverse_matrix(const BlockMatrixType& T_k, BasisValuesMatrixType& M) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      apply_inverse_matrix_block(jj, T_k.block(jj), M[jj]);
  }

  template <size_t domainDim = basis_dimDomain, class anything = void>
  struct helper
  {
#  if HAVE_QHULL
    static void calculate_plane_coefficients(const MomentBasis& basis_functions)
    {
      if (!basis_functions.plane_coefficients()[0].size())
        basis_functions.calculate_plane_coefficients();
    }

    static bool is_realizable(const BlockVectorType& u, const MomentBasis& basis_functions)
    {
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (const auto& coeff : basis_functions.plane_coefficients()[jj])
          if (!(u.block(jj) * coeff.first < coeff.second))
            return false;
      return true;
    }
#  else
    static void calculate_plane_coefficients(const MomentBasis& /*basis_functions*/)
    {
      DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
    }

    static bool is_realizable(const BlockVectorType& /*u*/, const MomentBasis& /*basis_functions*/)
    {
      DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
      return false;
    }
#  endif
  }; // class helper<...>

  template <class anything>
  struct helper<1, anything>
  {
    static void calculate_plane_coefficients(const MomentBasis& /*basis_functions*/) {}

    static bool is_realizable(const BlockVectorType& u, const MomentBasis& basis_functions)
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

  void row_jacobian(const size_t row,
                    const BasisValuesMatrixType& M,
                    BlockMatrixType& H,
                    DynamicRowDerivativeRangeType& ret,
                    bool L_calculated = false) const
  {
    assert(row < basis_dimDomain);
    calculate_J(M, ret, row);
    calculate_A_Binv(ret, H, L_calculated);
  } // void partial_u_col(...)

  // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
  static void calculate_A_Binv(DynamicRowDerivativeRangeType& A, BlockMatrixType& B, bool L_calculated = false)
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
  void calculate_J(const BasisValuesMatrixType& M, DynamicRowDerivativeRangeType& J_dd, const size_t dd) const
  {
    assert(dd < basis_dimDomain);
    J_dd.set_all_entries(0.);
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
          J_dd.set_entry(offset + mm, offset + nn, J_dd.get_entry(offset + nn, offset + mm));
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

  const MomentBasis& basis_functions_;
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
};
#endif // ENTROPY_FLUX_USE_PARTIAL_MOMENTS_SPECIALIZATION

#if ENTROPY_FLUX_USE_3D_HATFUNCTIONS_SPECIALIZATION
/**
 * Specialization of EntropyBasedFluxImplementation for 3D Hatfunctions
 */
template <class D, class R, size_t dimRange_or_refinements>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1>>
  : public XT::Functions::FunctionInterface<HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1>::dimRange,
                                            3,
                                            HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1>::dimRange,
                                            R>
{
public:
  using MomentBasis = HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1>;
  using BaseType = typename XT::Functions::
      FunctionInterface<MomentBasis::dimRange, MomentBasis::dimDomain, MomentBasis::dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;
  using BaseType::d;
  using BaseType::r;
  static const size_t basis_dimDomain = MomentBasis::dimDomain;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, 3>;
  using LocalMatrixType = XT::Common::FieldMatrix<RangeFieldType, 3, 3>;
  using BasisValuesMatrixType = std::vector<std::vector<LocalVectorType>>;
  using QuadraturePointsType = std::vector<std::vector<BasisDomainType>>;
  using QuadratureWeightsType = std::vector<std::vector<RangeFieldType>>;
#  if HAVE_EIGEN
  using SparseMatrixType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::eigen_sparse>::MatrixType;
  using VectorType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::eigen_sparse>::VectorType;
#  else
  using SparseMatrixType = typename XT::LA::Container<RangeFieldType, XT::LA::default_sparse_backend>::MatrixType;
  using VectorType = typename XT::LA::Container<RangeFieldType, XT::LA::default_sparse_backend>::VectorType;
#  endif
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
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
  {
    const auto& triangulation = basis_functions_.triangulation();
    const auto& vertices = triangulation.vertices();
    const auto& faces = triangulation.faces();
    assert(vertices.size() == basis_dimRange);
    // create pattern
    XT::LA::SparsityPatternDefault pattern(basis_dimRange);
    for (size_t vertex_index = 0; vertex_index < basis_dimRange; ++vertex_index) {
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
  } // constructor

  virtual int order(const XT::Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  VectorType get_isotropic_alpha(const DomainType& u) const
  {
    static const auto alpha_iso = basis_functions_.alpha_iso();
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    const auto ret_dynvector = alpha_iso + alpha_iso_prime * std::log(basis_functions_.density(u));
    VectorType ret(ret_dynvector.size());
    for (size_t ii = 0; ii < ret.size(); ++ii)
      ret[ii] = ret_dynvector[ii];
    return ret;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    LocalVectorType local_alpha, local_ret;
    const auto& triangulation = basis_functions_.triangulation();
    const auto& faces = triangulation.faces();
    for (size_t dd = 0; dd < basis_dimDomain; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t jj = 0; jj < faces.size(); ++jj) {
        local_ret *= 0.;
        const auto& face = faces[jj];
        const auto& vertices = face->vertices();
        for (size_t ii = 0; ii < 3; ++ii)
          local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          auto factor_ll = std::exp(local_alpha * basis_ll) * quad_points_[jj][ll][dd] * quad_weights_[jj][ll];
          for (size_t ii = 0; ii < 3; ++ii)
            local_ret[ii] += basis_ll[ii] * factor_ll;
        } // ll (quad points)
        for (size_t ii = 0; ii < 3; ++ii)
          ret[dd][vertices[ii]->index()] += local_ret[ii];
      } // jj (faces)
    } // dd
    return ret;
  } // void evaluate(...)

  virtual void jacobian(const DomainType& u,
                        DynamicDerivativeRangeType& result,
                        const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  virtual void jacobian_with_alpha(const VectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    thread_local SparseMatrixType H(basis_dimRange, basis_dimRange, pattern_, 0);
    thread_local SparseMatrixType J(basis_dimRange, basis_dimRange, pattern_, 0);
    calculate_hessian(alpha, M_, H);
    for (size_t dd = 0; dd < basis_dimDomain; ++dd) {
      calculate_J(M_, J, dd);
      calculate_J_Hinv(J, H, result[dd]);
    }
  } // ... jacobian(...)

  DomainType evaluate_kinetic_flux(const DomainType& u_i,
                                   const DomainType& u_j,
                                   const BasisDomainType& n_ij,
                                   const size_t dd) const
  {
    const auto alpha_i = get_alpha(u_i, get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const BasisDomainType& n_ij,
                                               const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    DomainType ret(0);
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
  } // DomainType evaluate_kinetic_flux(...)

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();

    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");

    VectorType u_prime(basis_dimRange, 0., 0);
    for (size_t ii = 0; ii < basis_dimRange; ++ii)
      u_prime.set_entry(ii, u[ii] / density);
    VectorType alpha_iso_prime(basis_dimRange, 0., 0);
    basis_functions_.alpha_iso_prime(alpha_iso_prime);

    // if value has already been calculated for these values, skip computation
    RangeFieldType tau_prime = std::min(
        tau_ / ((1 + std::sqrt(basis_dimRange) * u_prime.l2_norm()) * density + std::sqrt(basis_dimRange) * tau_),
        tau_);
    thread_local SparseMatrixType H(basis_dimRange, basis_dimRange, pattern_, 0);
    thread_local auto solver = XT::LA::make_solver(H);

    // calculate moment vector for isotropic distribution
    VectorType u_iso(basis_dimRange, 0., 0);
    basis_functions_.u_iso(u_iso);
    VectorType alpha_k = alpha_in - alpha_iso_prime * std::log(density);
    VectorType v(basis_dimRange, 0., 0), g_k(basis_dimRange, 0., 0), d_k(basis_dimRange, 0., 0),
        tmp_vec(basis_dimRange, 0., 0), alpha_prime(basis_dimRange);
    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence_) {
      // regularize u
      v = u_prime;
      if (r > 0) {
        alpha_k = get_isotropic_alpha(u);
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
        alpha_prime = alpha_iso_prime;
        alpha_prime *= -std::log(density_tilde);
        alpha_prime += alpha_tilde;
        auto& u_eps_diff = tmp_vec;
        calculate_u(alpha_prime, u_eps_diff); // store u_alpha_prime in u_eps_diff
        u_eps_diff *= -(1 - epsilon_gamma_);
        u_eps_diff += v;
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.l2_norm() < tau_prime && is_realizable(u_eps_diff)) {
          ret->first = alpha_iso_prime;
          ret->first *= std::log(density);
          ret->first += alpha_prime;
          auto v_ret_eig = v * density;
          DomainType v_ret;
          for (size_t ii = 0; ii < d; ++ii)
            v_ret[ii] = v_ret_eig[ii];
          ret->second = std::make_pair(v_ret, r);
          return ret;
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
    DUNE_THROW(MathError, "Failed to converge");
    return ret;
  } // ... get_alpha(...)

private:
  static bool is_realizable(const VectorType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  // temporary vectors to store inner products and exponentials
  std::vector<std::vector<RangeFieldType>>& get_work_vecs() const
  {
    thread_local std::vector<std::vector<RangeFieldType>> work_vecs;
    const auto& triangulation = basis_functions_.triangulation();
    const auto& faces = triangulation.faces();
    work_vecs.resize(faces.size());
    for (size_t jj = 0; jj < faces.size(); ++jj)
      work_vecs[jj].resize(quad_points_[jj].size());
    return work_vecs;
  }

private:
  // calculates ret = J H^{-1}. H is assumed to be symmetric positive definite, which gives ret^T = H^{-T} J^T =
  // H^{-1} J^T, so we just have to solve y = H^{-1} x for each row x of J
  void calculate_J_Hinv(SparseMatrixType& J, const SparseMatrixType& H, DynamicRowDerivativeRangeType& ret) const
  {
    thread_local VectorType solution(basis_dimRange, 0., 0), tmp_rhs(basis_dimRange, 0., 0);
#  if HAVE_EIGEN
    typedef ::Eigen::SparseMatrix<RangeFieldType, ::Eigen::ColMajor> ColMajorBackendType;
    ColMajorBackendType colmajor_copy(H.backend());
    colmajor_copy.makeCompressed();
    typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
    SolverType solver;
    solver.analyzePattern(colmajor_copy);
    solver.factorize(colmajor_copy);
#  else // HAVE_EIGEN
    auto solver = XT::LA::make_solver(H);
#  endif // HAVE_EIGEN
    for (size_t ii = 0; ii < basis_dimRange; ++ii) {
      // copy row to VectorType
      for (size_t kk = 0; kk < basis_dimRange; ++kk)
        tmp_rhs.set_entry(kk, J.get_entry(ii, kk));
        // solve
#  if HAVE_EIGEN
      solution.backend() = solver.solve(tmp_rhs.backend());
#  else // HAVE_EIGEN
      solver.apply(tmp_rhs, solution);
#  endif
      // copy result to C
      for (size_t kk = 0; kk < basis_dimRange; ++kk)
        ret.set_entry(ii, kk, solution.get_entry(kk));
    } // ii
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
    auto& work_vecs = get_work_vecs();
    for (size_t jj = 0; jj < faces.size(); ++jj) {
      const auto& face = faces[jj];
      const auto& vertices = face->vertices();
      local_u *= 0.;
      for (size_t ii = 0; ii < 3; ++ii)
        local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M_[jj][ll];
        work_vecs[jj][ll] = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          local_u[ii] += basis_ll[ii] * work_vecs[jj][ll];
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
    auto& work_vecs = get_work_vecs();
    for (size_t jj = 0; jj < faces.size(); ++jj) {
      H_local *= 0.;
      const auto& face = faces[jj];
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii)
        local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M[jj][ll];
        if (!use_work_vecs_results)
          work_vecs[jj][ll] = std::exp(local_alpha * basis_ll) * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          for (size_t kk = 0; kk < 3; ++kk)
            H_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * work_vecs[jj][ll];
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
    assert(dd < basis_dimDomain);
    J_dd *= 0.;
    LocalMatrixType J_local(0.);
    auto& work_vecs = get_work_vecs();
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
            J_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * work_vecs[jj][ll] * quad_points_[jj][ll][dd];
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        for (size_t kk = 0; kk < 3; ++kk)
          J_dd.add_to_entry(vertices[ii]->index(), vertices[kk]->index(), J_local[ii][kk]);
    } // jj (faces)
  } // void calculate_J(...)

  const MomentBasis& basis_functions_;
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
  XT::LA::SparsityPatternDefault pattern_;
};
#endif // ENTROPY_FLUX_USE_3D_HATFUNCTIONS_SPECIALIZATION

#if ENTROPY_FLUX_USE_1D_HATFUNCTIONS_SPECIALIZATION
#  if ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS
/**
 * Specialization of EntropyBasedFluxImplementation for 1D Hatfunctions (no change of basis, analytic integrals +
 * Taylor)
 */
template <class D, class R, size_t dimRange>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 1, R, dimRange, 1>>
  : public XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>
{
  using BaseType = typename XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using MomentBasis = HatFunctionMomentBasis<D, 1, R, dimRange, 1>;
  using BaseType::d;
  using BaseType::r;
  static const size_t basis_dimDomain = MomentBasis::dimDomain;
  static const size_t basis_dimRange = dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using VectorType = XT::Common::FieldVector<RangeFieldType, basis_dimRange>;
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon,
                                          const RangeFieldType taylor_tol = 0.1,
                                          const size_t max_taylor_order = 200)
    : basis_functions_(basis_functions)
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
  {}

  static bool is_realizable(const DomainType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  virtual int order(const XT::Common::Parameter& /*param*/) const override
  {
    return 1;
  }

  VectorType get_isotropic_alpha(const DomainType& u) const
  {
    static const auto alpha_iso = basis_functions_.alpha_iso();
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    return alpha_iso + alpha_iso_prime * std::log(basis_functions_.density(u));
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    // calculate < \mu m G_\alpha(u) >
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (std::abs(alpha[nn] - alpha[nn - 1]) > taylor_tol_) {
          ret[0][nn] +=
              2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                  * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
              + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                    * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
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
          ret[0][nn] += result * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn - 1]);
        }
      }
      if (nn < dimRange - 1) {
        if (std::abs(alpha[nn + 1] - alpha[nn]) > taylor_tol_) {
          ret[0][nn] +=
              -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                  * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
              + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                    * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
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
          ret[0][nn] += result * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn
    return ret;
  } // void evaluate_with_alpha(...)

  virtual void jacobian(const DomainType& u,
                        DynamicDerivativeRangeType& result,
                        const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  virtual void jacobian_with_alpha(const VectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    VectorType H_diag, J_diag;
    XT::Common::FieldVector<RangeFieldType, dimRange - 1> H_subdiag, J_subdiag;
    calculate_hessian(alpha, H_diag, H_subdiag);
    calculate_J(alpha, J_diag, J_subdiag);
    calculate_J_Hinv(result[0], J_diag, J_subdiag, H_diag, H_subdiag);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType evaluate_kinetic_flux(const DomainType& u_i,
                                   const DomainType& u_j,
                                   const BasisDomainType& n_ij,
                                   const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const BasisDomainType& n_ij,
                                               const size_t dd) const
  {
    assert(dd == 0);
    // calculate < \mu m G_\alpha(u) > * n_ij
    DomainType ret(0);
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
  } // DomainType evaluate_kinetic_flux(...)


  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();
    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    VectorType u_prime = u / density;
    VectorType alpha_initial = alpha_in - alpha_iso_prime * std::log(density);
    RangeFieldType tau_prime =
        std::min(tau_ / ((1 + std::sqrt(dimRange) * u_prime.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
    // The hessian H is always symmetric and tridiagonal, so we only need to store the diagonal and subdiagonal
    // elements
    VectorType H_diag;
    FieldVector<RangeFieldType, dimRange - 1> H_subdiag;

    // calculate moment vector for isotropic distribution
    VectorType u_iso = basis_functions_.u_iso();
    VectorType v;
    VectorType alpha_k = alpha_initial;
    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence_) {
      // regularize u
      v = u_prime;
      if (r > 0) {
        alpha_k = get_isotropic_alpha(u);
        VectorType r_times_u_iso(u_iso);
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
        VectorType g_k = calculate_gradient(alpha_k, v);
        // calculate Hessian H
        calculate_hessian(alpha_k, H_diag, H_subdiag);
        // calculate descent direction d_k;
        VectorType d_k(0), minus_g_k(g_k);
        minus_g_k *= -1;
        try {
          d_k = minus_g_k;
          XT::LA::solve_sym_tridiag_posdef(H_diag, H_subdiag, d_k);
        } catch (const Dune::MathError&) {
          if (r < r_max)
            break;
          else
            DUNE_THROW(Dune::MathError, "Failure to converge!");
        }

        const auto& alpha_tilde = alpha_k;
        const auto u_alpha_tilde = g_k + v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;
        const auto alpha_prime = alpha_tilde - alpha_iso_prime * std::log(density_tilde);
        const auto u_alpha_prime = calculate_u(alpha_prime);
        auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)) {
          ret->first = alpha_prime + alpha_iso_prime * std::log(density);
          ret->second = std::make_pair(v * density, r);
          return ret;
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
    DUNE_THROW(MathError, "Failed to converge");

    return ret;
  } // ... get_alpha(...)

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

private:
  RangeFieldType calculate_f(const VectorType& alpha_k, const VectorType& v) const
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

  VectorType calculate_u(const VectorType& alpha_k) const
  {
    VectorType u(0);
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
  } // VectorType calculate_u(...)

  VectorType calculate_gradient(const VectorType& alpha_k, const VectorType& v) const
  {
    return calculate_u(alpha_k) - v;
  }

  void calculate_hessian(const VectorType& alpha_k,
                         VectorType& diag,
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
  calculate_J(const VectorType& alpha_k, VectorType& diag, FieldVector<RangeFieldType, dimRange - 1>& subdiag) const
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
                    * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1])) / std::pow(alpha_k[nn - 1] - alpha_k[nn], 4);
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
  static void calculate_J_Hinv(DynamicRowDerivativeRangeType& ret,
                               const VectorType& J_diag,
                               const FieldVector<RangeFieldType, dimRange - 1>& J_subdiag,
                               VectorType& H_diag,
                               FieldVector<RangeFieldType, dimRange - 1>& H_subdiag)
  {
    // factorize H = LDL^T, where L is unit lower bidiagonal and D is diagonal
    // H_diag is overwritten by the diagonal elements of D
    // H_subdiag is overwritten by the subdiagonal elements of L
    XT::LA::tridiagonal_ldlt(H_diag, H_subdiag);

    // copy J to dense matrix
    ret.set_all_entries(0.);
    for (size_t ii = 0; ii < dimRange - 1; ++ii) {
      ret[ii][ii] = J_diag[ii];
      ret[ii + 1][ii] = J_subdiag[ii];
      ret[ii][ii + 1] = J_subdiag[ii];
    }
    ret[dimRange - 1][dimRange - 1] = J_diag[dimRange - 1];

    // Solve ret H = J which is equivalent to (as H and J are symmetric) to H ret^T = J;
    XT::LA::solve_tridiagonal_ldlt_factorized(H_diag, H_subdiag, ret);
    // transpose ret
    RangeFieldType* ret_ptr = &(ret[0][0]);
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < ii; ++jj)
        std::swap(ret_ptr[jj * dimRange + ii], ret_ptr[ii * dimRange + jj]);
  } // void calculate_J_Hinv(...)

  const MomentBasis& basis_functions_;
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
};

#  else // ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS

/**
 * Specialization of EntropyBasedFluxImplementation for 1D Hatfunctions (no change of basis, use structure)
 */
template <class D, class R, size_t dimRange>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 1, R, dimRange, 1>>
  : public XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>
{
  using BaseType = typename XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using MomentBasis = HatFunctionMomentBasis<D, 1, R, dimRange, 1>;
  using BaseType::d;
  using BaseType::r;
  static const size_t basis_dimDomain = MomentBasis::dimDomain;
  static const size_t basis_dimRange = dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using VectorType = XT::Common::FieldVector<RangeFieldType, basis_dimRange>;
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;
  static const size_t num_intervals = dimRange - 1;
  static const size_t block_size = 2;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;
  using BasisValuesMatrixType = FieldVector<std::vector<LocalVectorType>, num_intervals>;
  using QuadraturePointsType = FieldVector<std::vector<RangeFieldType>, num_intervals>;
  using QuadratureWeightsType = FieldVector<std::vector<RangeFieldType>, num_intervals>;

  explicit EntropyBasedFluxImplementation(
      const MomentBasis& basis_functions,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52))
    : basis_functions_(basis_functions)
    , grid_points_(basis_functions_.triangulation())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
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

  virtual int order(const XT::Common::Parameter& /*param*/) const override
  {
    return 1;
  }

  VectorType get_isotropic_alpha(const DomainType& u) const
  {
    static const auto alpha_iso = basis_functions_.alpha_iso();
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    return alpha_iso + alpha_iso_prime * std::log(basis_functions_.density(u));
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    LocalVectorType local_alpha;
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ii = 0; ii < 2; ++ii)
        local_alpha[ii] = alpha[jj + ii];
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M_[jj][ll];
        auto factor_ll = std::exp(local_alpha * basis_ll) * quad_points_[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 2; ++ii)
          ret[0][jj + ii] += basis_ll[ii] * factor_ll;
      } // ll (quad points)
    } // jj (intervals)
    return ret;
  } // void evaluate(...)

  virtual void jacobian(const DomainType& u,
                        DynamicDerivativeRangeType& result,
                        const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  virtual void jacobian_with_alpha(const VectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    VectorType H_diag, J_diag;
    FieldVector<RangeFieldType, dimRange - 1> H_subdiag, J_subdiag;
    calculate_hessian(alpha, M_, H_diag, H_subdiag);
    calculate_J(M_, J_diag, J_subdiag);
    calculate_J_Hinv(result[0], J_diag, J_subdiag, H_diag, H_subdiag);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType evaluate_kinetic_flux(const DomainType& u_i,
                                   const DomainType& u_j,
                                   const BasisDomainType& n_ij,
                                   const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const BasisDomainType& n_ij,
                                               const size_t dd) const
  {
    assert(dd == 0);
    // calculate < \mu m G_\alpha(u) > * n_ij
    DomainType ret(0);
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
  } // ... evaluate_kinetic_flux(...)

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();
    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    static const auto alpha_iso_prime = basis_functions_.alpha_iso_prime();
    VectorType u_prime = u / density;
    VectorType alpha_initial = alpha_in - alpha_iso_prime * std::log(density);
    RangeFieldType tau_prime =
        std::min(tau_ / ((1 + std::sqrt(dimRange) * u_prime.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
    // The hessian H is always symmetric and tridiagonal, so we only need to store the diagonal and subdiagonal
    // elements
    VectorType H_diag;
    FieldVector<RangeFieldType, dimRange - 1> H_subdiag;

    // calculate moment vector for isotropic distribution
    VectorType u_iso = basis_functions_.u_iso();
    VectorType v;
    VectorType alpha_k = alpha_initial;

    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& r : r_sequence_) {
      // regularize u
      v = u_prime;
      if (r > 0) {
        alpha_k = get_isotropic_alpha(u);
        VectorType r_times_u_iso(u_iso);
        r_times_u_iso *= r;
        v *= 1 - r;
        v += r_times_u_iso;
      }

      // calculate f_0
      RangeFieldType f_k = calculate_f(alpha_k, v);

      int pure_newton = 0;
      VectorType g_k, d_k, minus_g_k, u_alpha_prime;
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
          if (r < r_max)
            break;
          else
            DUNE_THROW(Dune::MathError, "Failure to converge!");
        }

        const auto& alpha_tilde = alpha_k;
        const auto u_alpha_tilde = g_k + v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;
        const auto alpha_prime = alpha_tilde - alpha_iso_prime * std::log(density_tilde);
        calculate_u(alpha_prime, u_alpha_prime);
        auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)) {
          ret->first = alpha_prime + alpha_iso_prime * std::log(density);
          ret->second = std::make_pair(v * density, r);
          return ret;
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
    DUNE_THROW(MathError, "Failed to converge");
    return ret;
  } // ... get_alpha(...)

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

private:
  static bool is_realizable(const DomainType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  FieldVector<std::vector<RangeFieldType>, num_intervals>& working_storage() const
  {
    thread_local FieldVector<std::vector<RangeFieldType>, num_intervals> work_vec;
    for (size_t jj = 0; jj < num_intervals; ++jj)
      work_vec[jj].resize(quad_points_[jj].size());
    return work_vec;
  }


  RangeFieldType calculate_f(const VectorType& alpha, const VectorType& v) const
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

  void calculate_u(const VectorType& alpha, VectorType& u) const
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

  void calculate_gradient(const VectorType& alpha, const VectorType& v, VectorType& g_k) const
  {
    calculate_u(alpha, g_k);
    g_k -= v;
  }

  void calculate_hessian(const VectorType& alpha,
                         const BasisValuesMatrixType& M,
                         VectorType& H_diag,
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
                   VectorType& J_diag,
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
  static void calculate_J_Hinv(DynamicRowDerivativeRangeType& ret,
                               const VectorType& J_diag,
                               const FieldVector<RangeFieldType, dimRange - 1>& J_subdiag,
                               VectorType& H_diag,
                               FieldVector<RangeFieldType, dimRange - 1>& H_subdiag)
  {
    // factorize H = LDL^T, where L is unit lower bidiagonal and D is diagonal
    // H_diag is overwritten by the diagonal elements of D
    // H_subdiag is overwritten by the subdiagonal elements of L
    XT::LA::tridiagonal_ldlt(H_diag, H_subdiag);

    // copy J to dense matrix
    ret.set_all_entries(0.);
    for (size_t ii = 0; ii < dimRange - 1; ++ii) {
      ret.set_entry(ii, ii, J_diag[ii]);
      ret.set_entry(ii + 1, ii, J_subdiag[ii]);
      ret.set_entry(ii, ii + 1, J_subdiag[ii]);
    }
    ret.set_entry(dimRange - 1, dimRange - 1, J_diag[dimRange - 1]);

    // Solve ret H = J which is equivalent to (as H and J are symmetric) to H ret^T = J;
    XT::LA::solve_tridiagonal_ldlt_factorized(H_diag, H_subdiag, ret);
    // transpose ret
    RangeFieldType* ret_ptr = &(ret[0][0]);
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < ii; ++jj)
        std::swap(ret_ptr[jj * dimRange + ii], ret_ptr[ii * dimRange + jj]);
  } // void calculate_J_Hinv(...)

private:
  const MomentBasis& basis_functions_;
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
};
#  endif // ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS
#endif // ENTROPY_FLUX_USE_1D_HATFUNCTIONS_SPECIALIZATION


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_ENTROPYFLUX_IMPLEMENTATIONS_HH
