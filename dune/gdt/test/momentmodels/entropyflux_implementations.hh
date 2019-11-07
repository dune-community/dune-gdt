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
#include <dune/xt/common/cblas.hh>
#include <dune/xt/common/mkl.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/algorithms/solve_sym_tridiag_posdef.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/container/pattern.hh>

#include <dune/xt/functions/interfaces/function.hh>

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/type_traits.hh>

#include "config.h"

#if HAVE_CLP
#  include <coin/ClpSimplex.hpp>
#endif // HAVE_CLP

namespace Dune {
namespace GDT {


template <class T>
std::enable_if_t<std::is_arithmetic<T>::value, T> superbee(const T first_slope, const T second_slope)
{
  return XT::Common::maxmod(XT::Common::minmod(first_slope, 2. * second_slope),
                            XT::Common::minmod(2. * first_slope, second_slope));
}


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


enum class SlopeLimiterType
{
  minmod,
  superbee,
  mc,
  no_slope
};


/**
 * Unspecialized implementation, should work with all bases
 */
template <class MomentBasisImp>
class EntropyBasedFluxImplementationUnspecializedBase
  : public XT::Functions::FunctionInterface<MomentBasisImp::dimRange,
                                            MomentBasisImp::dimFlux,
                                            MomentBasisImp::dimRange,
                                            typename MomentBasisImp::R>
{
  using BaseType = typename XT::Functions::FunctionInterface<MomentBasisImp::dimRange,
                                                             MomentBasisImp::dimFlux,
                                                             MomentBasisImp::dimRange,
                                                             typename MomentBasisImp::R>;
  using ThisType = EntropyBasedFluxImplementationUnspecializedBase;

public:
  using MomentBasis = MomentBasisImp;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using FluxDomainType = FieldVector<DomainFieldType, dimFlux>;
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
  using QuadraturePointsType = std::vector<BasisDomainType, boost::alignment::aligned_allocator<BasisDomainType, 64>>;
  using QuadratureWeightsType = std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>;
  static const EntropyType entropy = MomentBasis::entropy;

  explicit EntropyBasedFluxImplementationUnspecializedBase(const MomentBasis& basis_functions,
                                                           const RangeFieldType tau,
                                                           const bool disable_realizability_check,
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
    , disable_realizability_check_(disable_realizability_check)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , realizability_helper_(basis_functions_,
                            quad_points_,
                            disable_realizability_check_
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


  // ============================================================================================
  // ============================= FunctionInterface methods ====================================
  // ============================================================================================


  int order(const XT::Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const DomainType& alpha) const
  {
    RangeReturnType ret(0.);
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, M_, eta_ast_prime_vals);
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t ll = 0; ll < quad_weights_.size(); ++ll) {
        const auto factor = eta_ast_prime_vals[ll] * quad_weights_[ll] * quad_points_[ll][dd];
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
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  virtual void jacobian_with_alpha(const DomainType& alpha, DynamicDerivativeRangeType& result) const
  {
    thread_local auto H = XT::Common::make_unique<MatrixType>();
    calculate_hessian(alpha, M_, *H);
    for (size_t dd = 0; dd < dimFlux; ++dd)
      row_jacobian(dd, M_, *H, result[dd], dd > 0);
  }

  void row_jacobian(const size_t row,
                    const BasisValuesMatrixType& M,
                    MatrixType& H,
                    DynamicRowDerivativeRangeType& ret,
                    bool L_calculated = false) const
  {
    assert(row < dimFlux);
    calculate_J(M, ret, row);
    calculate_J_Hinv(ret, H, L_calculated);
  } // void partial_u_col(...)


  // ============================================================================================
  // ============ Evaluations of ansatz distribution, moments, hessian etc. =====================
  // ============================================================================================


  // Solves the minimum entropy optimization problem for u.
  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const DomainType& alpha_in, const bool regularize) const = 0;

  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // returns density rho = < eta_ast_prime(beta_in * b(v)) >
  RangeFieldType get_rho(const DomainType& beta_in, const BasisValuesMatrixType& M) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(beta_in, M, eta_ast_prime_vals);
    return std::inner_product(
        quad_weights_.begin(), quad_weights_.end(), eta_ast_prime_vals.begin(), RangeFieldType(0.));
  }

  // returns < eta_ast(beta_in * b(v)) >
  RangeFieldType get_eta_ast_integrated(const DomainType& beta_in, const BasisValuesMatrixType& M) const
  {
    auto& eta_ast_vals = working_storage();
    evaluate_eta_ast(beta_in, M, eta_ast_vals);
    return std::inner_product(quad_weights_.begin(), quad_weights_.end(), eta_ast_vals.begin(), RangeFieldType(0.));
  }

  // returns < b \eta_{\ast}^{\prime}(\alpha^T b(v)) >
  DomainType get_u(const DomainType& alpha) const
  {
    DomainType ret;
    calculate_u(alpha, M_, ret);
    return ret;
  }

  DomainType get_u(const QuadratureWeightsType& eta_ast_prime_vals) const
  {
    DomainType ret(0.);
    const size_t num_quad_points = quad_weights_.size();
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto factor_ll = eta_ast_prime_vals[ll] * quad_weights_[ll];
      const auto* basis_ll = &(M_.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < basis_dimRange; ++ii)
        ret[ii] += basis_ll[ii] * factor_ll;
    } // ll
    return ret;
  }

  // calculate ret = < b eta_ast_prime(beta_in * b) >
  void calculate_u(const DomainType& beta_in,
                   const BasisValuesMatrixType& M,
                   DomainType& ret,
                   bool same_beta = false,
                   bool only_first_component = false) const
  {
    auto& eta_ast_prime_vals = working_storage();
    if (!same_beta)
      evaluate_eta_ast_prime(beta_in, M, eta_ast_prime_vals);
    std::fill(ret.begin(), ret.end(), 0.);
    const size_t num_quad_points = quad_weights_.size();
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto factor_ll = eta_ast_prime_vals[ll] * quad_weights_[ll];
      const auto* basis_ll = &(M.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < (only_first_component ? 1 : basis_dimRange); ++ii)
        ret[ii] += basis_ll[ii] * factor_ll;
    } // ll
  }

  void calculate_hessian(const DomainType& alpha,
                         const BasisValuesMatrixType& M,
                         MatrixType& H,
                         const bool use_stored_data = false) const
  {
    auto& eta_ast_twoprime_values = working_storage();
    if (!use_stored_data)
      evaluate_eta_ast_twoprime(alpha, M, eta_ast_twoprime_values);
    calculate_hessian(eta_ast_twoprime_values, M, H);
  } // void calculate_hessian(...)

  void calculate_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                         const BasisValuesMatrixType& M,
                         MatrixType& H) const
  {
    std::fill(H.begin(), H.end(), 0.);
    const size_t num_quad_points = quad_weights_.size();
    // matrix is symmetric, we only use lower triangular part
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      auto factor_ll = eta_ast_twoprime_vals[ll] * quad_weights_[ll];
      const auto* basis_ll = &(M.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < basis_dimRange; ++ii) {
        auto* H_row = &(H[ii][0]);
        const auto factor_ll_ii = basis_ll[ii] * factor_ll;
        for (size_t kk = 0; kk <= ii; ++kk) {
          H_row[kk] += basis_ll[kk] * factor_ll_ii;
        } // kk
      } // ii
    } // ll
  } // void calculate_hessian(...)

  void apply_inverse_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                             const DomainType& u,
                             DomainType& Hinv_u) const
  {
    thread_local auto H = XT::Common::make_unique<MatrixType>();
    calculate_hessian(eta_ast_twoprime_vals, M_, *H);
    XT::LA::cholesky(*H);
    thread_local DomainType tmp_vec;
    XT::LA::solve_lower_triangular(*H, tmp_vec, u);
    XT::LA::solve_lower_triangular_transposed(*H, Hinv_u, tmp_vec);
  }

  // J = df/dalpha is the derivative of the flux with respect to alpha.
  // As F = (f_1, f_2, f_3) is matrix-valued
  // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
  // vector-valued),
  // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
  // this function returns the dd-th matrix df_dd/dalpha of J
  // assumes eta_ast_twoprime_values already contains the needed \eta_{\ast}^{\prime \prime} (\alpha * m) values
  void calculate_J(const BasisValuesMatrixType& M, DynamicRowDerivativeRangeType& J_dd, const size_t dd) const
  {
    assert(dd < dimFlux);
    const auto& eta_ast_twoprime_values = working_storage();
    J_dd.set_all_entries(0.);
    const size_t num_quad_points = quad_points_.size();
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto factor_ll = eta_ast_twoprime_values[ll] * quad_points_[ll][dd] * quad_weights_[ll];
      const auto* basis_ll = &(M.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < basis_dimRange; ++ii) {
        const auto factor_ll_ii = factor_ll * basis_ll[ii];
        if (!XT::Common::is_zero(factor_ll_ii)) {
          for (size_t kk = 0; kk <= ii; ++kk)
            J_dd.unsafe_add_to_entry(ii, kk, basis_ll[kk] * factor_ll_ii);
        }
      } // ii
    } // ll
    // symmetric update for upper triangular part of J
    for (size_t mm = 0; mm < basis_dimRange; ++mm)
      for (size_t nn = mm + 1; nn < basis_dimRange; ++nn)
        J_dd.set_entry(mm, nn, J_dd.get_entry(nn, mm));
  } // void calculate_J(...)

  // calculates J = J H^{-1}. H is assumed to be symmetric positive definite.
  static void calculate_J_Hinv(DynamicRowDerivativeRangeType& A, MatrixType& B, bool L_calculated = false)
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
  } // void calculate_J_Hinv(...)


  // ============================================================================================
  // ============================= Entropy evaluations ==========================================
  // ============================================================================================


  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast(ret);
  }

  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t ll = 0; ll < ret.size(); ++ll)
        ret[ll] = -std::log(1 - ret[ll]);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_prime(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast_prime(ret);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast_prime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t ll = 0; ll < ret.size(); ++ll)
        ret[ll] /= (1 - ret[ll]);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void
  evaluate_eta_ast_twoprime(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast_twoprime(ret);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already
  // contains exp(alpha^T b(v_i))
  void evaluate_eta_ast_twoprime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t ll = 0; ll < ret.size(); ++ll)
        ret[ll] /= std::pow(1 - ret[ll], 2);
  }

  // stores evaluations of exp(alpha^T b(v_i)) for all quadrature points v_i
  void store_exp_evaluations(QuadratureWeightsType& exp_evaluations, const DomainType& alpha) const
  {
    exp_evaluations.resize(quad_points_.size());
    this->calculate_scalar_products(alpha, M_, exp_evaluations);
    this->apply_exponential(exp_evaluations);
  }

  void store_eta_ast_prime_vals(const QuadratureWeightsType& exp_evaluations, QuadratureWeightsType& eta_ast_prime_vals)
  {
    eta_ast_prime_vals = exp_evaluations;
    evaluate_eta_ast_prime(eta_ast_prime_vals);
  }

  void store_eta_ast_twoprime_vals(const QuadratureWeightsType& exp_evaluations,
                                   QuadratureWeightsType& eta_ast_twoprime_vals)
  {
    eta_ast_twoprime_vals = exp_evaluations;
    evaluate_eta_ast_twoprime(eta_ast_twoprime_vals);
  }

  // stores evaluations of a given boundary distribution psi(v) at all quadrature points v_i
  void store_boundary_distribution_evaluations(
      QuadratureWeightsType& boundary_distribution_evaluations,
      const std::function<RangeFieldType(const FluxDomainType&)>& boundary_distribution) const
  {
    boundary_distribution_evaluations.resize(quad_points_.size());
    for (size_t ll = 0; ll < quad_points_.size(); ++ll)
      boundary_distribution_evaluations[ll] = boundary_distribution(quad_points_[ll]);
  }


  // ============================================================================================
  // =============================== Kinetic fluxes =============================================
  // ============================================================================================


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, \psi is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType
  evaluate_kinetic_flux(const DomainType& u_i, const DomainType& u_j, const FluxDomainType& n_ij, const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const DomainType& alpha_i,
                                               const DomainType& alpha_j,
                                               const FluxDomainType& n_ij,
                                               const size_t dd) const

  {
    thread_local FieldVector<QuadratureWeightsType, 2> eta_ast_prime_vals;
    eta_ast_prime_vals[0].resize(quad_points_.size());
    eta_ast_prime_vals[1].resize(quad_points_.size());
    evaluate_eta_ast_prime(alpha_i, M_, eta_ast_prime_vals[0]);
    evaluate_eta_ast_prime(alpha_j, M_, eta_ast_prime_vals[1]);
    DomainType ret(0);
    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      const auto position = quad_points_[ll][dd];
      RangeFieldType factor = position * n_ij[dd] > 0. ? eta_ast_prime_vals[0][ll] : eta_ast_prime_vals[1][ll];
      factor *= quad_weights_[ll] * position;
      const auto* basis_ll = &(M_.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < basis_dimRange; ++ii)
        ret[ii] += basis_ll[ii] * factor;
    } // ll
    ret *= n_ij[dd];
    return ret;
  } // DomainType evaluate_kinetic_flux_with_alphas(...)

  // Calculates left and right kinetic flux with reconstructed densities. Ansatz distribution values contains
  // evaluations of the ansatz distribution at each quadrature point for a stencil of three entities. The distributions
  // are reconstructed pointwise for each quadrature point and the resulting (part of) the kinetic flux is <
  // psi_reconstr * b * v>_{+/-}.
  template <SlopeLimiterType slope_type, class FluxesMapType>
  void calculate_reconstructed_fluxes(const FieldVector<const QuadratureWeightsType*, 3>& ansatz_distribution_values,
                                      FluxesMapType& flux_values,
                                      const size_t dd) const
  {
    // get left and right reconstructed values for each quadrature point v_i
    thread_local XT::Common::FieldVector<QuadratureWeightsType, 2> reconstructed_values(
        QuadratureWeightsType(quad_points_.size()));
    auto& vals_left = reconstructed_values[0];
    auto& vals_right = reconstructed_values[1];
    if (slope_type == SlopeLimiterType::no_slope) {
      for (size_t ll = 0; ll < quad_points_.size(); ++ll)
        vals_left[ll] = vals_right[ll] = (*ansatz_distribution_values[1])[ll];
    } else {
      const auto slope_func =
          (slope_type == SlopeLimiterType::minmod) ? XT::Common::minmod<RangeFieldType> : superbee<RangeFieldType>;
      for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
        const auto slope = slope_func((*ansatz_distribution_values[1])[ll] - (*ansatz_distribution_values[0])[ll],
                                      (*ansatz_distribution_values[2])[ll] - (*ansatz_distribution_values[1])[ll]);
        vals_left[ll] = (*ansatz_distribution_values[1])[ll] - 0.5 * slope;
        vals_right[ll] = (*ansatz_distribution_values[1])[ll] + 0.5 * slope;
      } // ll
    }

    BasisDomainType coord(0.5);
    coord[dd] = 0;
    auto& left_flux_value = flux_values[coord];
    coord[dd] = 1;
    auto& right_flux_value = flux_values[coord];
    right_flux_value = left_flux_value = DomainType(0.);

    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      const auto position = quad_points_[ll][dd];
      RangeFieldType factor = position > 0. ? vals_right[ll] : vals_left[ll];
      factor *= quad_weights_[ll] * position;
      auto& val = position > 0. ? right_flux_value : left_flux_value;
      const auto* basis_ll = &(M_.get_entry_ref(ll, 0.));
      for (size_t ii = 0; ii < basis_dimRange; ++ii)
        val[ii] += basis_ll[ii] * factor;
    } // ll
  } // void calculate_reconstructed_fluxes(...)


  // ============================================================================================
  // ================================== Helper functions ========================================
  // ============================================================================================


  // get permutation instead of sorting directly to be able to sort two vectors the same way
  // see
  // https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  template <class T, class Alloc, class Compare>
  static std::vector<std::size_t> get_sort_permutation(const std::vector<T, Alloc>& vec, const Compare& compare)
  {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
    return p;
  }

  template <class T, class Alloc>
  static void apply_permutation_in_place(std::vector<T, Alloc>& vec, const std::vector<std::size_t>& p)
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
  static void join_duplicate_quadpoints(QuadraturePointsType& quad_points, QuadratureWeightsType& quad_weights)
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

  // temporary vectors to store inner products and exponentials
  QuadratureWeightsType& working_storage() const
  {
    thread_local QuadratureWeightsType work_vec;
    work_vec.resize(quad_points_.size());
    return work_vec;
  }

  bool all_positive(const QuadratureWeightsType& vals) const
  {
    for (size_t ll = 0; ll < quad_points_.size(); ++ll) {
      const auto val = vals[ll];
      if (val < 0. || std::isinf(val) || std::isnan(val))
        return false;
    }
    return true;
  }

  void copy_transposed(const MatrixType& T_k, MatrixType& T_k_trans) const
  {
    for (size_t ii = 0; ii < basis_dimRange; ++ii)
      for (size_t jj = 0; jj <= ii; ++jj)
        T_k_trans[jj][ii] = T_k[ii][jj];
  }

  // calculates alpha^T b(v_i) for all quadrature points v_i
  void calculate_scalar_products(const DomainType& beta_in,
                                 const BasisValuesMatrixType& M,
                                 QuadratureWeightsType& scalar_products) const
  {
#if HAVE_MKL
    XT::Common::Cblas::dgemv(XT::Common::Cblas::row_major(),
                             XT::Common::Cblas::no_trans(),
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
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto* basis_ll = &(M.get_entry_ref(ll, 0.));
      scalar_products[ll] = std::inner_product(beta_in.begin(), beta_in.end(), basis_ll, 0.);
    }
#endif
  }

  // calculates exp(val) for all vals in values
  void apply_exponential(QuadratureWeightsType& values) const
  {
    assert(values.size() < std::numeric_limits<int>::max());
    XT::Common::Mkl::exp(static_cast<int>(values.size()), values.data(), values.data());
  }

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const RangeFieldType density) const
  {
    return std::make_unique<VectorType>(basis_functions_.alpha_iso(density));
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const DomainType& u) const
  {
    return get_isotropic_alpha(basis_functions_.density(u));
  }

#if HAVE_CLP
  template <class BasisFuncImp = MomentBasis, bool anything = true>
  struct RealizabilityHelper
  {
    static_assert(std::is_same<BasisFuncImp, MomentBasis>::value, "BasisFuncImp has to be MomentBasis!");

    RealizabilityHelper(const MomentBasis& basis_functions,
                        const QuadraturePointsType& quad_points,
                        const bool /*disable_realizability_check*/,
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
      const auto phi = u / density;
      setup_linear_program(reinitialize);
      auto& lp = **lp_;
      constexpr int num_rows = static_cast<int>(basis_dimRange);
      // set rhs (equality constraints, so set both bounds equal
      for (int ii = 0; ii < num_rows; ++ii) {
        size_t uii = static_cast<size_t>(ii);
        lp.setRowLower(ii, phi[uii]);
        lp.setRowUpper(ii, phi[uii]);
      }
      // set maximal wall time. If this is not set, in rare cases the primal method never returns
      lp.setMaximumWallSeconds(60);
      // Now check solvability
      lp.primal();
      return lp.primalFeasible();
    }

    const MomentBasis& basis_functions_;
    const QuadraturePointsType& quad_points_;
    XT::Common::PerThreadValue<std::unique_ptr<ClpSimplex>>& lp_;
  }; // struct RealizabilityHelper<...>
#else // HAVE_CLP
  template <class BasisFuncImp = MomentBasis, bool anything = true>
  struct RealizabilityHelper
  {
    RealizabilityHelper(const MomentBasis& /*basis_functions*/,
                        const QuadraturePointsType& /*quad_points*/,
                        const bool disable_realizability_check)
    {
      if (!disable_realizability_check)
        std::cerr << "Warning: You are missing Clp, realizability stopping condition will not be checked!" << std::endl;
    }

    bool is_realizable(const DomainType& /*u*/, const bool /*reinitialize*/) const
    {
      return true;
    }
  }; // struct RealizabilityHelper<...>
#endif // HAVE_CLP

  // specialization for hatfunctions
  template <size_t dimRange_or_refinements, bool anything>
  struct RealizabilityHelper<
      HatFunctionMomentBasis<DomainFieldType, dimFlux, RangeFieldType, dimRange_or_refinements, 1, dimFlux>,
      anything>
  {
    RealizabilityHelper(const MomentBasis& /*basis_functions*/,
                        const std::vector<BasisDomainType>& /*quad_points*/,
                        const bool /*disable_realizability_check*/
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

  // For each basis evaluation b, calculates T_k^{-1} b. As the basis evaluations are the rows of M, we want to
  // calculate (T_k^{-1} M^T)^T = M T_k^{-T}
  void apply_inverse_matrix(const MatrixType& T_k, BasisValuesMatrixType& M) const
  {
#if HAVE_MKL
    // Calculate the transpose here first as this is much faster than passing the matrix to dtrsm and using CblasTrans
    thread_local auto T_k_trans = std::make_unique<MatrixType>(0.);
    copy_transposed(T_k, *T_k_trans);
    assert(quad_points_.size() < std::numeric_limits<int>::max());
    XT::Common::Cblas::dtrsm(XT::Common::Cblas::row_major(),
                             XT::Common::Cblas::right(),
                             XT::Common::Cblas::upper(),
                             XT::Common::Cblas::no_trans(),
                             XT::Common::Cblas::non_unit(),
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
      auto* M_row = &(M.get_entry_ref(ll, 0.));
      std::copy_n(M_row, basis_dimRange, tmp_vec.begin());
      XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
      std::copy_n(tmp_vec2.begin(), basis_dimRange, M_row);
    }
#endif
  }

  const MomentBasis& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const bool disable_realizability_check_;
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
  using BaseType::basis_dimRange;
  using BaseType::dimFlux;
  using BaseType::entropy;
  using typename BaseType::AlphaReturnType;
  using typename BaseType::BasisDomainType;
  using typename BaseType::BasisValuesMatrixType;
  using typename BaseType::DomainType;
  using typename BaseType::FluxDomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::MomentBasis;
  using typename BaseType::QuadratureWeightsType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::VectorType;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const bool disable_realizability_check,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : BaseType(
          basis_functions, tau, disable_realizability_check, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon)
    , T_minus_one_(std::make_unique<MatrixType>())
  {
    XT::LA::eye_matrix(*T_minus_one_);
  }

  using BaseType::get_alpha;

  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const DomainType& alpha_in, const bool regularize) const override final
  {
    auto ret = std::make_unique<AlphaReturnType>();

    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_one = basis_functions_.alpha_one();
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");

    constexpr bool rescale = (entropy == EntropyType::MaxwellBoltzmann);

    VectorType phi = rescale ? u / density : u;
    VectorType alpha_initial = rescale ? alpha_in - alpha_one * std::log(density) : alpha_in;
    VectorType beta_in = alpha_initial;
    VectorType v, u_eps_diff, g_k, beta_out;
    RangeFieldType first_error_cond, second_error_cond, tau_prime;

    const auto u_iso = rescale ? basis_functions_.u_iso() : basis_functions_.u_iso() * density;
    const RangeFieldType dim_factor = is_full_moment_basis<MomentBasis>::value ? 1. : std::sqrt(dimFlux);
    tau_prime =
        rescale ? std::min(tau_ / ((1 + dim_factor * phi.two_norm()) * density + dim_factor * tau_), tau_) : tau_;

    thread_local auto T_k = XT::Common::make_unique<MatrixType>();

    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& rr : r_sequence) {
      // regularize u
      v = phi;
      if (rr > 0) {
        beta_in = *get_isotropic_alpha(v);
        VectorType r_times_u_iso = u_iso;
        r_times_u_iso *= rr;
        v *= 1 - rr;
        v += r_times_u_iso;
      }
      *T_k = *T_minus_one_;
      // calculate T_k u
      VectorType v_k = v;
      // calculate values of basis p = S_k m
      thread_local BasisValuesMatrixType P_k(M_.backend(), false, 0., 0);
      std::copy_n(M_.data(), M_.rows() * M_.cols(), P_k.data());
      // calculate f_0
      RangeFieldType f_k = get_eta_ast_integrated(beta_in, P_k);
      f_k -= beta_in * v_k;

      thread_local auto H = XT::Common::make_unique<MatrixType>(0.);

      int backtracking_failed = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase rr if too many iterations are used or cholesky decomposition fails
        if (kk > k_0_ && rr < r_max)
          break;
        try {
          change_basis(beta_in, v_k, P_k, *T_k, g_k, beta_out, *H);
        } catch (const Dune::MathError&) {
          if (rr < r_max)
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

        // Calculate values for stopping criteria (in original basis).
        VectorType alpha_tilde, d_alpha_tilde;
        XT::LA::solve_lower_triangular_transposed(*T_k, alpha_tilde, beta_out);
        XT::LA::solve_lower_triangular_transposed(*T_k, d_alpha_tilde, d_k);
        VectorType u_alpha_tilde;
        calculate_u(alpha_tilde, M_, u_alpha_tilde);
        VectorType g_alpha_tilde = u_alpha_tilde - v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;

        // if rescale == true, we ensure that the ansatz density exactly matches the density of u
        const auto alpha_prime = rescale ? alpha_tilde - alpha_one * std::log(density_tilde) : alpha_tilde;
        VectorType u_alpha_prime;
        if (rescale)
          calculate_u(alpha_prime, M_, u_alpha_prime);
        else
          u_alpha_prime = u_alpha_tilde;

        // calculate stopping criteria
        u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        first_error_cond = g_alpha_tilde.two_norm();
        second_error_cond = std::exp(
            -(rescale ? d_alpha_tilde.one_norm() + std::abs(std::log(density_tilde)) : d_alpha_tilde.one_norm()));
        // working_storage contains eta_ast_prime evaluations due to the calculate_u call above
        const auto& eta_ast_prime_vals = working_storage();
        if (first_error_cond < tau_prime && 1 - epsilon_gamma_ < second_error_cond
            && (entropy == EntropyType::MaxwellBoltzmann || all_positive(eta_ast_prime_vals))
            && (disable_realizability_check_
                || realizability_helper_.is_realizable(u_eps_diff, kk == static_cast<size_t>(0.8 * k_0_)))) {
          ret->first = rescale ? alpha_prime + alpha_one * std::log(density) : alpha_prime;
          ret->second = std::make_pair(rescale ? v * density : v, rr);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          beta_in = beta_out;
          // backtracking line search
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm()) {
            VectorType beta_new = d_k;
            beta_new *= zeta_k;
            beta_new += beta_out;
            RangeFieldType f = get_eta_ast_integrated(beta_new, P_k);
            f -= beta_new * v_k;
            if (backtracking_failed >= 2 || XT::Common::FloatCmp::le(f, f_k + xi_ * zeta_k * (g_k * d_k))) {
              beta_in = beta_new;
              f_k = f;
              backtracking_failed = 0;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          if (zeta_k <= epsilon_ * beta_out.two_norm() / d_k.two_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    const std::string err_msg = "Failed to converge for " + XT::Common::to_string(u) + " with density "
                                + XT::Common::to_string(density) + " and multiplier " + XT::Common::to_string(beta_in)
                                + " due to too many iterations! Last u_eps_diff = " + XT::Common::to_string(u_eps_diff)
                                + ", first_error_cond = " + XT::Common::to_string(first_error_cond)
                                + ", second_error_cond = " + XT::Common::to_string(second_error_cond)
                                + ", tau_prime = " + XT::Common::to_string(tau_prime);
    DUNE_THROW(MathError, err_msg);

    return ret;
  }

  using BaseType::all_positive;
  using BaseType::apply_inverse_matrix;
  using BaseType::calculate_hessian;
  using BaseType::calculate_u;
  using BaseType::get_eta_ast_integrated;
  using BaseType::get_isotropic_alpha;
  using BaseType::working_storage;

  void change_basis(const DomainType& beta_in,
                    DomainType& v_k,
                    BasisValuesMatrixType& P_k,
                    MatrixType& T_k,
                    DomainType& g_k,
                    DomainType& beta_out,
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
    calculate_u(beta_out, P_k, g_k, false);
    g_k -= v_k;
  } // void change_basis(...)

  using BaseType::basis_functions_;
  using BaseType::chi_;
  using BaseType::disable_realizability_check_;
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
  using BaseType::dimFlux;
  using BaseType::entropy;
  using typename BaseType::AlphaReturnType;
  using typename BaseType::BasisValuesMatrixType;
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::VectorType;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const bool disable_realizability_check,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : BaseType(
          basis_functions, tau, disable_realizability_check, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon)
  {}

  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  virtual std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const DomainType& alpha_in, const bool regularize) const override final
  {
    auto ret = std::make_unique<AlphaReturnType>();
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_one = basis_functions_.alpha_one();
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");

    constexpr bool rescale = (entropy == EntropyType::MaxwellBoltzmann);

    VectorType phi = rescale ? u / density : u;
    VectorType alpha_initial = rescale ? alpha_in - alpha_one * std::log(density) : alpha_in;
    VectorType v, g_k, d_k, tmp_vec, alpha_prime;
    RangeFieldType first_error_cond, second_error_cond, tau_prime;
    auto u_iso = rescale ? basis_functions_.u_iso() : basis_functions_.u_iso() * density;
    const RangeFieldType dim_factor = is_full_moment_basis<MomentBasis>::value ? 1. : std::sqrt(dimFlux);
    tau_prime =
        rescale ? std::min(tau_ / ((1 + dim_factor * phi.two_norm()) * density + dim_factor * tau_), tau_) : tau_;
    VectorType alpha_k = alpha_initial;
    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& rr : r_sequence) {
      // regularize u
      v = phi;
      if (rr > 0) {
        alpha_k = *get_isotropic_alpha(v);
        VectorType r_times_u_iso = u_iso;
        r_times_u_iso *= rr;
        v *= 1 - rr;
        v += r_times_u_iso;
      }
      // calculate T_k u
      VectorType v_k = v;
      // calculate f_0
      RangeFieldType f_k = get_eta_ast_integrated(alpha_k, M_);
      f_k -= alpha_k * v_k;

      thread_local auto H = XT::Common::make_unique<MatrixType>(0.);

      int backtracking_failed = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase rr if too many iterations are used
        if (kk > k_0_ && rr < r_max)
          break;
        // calculate gradient g
        calculate_u(alpha_k, M_, g_k);
        g_k -= v_k;
        // calculate Hessian H
        calculate_hessian(alpha_k, M_, *H, entropy == EntropyType::MaxwellBoltzmann);
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
          if (rr < r_max)
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
        auto& u_eps_diff = tmp_vec;
        alpha_prime = alpha_tilde;
        if (rescale) {
          alpha_prime -= alpha_one * std::log(density_tilde);
          calculate_u(alpha_prime, M_, u_eps_diff);
        } else {
          u_eps_diff = u_alpha_tilde;
        }
        u_eps_diff *= -(1 - epsilon_gamma_);
        u_eps_diff += v;

        first_error_cond = g_k.two_norm();
        second_error_cond = std::exp(-(rescale ? d_k.one_norm() + std::abs(std::log(density_tilde)) : d_k.one_norm()));
        auto& eta_ast_prime_vals = working_storage();
        // if rescale is true, working storage already contains the eta_ast_prime evaluations due to the call to
        // calculate_u above
        if (!rescale)
          evaluate_eta_ast_prime(alpha_prime, M_, eta_ast_prime_vals);
        if (first_error_cond < tau_prime && 1 - epsilon_gamma_ < second_error_cond
            && (entropy == EntropyType::MaxwellBoltzmann || all_positive(eta_ast_prime_vals))
            && (disable_realizability_check_
                || realizability_helper_.is_realizable(u_eps_diff, kk == static_cast<size_t>(0.8 * k_0_)))) {
          ret->first = rescale ? alpha_prime + alpha_one * std::log(density) : alpha_prime;
          ret->second = std::make_pair(rescale ? v * density : v, rr);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          auto& alpha_new = tmp_vec;
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
            // calculate alpha_new = alpha_k + zeta_k d_k
            alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;
            // calculate f(alpha_new)
            RangeFieldType f_new = get_eta_ast_integrated(alpha_new, M_);
            f_new -= alpha_new * v_k;
            if (backtracking_failed >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              backtracking_failed = 0;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
          if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    const std::string err_msg =
        "Failed to converge for " + XT::Common::to_string(u) + " with density " + XT::Common::to_string(density);
    DUNE_THROW(MathError, err_msg);
    return ret;
  }

  using BaseType::all_positive;
  using BaseType::calculate_hessian;
  using BaseType::calculate_u;
  using BaseType::evaluate_eta_ast_prime;
  using BaseType::get_eta_ast_integrated;
  using BaseType::get_isotropic_alpha;
  using BaseType::working_storage;

  using BaseType::basis_functions_;
  using BaseType::chi_;
  using BaseType::disable_realizability_check_;
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
template <class D, size_t d, class R, size_t dimRange_or_refinements, size_t fluxDim, EntropyType entropy>
class EntropyBasedFluxImplementation<PartialMomentBasis<D, d, R, dimRange_or_refinements, 1, fluxDim, 1, entropy>>
  : public XT::Functions::FunctionInterface<
        PartialMomentBasis<D, d, R, dimRange_or_refinements, 1, fluxDim, 1, entropy>::dimRange,
        fluxDim,
        PartialMomentBasis<D, d, R, dimRange_or_refinements, 1, fluxDim, 1, entropy>::dimRange,
        R>
{
public:
  using MomentBasis = PartialMomentBasis<D, d, R, dimRange_or_refinements, 1, fluxDim, 1, entropy>;
  using BaseType =
      typename XT::Functions::FunctionInterface<MomentBasis::dimRange, MomentBasis::dimFlux, MomentBasis::dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using FluxDomainType = FieldVector<DomainFieldType, dimFlux>;
  static const size_t block_size = (dimFlux == 1) ? 2 : 4;
  static const size_t num_blocks = basis_dimRange / block_size;
  using BlockMatrixType = XT::Common::BlockedFieldMatrix<RangeFieldType, num_blocks, block_size>;
  using LocalMatrixType = typename BlockMatrixType::BlockType;
  using BlockVectorType = XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>;
  using VectorType = BlockVectorType;
  using LocalVectorType = typename BlockVectorType::BlockType;
  using BasisValuesMatrixType = std::vector<XT::LA::CommonDenseMatrix<RangeFieldType>>;
  using QuadraturePointsType =
      std::vector<std::vector<BasisDomainType, boost::alignment::aligned_allocator<BasisDomainType, 64>>>;
  using BlockQuadratureWeightsType =
      std::vector<RangeFieldType, boost::alignment::aligned_allocator<RangeFieldType, 64>>;
  using QuadratureWeightsType = std::vector<BlockQuadratureWeightsType>;
  using AlphaReturnType = std::pair<BlockVectorType, std::pair<DomainType, RangeFieldType>>;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const bool disable_realizability_check,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , quad_points_(num_blocks)
    , quad_weights_(num_blocks)
    , M_(num_blocks, XT::LA::CommonDenseMatrix<RangeFieldType>())
    , tau_(tau)
    , disable_realizability_check_(disable_realizability_check)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
  {
    XT::LA::eye_matrix(T_minus_one_);
    if (!disable_realizability_check_)
      helper<dimFlux>::calculate_plane_coefficients(basis_functions_);
    const auto& quadratures = basis_functions_.quadratures();
    assert(quadratures.size() == num_blocks);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (const auto& quad_point : quadratures[jj]) {
        quad_points_[jj].emplace_back(quad_point.position());
        quad_weights_[jj].emplace_back(quad_point.weight());
      }
    } // jj
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      M_[jj] = XT::LA::CommonDenseMatrix<RangeFieldType>(quad_points_[jj].size(), block_size, 0., 0);
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = basis_functions_.evaluate_on_face(quad_points_[jj][ll], jj);
        for (size_t ii = 0; ii < block_size; ++ii)
          M_[jj].set_entry(ll, ii, val[ii]);
      } // ll
    } // jj
  }


  // ============================================================================================
  // ============================= FunctionInterface methods ====================================
  // ============================================================================================


  int order(const XT::Common::Parameter& /*param*/) const override
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
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, M_, eta_ast_prime_vals);
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = block_size * jj;
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto factor = eta_ast_prime_vals[jj][ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
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
    for (size_t dd = 0; dd < dimFlux; ++dd)
      row_jacobian(dd, M_, *H, result[dd], dd > 0);
  }

  void row_jacobian(const size_t row,
                    const BasisValuesMatrixType& M,
                    BlockMatrixType& H,
                    DynamicRowDerivativeRangeType& ret,
                    bool L_calculated = false) const
  {
    assert(row < dimFlux);
    calculate_J(M, ret, row);
    calculate_J_Hinv(ret, H, L_calculated);
  } // void partial_u_col(...)


  // ============================================================================================
  // ============ Evaluations of ansatz distribution, moments, hessian etc. =====================
  // ============================================================================================


  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // Solves the minimum entropy optimization problem for u.
  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();

    constexpr bool rescale = (entropy == EntropyType::MaxwellBoltzmann);

    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    static const auto alpha_one = std::make_unique<BlockVectorType>(basis_functions_.alpha_one());
    auto alpha_initial = std::make_unique<BlockVectorType>(*alpha_one);
    if (rescale) {
      *alpha_initial *= -std::log(density);
      *alpha_initial += alpha_in;
    } else {
      *alpha_initial = alpha_in;
    }

    if (!(density > 0. || !(basis_functions_.min_density(u) > 0.)) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    auto phi = std::make_unique<BlockVectorType>(u);
    if (rescale)
      *phi *= 1. / density;

    // if value has already been calculated for these values, skip computation
    RangeFieldType tau_prime = rescale ? std::min(tau_
                                                      / ((1 + std::sqrt(basis_dimRange) * phi->two_norm()) * density
                                                         + std::sqrt(basis_dimRange) * tau_),
                                                  tau_)
                                       : tau_;

    // calculate moment vector for isotropic distribution
    auto u_iso = std::make_unique<BlockVectorType>(basis_functions_.u_iso());
    if (!rescale)
      *u_iso *= density;

    // define further variables
    auto g_k = std::make_unique<BlockVectorType>();
    auto beta_out = std::make_unique<BlockVectorType>();
    auto v = std::make_unique<BlockVectorType>();
    thread_local auto T_k = XT::Common::make_unique<BlockMatrixType>();
    auto beta_in = std::make_unique<BlockVectorType>(*alpha_initial);

    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& rr : r_sequence) {
      // regularize u
      *v = *phi;
      if (rr > 0.) {
        *beta_in = *get_isotropic_alpha(*v);
        // calculate v = (1-r) u + r u_iso
        // use beta_out as storage for u_iso_in * r
        *v *= (1 - rr);
        *beta_out = *u_iso;
        *beta_out *= rr;
        *v += *beta_out;
      }
      for (size_t jj = 0; jj < num_blocks; ++jj)
        T_k->block(jj) = T_minus_one_;
      // calculate T_k u
      auto v_k = std::make_unique<BlockVectorType>(*v);
      // calculate values of basis p = S_k m
      thread_local BasisValuesMatrixType P_k(num_blocks, XT::LA::CommonDenseMatrix<RangeFieldType>(0, 0, 0., 0));
      copy_basis_matrix(M_, P_k);
      // calculate f_0
      RangeFieldType f_k = get_eta_ast_integrated(*beta_in, P_k) - *beta_in * *v_k;

      thread_local auto H = XT::Common::make_unique<BlockMatrixType>(0.);

      int backtracking_failed = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase r if too many iterations are used or cholesky decomposition fails
        if (kk > k_0_ && rr < r_max)
          break;
        try {
          change_basis(*beta_in, *v_k, P_k, *T_k, *g_k, *beta_out, *H);
        } catch (const Dune::MathError&) {
          if (rr < r_max)
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
        calculate_u(*alpha_tilde, M_, *u_alpha_tilde);
        *g_alpha_tilde = *u_alpha_tilde;
        *g_alpha_tilde -= *v;
        auto density_tilde = basis_functions_.density(*u_alpha_tilde);
        if (!(density_tilde > 0.) || !(basis_functions_.min_density(*u_alpha_tilde) > 0.) || std::isinf(density_tilde))
          break;

        // if rescale == true, we ensure that the ansatz density exactly matches the density of u
        if (rescale) {
          *alpha_prime = *alpha_one;
          *alpha_prime *= -std::log(density_tilde);
          *alpha_prime += *alpha_tilde;
          calculate_u(*alpha_prime, M_, *u_alpha_prime);
        } else {
          *alpha_prime = *alpha_tilde;
          *u_alpha_prime = *u_alpha_tilde;
        }

        // calculate stopping criteria
        *u_eps_diff = *u_alpha_prime;
        *u_eps_diff *= -(1 - epsilon_gamma_);
        *u_eps_diff += *v;

        // working_storage contains eta_ast_prime evaluations due to the calculate_u call above
        auto& eta_ast_prime_vals = working_storage();
        if (g_alpha_tilde->two_norm() < tau_prime
            && 1 - epsilon_gamma_ < std::exp(-(rescale ? d_alpha_tilde->one_norm() + std::abs(std::log(density_tilde))
                                                       : d_alpha_tilde->one_norm()))
            && (entropy == EntropyType::MaxwellBoltzmann || all_positive(eta_ast_prime_vals))
            && (disable_realizability_check_ || helper<dimFlux>::is_realizable(*u_eps_diff, basis_functions_))) {
          ret->second.first = *v;
          ret->second.second = rr;
          if (rescale) {
            ret->first = *alpha_one;
            ret->first *= std::log(density);
            ret->first += *alpha_prime;
            ret->second.first *= density;
          } else {
            ret->first = *alpha_prime;
          }
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          *beta_in = *beta_out;
          // backtracking line search
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * beta_out->two_norm() / d_k->two_norm()) {
            thread_local auto beta_new = std::make_unique<BlockVectorType>();
            *beta_new = *d_k;
            *beta_new *= zeta_k;
            *beta_new += *beta_out;
            RangeFieldType f = get_eta_ast_integrated(*beta_new, P_k) - *beta_new * *v_k;
            if (backtracking_failed >= 2 || f <= f_k + xi_ * zeta_k * (*g_k * *d_k)) {
              *beta_in = *beta_new;
              f_k = f;
              backtracking_failed = 0;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          if (zeta_k <= epsilon_ * beta_out->two_norm() / d_k->two_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    DUNE_THROW(MathError, "Failed to converge");

    return ret;
  }

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
    calculate_u(beta_out, P_k, g_k);
    g_k -= v_k;
  } // void change_basis(...)

  // returns density rho = < eta_ast_prime(beta_in * b(v)) >
  RangeFieldType get_rho(const BlockVectorType& beta_in, const BasisValuesMatrixType& M) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(beta_in, M, eta_ast_prime_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_prime_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  // returns < eta_ast(beta_in * b(v)) >
  RangeFieldType get_eta_ast_integrated(const BlockVectorType& beta_in, const BasisValuesMatrixType& M) const
  {
    auto& eta_ast_vals = working_storage();
    evaluate_eta_ast(beta_in, M, eta_ast_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  // returns < b \eta_{\ast}^{\prime}(\alpha^T b(v)) >
  DomainType get_u(const DomainType& alpha) const
  {
    BlockVectorType u_block;
    calculate_u(alpha, M_, u_block);
    return u_block;
  }

  DomainType get_u(const QuadratureWeightsType& eta_ast_prime_vals) const
  {
    DomainType ret;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto factor = eta_ast_prime_vals[jj][ll] * quad_weights_[jj][ll];
        const auto* basis_ll = &(M_[jj].get_entry_ref(ll, 0.));
        for (size_t ii = 0; ii < block_size; ++ii)
          ret[offset + ii] += basis_ll[ii] * factor;
      } // ll
    } // jj (intervals)
    return ret;
  }

  // calculate ret = \int (b eta_ast_prime(beta_in * m))
  void calculate_u(const BlockVectorType& beta_in, const BasisValuesMatrixType& M, BlockVectorType& ret) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(beta_in, M, eta_ast_prime_vals);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      auto& ret_jj = ret.block(jj);
      std::fill(ret_jj.begin(), ret_jj.end(), 0.);
      const size_t num_quad_points = quad_weights_[jj].size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        const auto factor = eta_ast_prime_vals[jj][ll] * quad_weights_[jj][ll];
        const auto* basis_ll = &(M[jj].get_entry_ref(ll, 0.));
        for (size_t ii = 0; ii < block_size; ++ii)
          ret_jj[ii] += basis_ll[ii] * factor;
      } // ll
    } // jj
  }

  void calculate_hessian(const BlockVectorType& alpha, const BasisValuesMatrixType& M, BlockMatrixType& H) const
  {
    auto& eta_ast_twoprime_vals = working_storage();
    evaluate_eta_ast_twoprime(alpha, M, eta_ast_twoprime_vals);
    calculate_hessian(eta_ast_twoprime_vals, M, H);
  } // void calculate_hessian(...)

  void calculate_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                         const BasisValuesMatrixType& M,
                         BlockMatrixType& H) const
  {
    // matrix is symmetric, we only use lower triangular part
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      std::fill(H.block(jj).begin(), H.block(jj).end(), 0.);
      const size_t num_quad_points = quad_weights_[jj].size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        auto factor_ll = eta_ast_twoprime_vals[jj][ll] * quad_weights_[jj][ll];
        const auto* basis_ll = &(M[jj].get_entry_ref(ll, 0.));
        for (size_t ii = 0; ii < block_size; ++ii) {
          auto* H_row = &(H.block(jj)[ii][0]);
          const auto factor_ll_ii = basis_ll[ii] * factor_ll;
          for (size_t kk = 0; kk <= ii; ++kk)
            H_row[kk] += basis_ll[kk] * factor_ll_ii;
        } // ii
      } // ll
    } // jj
  } // void calculate_hessian(...)

  void apply_inverse_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                             const DomainType& u,
                             DomainType& Hinv_u) const
  {
    thread_local auto H = XT::Common::make_unique<BlockMatrixType>();
    calculate_hessian(eta_ast_twoprime_vals, M_, *H);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      XT::LA::cholesky(H->block(jj));
    FieldVector<RangeFieldType, block_size> tmp_vec;
    FieldVector<RangeFieldType, block_size> block_u;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      // calculate C = A (L^T)^{-1}
      const auto offset = block_size * jj;
      for (size_t ii = 0; ii < block_size; ++ii)
        block_u[ii] = u[offset + ii];
      XT::LA::solve_lower_triangular(H->block(jj), tmp_vec, block_u);
      XT::LA::solve_lower_triangular_transposed(H->block(jj), block_u, tmp_vec);
      for (size_t ii = 0; ii < block_size; ++ii)
        Hinv_u[offset + ii] = block_u[ii];
    } // jj
  }

  // J = df/dalpha is the derivative of the flux with respect to alpha.
  // As F = (f_1, f_2, f_3) is matrix-valued
  // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
  // vector-valued),
  // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
  // this function returns the dd-th matrix df_dd/dalpha of J
  // assumes work_vecs already contains the needed eta_ast_twoprime(alpha * m) values
  void calculate_J(const BasisValuesMatrixType& M, DynamicRowDerivativeRangeType& J_dd, const size_t dd) const
  {
    assert(dd < dimFlux);
    const auto& eta_ast_twoprime_values = working_storage();
    // matrix is symmetric, we only use lower triangular part
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = jj * block_size;
      // set entries in this block to zero
      for (size_t ii = 0; ii < block_size; ++ii)
        std::fill_n(&(J_dd.get_entry_ref(offset + ii, offset)), ii + 1, 0.);
      // now calculate integral
      const size_t num_quad_points = quad_weights_[jj].size();
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        auto factor_ll = eta_ast_twoprime_values[jj][ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
        const auto* basis_ll = &(M[jj].get_entry_ref(ll, 0.));
        for (size_t ii = 0; ii < block_size; ++ii) {
          auto* J_row = &(J_dd[offset + ii][offset]);
          const auto factor_ll_ii = basis_ll[ii] * factor_ll;
          for (size_t kk = 0; kk <= ii; ++kk)
            J_row[kk] += basis_ll[kk] * factor_ll_ii;
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

  // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
  static void calculate_J_Hinv(DynamicRowDerivativeRangeType& A, BlockMatrixType& B, bool L_calculated = false)
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
  } // void calculate_J_Hinv(...)


  // ============================================================================================
  // ============================= Entropy evaluations ==========================================
  // ============================================================================================


  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast(ret);
  }

  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] = -std::log(1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_prime(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast_prime(ret);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast_prime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= (1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void
  evaluate_eta_ast_twoprime(const DomainType& alpha, const BasisValuesMatrixType& M, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, M, ret);
    apply_exponential(ret);
    evaluate_eta_ast_twoprime(ret);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already
  // contains exp(alpha^T b(v_i))
  void evaluate_eta_ast_twoprime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_blocks; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= std::pow(1 - ret[jj][ll], 2);
  }

  // stores evaluations of exp(alpha^T b(v_i)) for all quadrature points v_i
  void store_exp_evaluations(QuadratureWeightsType& exp_evaluations, const DomainType& alpha) const
  {
    this->calculate_scalar_products(alpha, M_, exp_evaluations);
    this->apply_exponential(exp_evaluations);
  }

  void store_eta_ast_prime_vals(const QuadratureWeightsType& exp_evaluations, QuadratureWeightsType& eta_ast_prime_vals)
  {
    eta_ast_prime_vals = exp_evaluations;
    evaluate_eta_ast_prime(eta_ast_prime_vals);
  }

  void store_eta_ast_twoprime_vals(const QuadratureWeightsType& exp_evaluations,
                                   QuadratureWeightsType& eta_ast_twoprime_vals)
  {
    eta_ast_twoprime_vals = exp_evaluations;
    evaluate_eta_ast_twoprime(eta_ast_twoprime_vals);
  }

  // stores evaluations of a given boundary distribution psi(v) at all quadrature points v_i
  void store_boundary_distribution_evaluations(
      QuadratureWeightsType& boundary_distribution_evaluations,
      const std::function<RangeFieldType(const FluxDomainType&)>& boundary_distribution) const
  {
    boundary_distribution_evaluations.resize(num_blocks);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      boundary_distribution_evaluations[jj].resize(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        boundary_distribution_evaluations[jj][ll] = boundary_distribution(quad_points_[jj][ll]);
    }
  }


  // ============================================================================================
  // =============================== Kinetic fluxes =============================================
  // ============================================================================================


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, \psi is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType
  evaluate_kinetic_flux(const DomainType& u_i, const DomainType& u_j, const FluxDomainType& n_ij, const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = std::make_unique<BlockVectorType>(get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first);
    const auto alpha_j = std::make_unique<BlockVectorType>(get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first);
    evaluate_kinetic_flux_with_alphas(*alpha_i, *alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const BlockVectorType& alpha_i,
                                               const BlockVectorType& alpha_j,
                                               const FluxDomainType& n_ij,
                                               const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    thread_local FieldVector<QuadratureWeightsType, 2> eta_ast_prime_vals{
        {QuadratureWeightsType(num_blocks), QuadratureWeightsType(num_blocks)}};
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      eta_ast_prime_vals[0][jj].resize(quad_points_[jj].size());
      eta_ast_prime_vals[1][jj].resize(quad_points_[jj].size());
    }
    evaluate_eta_ast_prime(alpha_i, M_, eta_ast_prime_vals[0]);
    evaluate_eta_ast_prime(alpha_j, M_, eta_ast_prime_vals[1]);
    DomainType ret(0);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto position = quad_points_[jj][ll][dd];
        RangeFieldType factor =
            position * n_ij[dd] > 0. ? eta_ast_prime_vals[0][jj][ll] : eta_ast_prime_vals[1][jj][ll];
        factor *= quad_weights_[jj][ll] * position;
        for (size_t ii = 0; ii < block_size; ++ii)
          ret[offset + ii] += M_[jj].get_entry(ll, ii) * factor;
      } // ll
    } // jj
    ret *= n_ij[dd];
    return ret;
  } // DomainType evaluate_kinetic_flux(...)

  // Calculates left and right kinetic flux with reconstructed densities. Ansatz distribution values contains
  // evaluations of the ansatz distribution at each quadrature point for a stencil of three entities. The distributions
  // are reconstructed pointwise for each quadrature point and the resulting (part of) the kinetic flux is <
  // psi_reconstr * b * v>_{+/-}.
  template <SlopeLimiterType slope_type, class FluxesMapType>
  void calculate_reconstructed_fluxes(const FieldVector<const QuadratureWeightsType*, 3>& ansatz_distribution_values,
                                      FluxesMapType& flux_values,
                                      const size_t dd) const
  {
    // get flux storage
    BasisDomainType coord(0.5);
    coord[dd] = 0;
    auto& left_flux_value = flux_values[coord];
    coord[dd] = 1;
    auto& right_flux_value = flux_values[coord];
    right_flux_value = left_flux_value = DomainType(0.);

    // evaluate exp(alpha^T b(v_i)) at all quadratures points v_i for all three alphas
    thread_local XT::Common::FieldVector<BlockQuadratureWeightsType, 2> reconstructed_values(
        BlockQuadratureWeightsType(quad_points_[0].size()));

    // get left and right reconstructed values for each quadrature point v_i
    auto& vals_left = reconstructed_values[0];
    auto& vals_right = reconstructed_values[1];

    const auto slope_func =
        (slope_type == SlopeLimiterType::minmod) ? XT::Common::minmod<RangeFieldType> : superbee<RangeFieldType>;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      // reconstruct densities
      if (slope_type == SlopeLimiterType::no_slope) {
        for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
          vals_left[ll] = vals_right[ll] = (*ansatz_distribution_values[1])[jj][ll];
      } else {
        for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
          const auto slope =
              slope_func((*ansatz_distribution_values[1])[jj][ll] - (*ansatz_distribution_values[0])[jj][ll],
                         (*ansatz_distribution_values[2])[jj][ll] - (*ansatz_distribution_values[1])[jj][ll]);
          vals_left[ll] = (*ansatz_distribution_values[1])[jj][ll] - 0.5 * slope;
          vals_right[ll] = (*ansatz_distribution_values[1])[jj][ll] + 0.5 * slope;
        } // ll
      }
      // calculate fluxes
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto position = quad_points_[jj][ll][dd];
        RangeFieldType factor = position > 0. ? vals_right[ll] : vals_left[ll];
        factor *= quad_weights_[jj][ll] * position;
        auto& val = position > 0. ? right_flux_value : left_flux_value;
        for (size_t ii = 0; ii < block_size; ++ii)
          val[offset + ii] += M_[jj].get_entry(ll, ii) * factor;
      } // ll
    } // jj
  } // void calculate_reconstructed_fluxes(...)


  // ============================================================================================
  // ================================== Helper functions ========================================
  // ============================================================================================


  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  // temporary vectors to store inner products and exponentials
  QuadratureWeightsType& working_storage() const
  {
    thread_local QuadratureWeightsType work_vecs(num_blocks);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      work_vecs[jj].resize(quad_points_[jj].size());
    return work_vecs;
  }

  bool all_positive(const QuadratureWeightsType& vals) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = vals[jj][ll];
        if (val < 0. || std::isinf(val) || std::isnan(val))
          return false;
      }
    return true;
  }

  std::unique_ptr<BlockVectorType> get_isotropic_alpha(const RangeFieldType density) const
  {
    return std::make_unique<BlockVectorType>(basis_functions_.alpha_iso(density));
  }

  std::unique_ptr<BlockVectorType> get_isotropic_alpha(const DomainType& u) const
  {
    return get_isotropic_alpha(basis_functions_.density(u));
  }

  void copy_basis_matrix(const BasisValuesMatrixType& source_mat, BasisValuesMatrixType& range_mat) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      range_mat[jj].backend() = source_mat[jj].backend();
  }

  void calculate_scalar_products_block(const size_t jj,
                                       const LocalVectorType& beta_in,
                                       const XT::LA::CommonDenseMatrix<RangeFieldType>& M,
                                       BlockQuadratureWeightsType& scalar_products) const
  {
    const size_t num_quad_points = quad_points_[jj].size();
    for (size_t ll = 0; ll < num_quad_points; ++ll) {
      const auto* basis_ll = &(M.get_entry_ref(ll, 0.));
      scalar_products[ll] = std::inner_product(beta_in.begin(), beta_in.end(), basis_ll, 0.);
    } // ll
  }

  void calculate_scalar_products(const BlockVectorType& beta_in,
                                 const BasisValuesMatrixType& M,
                                 QuadratureWeightsType& scalar_products) const
  {
    scalar_products.resize(num_blocks);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      scalar_products[jj].resize(quad_weights_[jj].size());
      calculate_scalar_products_block(jj, beta_in.block(jj), M[jj], scalar_products[jj]);
    }
  }

  void apply_exponential(BlockQuadratureWeightsType& values) const
  {
    assert(values.size() < std::numeric_limits<int>::max());
    XT::Common::Mkl::exp(static_cast<int>(values.size()), values.data(), values.data());
  }

  void apply_exponential(QuadratureWeightsType& values) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      apply_exponential(values[jj]);
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
        auto* basis_ll = &(M.get_entry_ref(ll, 0.));
        basis_ll[0] *= T_00_inv;
        basis_ll[1] = (basis_ll[1] - T_k[1][0] * basis_ll[0]) * T_11_inv;
      }
    } else if (block_size == 4) {
      FieldVector<RangeFieldType, 4> diag_inv;
      for (size_t ii = 0; ii < 4; ++ii)
        diag_inv[ii] = 1. / T_k[ii][ii];
      for (size_t ll = 0; ll < num_quad_points; ++ll) {
        auto* basis_ll = &(M.get_entry_ref(ll, 0.));
        basis_ll[0] *= diag_inv[0];
        basis_ll[1] = (basis_ll[1] - T_k[1][0] * basis_ll[0]) * diag_inv[1];
        basis_ll[2] = (basis_ll[2] - T_k[2][0] * basis_ll[0] - T_k[2][1] * basis_ll[1]) * diag_inv[2];
        basis_ll[3] =
            (basis_ll[3] - T_k[3][0] * basis_ll[0] - T_k[3][1] * basis_ll[1] - T_k[3][2] * basis_ll[2]) * diag_inv[3];
      }
    } else {
#  if HAVE_MKL
      thread_local LocalMatrixType T_k_trans(0.);
      assert(num_quad_points < std::numeric_limits<int>::max());
      // Calculate the transpose here first as this is much faster than passing the matrix to dtrsm and using
      // CblasTrans
      copy_transposed(T_k, T_k_trans);
      XT::Common::Cblas::dtrsm(XT::Common::Cblas::row_major(),
                               XT::Common::Cblas::right(),
                               XT::Common::Cblas::upper(),
                               XT::Common::Cblas::no_trans(),
                               XT::Common::Cblas::non_unit(),
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
        auto* M_row = &(M.get_entry_ref(ll, 0.));
        std::copy_n(M_row, block_size, tmp_vec.begin());
        XT::LA::solve_lower_triangular(T_k, tmp_vec2, tmp_vec);
        std::copy_n(tmp_vec2.begin(), block_size, M_row);
      }
#  endif
    }
  }

  void apply_inverse_matrix(const BlockMatrixType& T_k, BasisValuesMatrixType& M) const
  {
    for (size_t jj = 0; jj < num_blocks; ++jj)
      apply_inverse_matrix_block(jj, T_k.block(jj), M[jj]);
  }

  template <size_t domainDim = dimFlux, class anything = void>
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
      std::cerr << "Warning: You are missing Qhull, realizability stopping condition will not be checked!" << std::endl;
    }

    static bool is_realizable(const BlockVectorType& /*u*/, const MomentBasis& /*basis_functions*/)
    {
      return true;
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
        const auto& v0 = basis_functions.partitioning()[jj];
        const auto& v1 = basis_functions.partitioning()[jj + 1];
        bool ret = (u0 >= 0) && (u1 <= v1 * u0) && (v0 * u0 <= u1);
        if (!ret)
          return false;
      } // jj
      return true;
    }
  }; // class helper<1, ...>

  const MomentBasis& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const bool disable_realizability_check_;
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
template <class D, class R, size_t dimRange_or_refinements, size_t fluxDim, EntropyType entropy>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1, fluxDim, entropy>>
  : public XT::Functions::FunctionInterface<
        HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1, fluxDim, entropy>::dimRange,
        fluxDim,
        HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1, fluxDim, entropy>::dimRange,
        R>
{
public:
  using MomentBasis = HatFunctionMomentBasis<D, 3, R, dimRange_or_refinements, 1, fluxDim, entropy>;
  using BaseType =
      typename XT::Functions::FunctionInterface<MomentBasis::dimRange, MomentBasis::dimFlux, MomentBasis::dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using FluxDomainType = FieldVector<DomainFieldType, dimFlux>;
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
                                          const bool /*disable_realizability_check*/,
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
    , v_positive_(basis_functions_.triangulation().faces().size())
    , M_(basis_functions_.triangulation().faces().size())
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , num_faces_(basis_functions_.triangulation().faces().size())
  {
    const auto& triangulation = basis_functions_.triangulation();
    const auto& faces = triangulation.faces();
    assert(triangulation.vertices().size() == basis_dimRange);
    // create pattern
    XT::LA::SparsityPatternDefault pattern(basis_dimRange);
    for (size_t vertex_index = 0; vertex_index < basis_dimRange; ++vertex_index) {
      const auto& vertex = triangulation.vertices()[vertex_index];
      const auto& adjacent_faces = triangulation.get_face_indices(vertex->position());
      for (const auto& face_index : adjacent_faces) {
        const auto& face_vertices = faces[face_index]->vertices();
        assert(face_vertices.size() == 3);
        for (size_t jj = 0; jj < 3; ++jj)
          pattern.insert(vertex_index, face_vertices[jj]->index());
      }
    }
    pattern.sort();
    pattern_ = pattern;
    // store basis evaluations
    const auto& quadratures = basis_functions_.quadratures();
    assert(quadratures.size() == num_faces_);
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      const auto face_center = faces[jj]->center();
      for (size_t dd = 0; dd < 3; ++dd)
        v_positive_[jj][dd] = face_center[dd] > 0.;
      quad_points_[jj].resize(quadratures[jj].size());
      quad_weights_[jj].resize(quadratures[jj].size());
      for (size_t ll = 0; ll < quadratures[jj].size(); ++ll) {
        const auto& quad_point = quadratures[jj][ll];
        quad_points_[jj][ll] = quad_point.position();
        quad_weights_[jj][ll] = quad_point.weight();
      }
    } // jj
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      M_[jj].resize(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        M_[jj][ll] = basis_functions_.evaluate_on_face(quad_points_[jj][ll], jj);
    } // jj

  } // constructor


  // ============================================================================================
  // ============================= FunctionInterface methods ====================================
  // ============================================================================================


  int order(const XT::Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  RangeReturnType evaluate_with_alpha(const DomainType& alpha) const
  {
    const auto alpha_vec = XT::LA::convert_to<VectorType>(alpha);
    return evaluate_with_alpha(alpha_vec);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    LocalVectorType local_ret;
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    const auto& faces = basis_functions_.triangulation().faces();
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      // calculate ret[dd] = < omega[dd] m G_\alpha(u) >
      for (size_t jj = 0; jj < num_faces_; ++jj) {
        local_ret *= 0.;
        const auto& vertices = faces[jj]->vertices();
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto& basis_ll = M_[jj][ll];
          auto factor_ll = eta_ast_prime_vals[jj][ll] * quad_points_[jj][ll][dd] * quad_weights_[jj][ll];
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
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
    jacobian_with_alpha(alpha, result);
  }

  void jacobian_with_alpha(const DomainType& alpha, DynamicDerivativeRangeType& result) const
  {
    const auto alpha_vec = XT::LA::convert_to<VectorType>(alpha);
    jacobian_with_alpha(alpha_vec, result);
  }

  void jacobian_with_alpha(const VectorType& alpha, DynamicDerivativeRangeType& result) const
  {
    thread_local SparseMatrixType H(basis_dimRange, basis_dimRange, pattern_, 0);
    thread_local SparseMatrixType J(basis_dimRange, basis_dimRange, pattern_, 0);
    calculate_hessian(alpha, M_, H);
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      calculate_J(M_, J, dd);
      calculate_J_Hinv(J, H, result[dd]);
    }
  } // ... jacobian(...)


  // ============================================================================================
  // ============ Evaluations of ansatz distribution, moments, hessian etc. =====================
  // ============================================================================================


  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // Solves the minimum entropy optimization problem for u.
  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();

    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");

    constexpr bool rescale = (entropy == EntropyType::MaxwellBoltzmann);

    // rescale u such that the density <psi> is 1 if rescale is true
    VectorType phi(basis_dimRange, 0., 0);
    for (size_t ii = 0; ii < basis_dimRange; ++ii)
      phi.set_entry(ii, rescale ? u[ii] / density : u[ii]);
    VectorType alpha_one(basis_dimRange, 0., 0);
    basis_functions_.alpha_one(alpha_one);

    RangeFieldType tau_prime =
        rescale
            ? std::min(
                  tau_ / ((1 + std::sqrt(basis_dimRange) * phi.l2_norm()) * density + std::sqrt(basis_dimRange) * tau_),
                  tau_)
            : tau_;
    thread_local SparseMatrixType H(basis_dimRange, basis_dimRange, pattern_, 0);
#  if HAVE_EIGEN
    typedef ::Eigen::SparseMatrix<RangeFieldType, ::Eigen::ColMajor> ColMajorBackendType;
    typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
    thread_local SolverType solver;
    ColMajorBackendType colmajor_copy(H.backend());
#  else // HAVE_EIGEN
    thread_local auto solver = XT::LA::make_solver(H);
#  endif // HAVE_EIGEN

    // calculate moment vector for isotropic distribution
    VectorType u_iso(basis_dimRange, 0., 0);
    basis_functions_.u_iso(u_iso);
    if (!rescale)
      u_iso *= density;
    VectorType alpha_k = alpha_in;
    if (rescale)
      alpha_k -= alpha_one * std::log(density);
    VectorType v(basis_dimRange, 0., 0), g_k(basis_dimRange, 0., 0), d_k(basis_dimRange, 0., 0),
        tmp_vec(basis_dimRange, 0., 0), alpha_prime(basis_dimRange);
    const auto& r_sequence = regularize ? r_sequence_ : std::vector<RangeFieldType>{0.};
    const auto r_max = r_sequence.back();
    for (const auto& rr : r_sequence_) {
      // regularize u
      v = phi;
      if (rr > 0) {
        alpha_k = *get_isotropic_alpha(v);
        tmp_vec = u_iso;
        tmp_vec *= rr;
        v *= 1 - rr;
        v += tmp_vec;
      }

      // calculate f_0
      RangeFieldType f_k = calculate_f(alpha_k, v);

      int backtracking_failed = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase r if too many iterations are used
        if (kk > k_0_ && rr < r_max)
          break;
        // calculate gradient g
        calculate_gradient(alpha_k, v, g_k);
        // calculate Hessian H
        calculate_hessian(alpha_k, M_, H, entropy == EntropyType::MaxwellBoltzmann);
        // calculate descent direction d_k;
        tmp_vec = g_k;
        tmp_vec *= -1;
        try {
#  if HAVE_EIGEN
          colmajor_copy = H.backend();
          solver.analyzePattern(colmajor_copy);
          solver.factorize(colmajor_copy);
          d_k.backend() = solver.solve(tmp_vec.backend());
#  else // HAVE_EIGEN
          solver.apply(tmp_vec, d_k);
#  endif // HAVE_EIGEN
        } catch (const XT::LA::Exceptions::linear_solver_failed& error) {
          if (rr < r_max) {
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
        auto& u_eps_diff = tmp_vec;
        if (rescale) {
          alpha_prime = alpha_one;
          alpha_prime *= -std::log(density_tilde);
          alpha_prime += alpha_tilde;
          calculate_u(alpha_prime, u_eps_diff);
        } else {
          alpha_prime = alpha_tilde;
          u_eps_diff = u_alpha_tilde;
        }
        u_eps_diff *= -(1 - epsilon_gamma_);
        u_eps_diff += v;
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.l2_norm() < tau_prime && is_realizable(u_eps_diff)) {
          if (rescale) {
            ret->first = alpha_one;
            ret->first *= std::log(density);
            ret->first += alpha_prime;
          } else {
            ret->first = alpha_prime;
          }
          const auto v_ret_eig = rescale ? v * density : v;
          ret->second = std::make_pair(XT::LA::convert_to<DomainType>(v_ret_eig), rr);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          auto& alpha_new = tmp_vec;
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * alpha_k.l2_norm() / d_k.l2_norm()) {
            // calculate alpha_new = alpha_k + zeta_k d_k
            alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;
            // calculate f(alpha_new)
            RangeFieldType f_new = calculate_f(alpha_new, v);
            if (backtracking_failed >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              backtracking_failed = 0;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
          if (zeta_k <= epsilon_ * alpha_k.l2_norm() / d_k.l2_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    DUNE_THROW(MathError, "Failed to converge");
    return ret;
  } // ... get_alpha(...)

  // returns density rho = < eta_ast_prime(alpha * b(v)) >
  RangeFieldType get_rho(const DomainType& alpha) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_faces_; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_prime_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  // returns < eta_ast(alpha * b(v)) >
  RangeFieldType get_eta_ast_integrated(const DomainType& alpha) const
  {
    return get_eta_ast_integrated(XT::LA::convert_to<VectorType>(alpha));
  }

  // returns < eta_ast(alpha * b(v)) >
  RangeFieldType get_eta_ast_integrated(const VectorType& alpha) const
  {
    auto& eta_ast_vals = working_storage();
    evaluate_eta_ast(alpha, eta_ast_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_faces_; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  DomainType get_u(const DomainType& alpha) const
  {
    VectorType u(basis_dimRange, 0., 0);
    calculate_u(XT::LA::convert_to<VectorType>(alpha), u);
    return XT::LA::convert_to<DomainType>(u);
  }

  DomainType get_u(const QuadratureWeightsType& eta_ast_prime_vals) const
  {
    VectorType u(basis_dimRange, 0., 0);
    calculate_u(eta_ast_prime_vals, u);
    return XT::LA::convert_to<DomainType>(u);
  }

  // calculate ret = < b eta_ast_prime(alpha * b) >
  void calculate_u(const VectorType& alpha, VectorType& u) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    calculate_u(eta_ast_prime_vals, u);
  } // void calculate_u(...)

  // calculate ret = < b eta_ast_prime(alpha * b) >
  void calculate_u(const QuadratureWeightsType& eta_ast_prime_vals, VectorType& u) const
  {
    u *= 0.;
    LocalVectorType local_u;
    const auto& faces = basis_functions_.triangulation().faces();
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      const auto& vertices = faces[jj]->vertices();
      std::fill(local_u.begin(), local_u.end(), 0.);
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto factor = eta_ast_prime_vals[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          local_u[ii] += M_[jj][ll][ii] * factor;
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        u.add_to_entry(vertices[ii]->index(), local_u[ii]);
    } // jj (faces)
  } // void calculate_u(...)

  RangeFieldType calculate_f(const VectorType& alpha, const VectorType& v) const
  {
    return get_eta_ast_integrated(alpha) - alpha * v;
  } // void calculate_f(...)

  void calculate_gradient(const VectorType& alpha, const VectorType& v, VectorType& g_k) const
  {
    calculate_u(alpha, g_k);
    g_k -= v;
  }

  void calculate_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                         const BasisValuesMatrixType& M,
                         SparseMatrixType& H) const
  {
    H *= 0.;
    LocalMatrixType H_local(0.);
    const auto& faces = basis_functions_.triangulation().faces();
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      H_local *= 0.;
      const auto& vertices = faces[jj]->vertices();
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M[jj][ll];
        const auto factor = eta_ast_twoprime_vals[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          for (size_t kk = 0; kk < 3; ++kk)
            H_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * factor;
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        for (size_t kk = 0; kk < 3; ++kk)
          H.add_to_entry(vertices[ii]->index(), vertices[kk]->index(), H_local[ii][kk]);
    } // jj (faces)
  } // void calculate_hessian(...)

  void calculate_hessian(const VectorType& alpha,
                         const BasisValuesMatrixType& M,
                         SparseMatrixType& H,
                         const bool use_working_storage = false) const
  {
    auto& eta_ast_twoprime_vals = working_storage();
    if (!use_working_storage)
      evaluate_eta_ast_twoprime(alpha, eta_ast_twoprime_vals);
    calculate_hessian(eta_ast_twoprime_vals, M, H);
  } // void calculate_hessian(...)

  void apply_inverse_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                             const DomainType& u,
                             DomainType& Hinv_u) const
  {
    thread_local SparseMatrixType H(basis_dimRange, basis_dimRange, pattern_, 0);
    calculate_hessian(eta_ast_twoprime_vals, M_, H);
#  if HAVE_EIGEN
    thread_local VectorType u_vec(basis_dimRange, 0., 0), Hinv_u_vec(basis_dimRange, 0., 0);
    std::copy(u.begin(), u.end(), u_vec.begin());
    typedef ::Eigen::SparseMatrix<RangeFieldType, ::Eigen::ColMajor> ColMajorBackendType;
    ColMajorBackendType colmajor_copy(H.backend());
    typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
    SolverType solver;
    solver.analyzePattern(colmajor_copy);
    solver.factorize(colmajor_copy);
    Hinv_u_vec.backend() = solver.solve(u_vec.backend());
    std::copy(Hinv_u_vec.begin(), Hinv_u_vec.end(), Hinv_u.begin());
#  else // HAVE_EIGEN
    auto solver = XT::LA::make_solver(H);
    VectorType Hinv_u_la(basis_dimRange);
    VectorType u_la = XT::Common::convert_to<VectorType>(u);
    solver.apply(u_la, Hinv_u_la);
    Hinv_u = XT::Common::convert_to<DomainType>(Hinv_u_la);
#  endif
  } // void apply_inverse_hessian(..)

  // J = df/dalpha is the derivative of the flux with respect to alpha.
  // As F = (f_1, f_2, f_3) is matrix-valued
  // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
  // vector-valued),
  // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
  // this function returns the dd-th matrix df_dd/dalpha of J
  // assumes work_vecs already contains the needed \eta_{\ast}^{\prime \prime} (\alpha * b) values
  void calculate_J(const BasisValuesMatrixType& M, SparseMatrixType& J_dd, const size_t dd) const
  {
    assert(dd < dimFlux);
    J_dd *= 0.;
    LocalMatrixType J_local(0.);
    auto& eta_ast_twoprime_vals = working_storage();
    const auto& faces = basis_functions_.triangulation().faces();
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      J_local *= 0.;
      const auto& vertices = faces[jj]->vertices();
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M[jj][ll];
        const auto factor = eta_ast_twoprime_vals[jj][ll] * quad_points_[jj][ll][dd] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          for (size_t kk = 0; kk < 3; ++kk)
            J_local[ii][kk] += basis_ll[ii] * basis_ll[kk] * factor;
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        for (size_t kk = 0; kk < 3; ++kk)
          J_dd.add_to_entry(vertices[ii]->index(), vertices[kk]->index(), J_local[ii][kk]);
    } // jj (faces)
  } // void calculate_J(...)

  // calculates ret = J H^{-1}. H is assumed to be symmetric positive definite, which gives ret^T = H^{-T} J^T =
  // H^{-1} J^T, so we just have to solve y = H^{-1} x for each row x of J
  void calculate_J_Hinv(SparseMatrixType& J, const SparseMatrixType& H, DynamicRowDerivativeRangeType& ret) const
  {
    thread_local VectorType solution(basis_dimRange, 0., 0), tmp_rhs(basis_dimRange, 0., 0);
#  if HAVE_EIGEN
    typedef ::Eigen::SparseMatrix<RangeFieldType, ::Eigen::ColMajor> ColMajorBackendType;
    ColMajorBackendType colmajor_copy(H.backend());
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


  // ============================================================================================
  // ============================= Entropy evaluations ==========================================
  // ============================================================================================


  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast(ret);
  }

  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_faces_; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] = -std::log(1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_prime(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast_prime(ret);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast_prime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_faces_; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= (1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_twoprime(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast_twoprime(ret);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already
  // contains exp(alpha^T b(v_i))
  void evaluate_eta_ast_twoprime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_faces_; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= std::pow(1 - ret[jj][ll], 2);
  }

  // stores evaluations of exp(alpha^T b(v_i)) for all quadrature points v_i
  void store_exp_evaluations(QuadratureWeightsType& exp_evaluations, const DomainType& alpha) const
  {
    this->calculate_scalar_products(XT::LA::convert_to<VectorType>(alpha), exp_evaluations);
    this->apply_exponential(exp_evaluations);
  }

  void store_eta_ast_prime_vals(const QuadratureWeightsType& exp_evaluations, QuadratureWeightsType& eta_ast_prime_vals)
  {
    eta_ast_prime_vals = exp_evaluations;
    evaluate_eta_ast_prime(eta_ast_prime_vals);
  }

  void store_eta_ast_twoprime_vals(const QuadratureWeightsType& exp_evaluations,
                                   QuadratureWeightsType& eta_ast_twoprime_vals)
  {
    eta_ast_twoprime_vals = exp_evaluations;
    evaluate_eta_ast_twoprime(eta_ast_twoprime_vals);
  }

  // stores evaluations of a given boundary distribution psi(v) at all quadrature points v_i
  void store_boundary_distribution_evaluations(
      QuadratureWeightsType& boundary_distribution_evaluations,
      const std::function<RangeFieldType(const FluxDomainType&)>& boundary_distribution) const
  {
    boundary_distribution_evaluations.resize(num_faces_);
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      boundary_distribution_evaluations[jj].resize(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        boundary_distribution_evaluations[jj][ll] = boundary_distribution(quad_points_[jj][ll]);
    }
  }


  // ============================================================================================
  // =============================== Kinetic fluxes =============================================
  // ============================================================================================


  DomainType
  evaluate_kinetic_flux(const DomainType& u_i, const DomainType& u_j, const FluxDomainType& n_ij, const size_t dd) const
  {
    const auto alpha_i = get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first;
    return evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const DomainType& alpha_i,
                                               const DomainType& alpha_j,
                                               const FluxDomainType& n_ij,
                                               const size_t dd) const
  {
    return evaluate_kinetic_flux_with_alphas(
        XT::LA::convert_to<VectorType>(alpha_i), XT::LA::convert_to<VectorType>(alpha_j), n_ij, dd);
  }

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const FluxDomainType& n_ij,
                                               const size_t dd) const
  {
    thread_local FieldVector<QuadratureWeightsType, 2> eta_ast_prime_vals;
    eta_ast_prime_vals[0].resize(num_faces_);
    eta_ast_prime_vals[1].resize(num_faces_);
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      eta_ast_prime_vals[0][jj].resize(quad_points_[jj].size());
      eta_ast_prime_vals[1][jj].resize(quad_points_[jj].size());
    }
    evaluate_eta_ast_prime(alpha_i, eta_ast_prime_vals[0]);
    evaluate_eta_ast_prime(alpha_j, eta_ast_prime_vals[1]);
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    DomainType ret(0);
    const auto& faces = basis_functions_.triangulation().faces();
    LocalVectorType local_ret;
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      const bool positive_dir = v_positive_[jj][dd];
      const auto& eta_ast_prime = ((n_ij[dd] > 0. && positive_dir) || (n_ij[dd] < 0. && !positive_dir))
                                      ? eta_ast_prime_vals[0][jj]
                                      : eta_ast_prime_vals[1][jj];
      local_ret *= 0.;
      const auto& vertices = faces[jj]->vertices();
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        RangeFieldType factor = eta_ast_prime[ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
        for (size_t ii = 0; ii < 3; ++ii)
          local_ret[ii] += M_[jj][ll][ii] * factor;
      } // ll (quad points)
      for (size_t ii = 0; ii < 3; ++ii)
        ret[vertices[ii]->index()] += local_ret[ii];
    } // jj (faces)
    ret *= n_ij[dd];
    return ret;
  } // DomainType evaluate_kinetic_flux(...)

  // Calculates left and right kinetic flux with reconstructed densities. Ansatz distribution values contains
  // evaluations of the ansatz distribution at each quadrature point for a stencil of three entities. The distributions
  // are reconstructed pointwise for each quadrature point and the resulting (part of) the kinetic flux is <
  // psi_reconstr * b * v>_{+/-}.
  template <SlopeLimiterType slope_type, class FluxesMapType>
  void calculate_reconstructed_fluxes(const FieldVector<const QuadratureWeightsType*, 3>& ansatz_distribution_values,
                                      FluxesMapType& flux_values,
                                      const size_t dd) const
  {
    // get flux storage
    BasisDomainType coord(0.5);
    coord[dd] = 0;
    auto& left_flux_value = flux_values[coord];
    coord[dd] = 1;
    auto& right_flux_value = flux_values[coord];
    right_flux_value = left_flux_value = DomainType(0.);
    thread_local XT::Common::FieldVector<std::vector<RangeFieldType>, 2> reconstructed_values(
        std::vector<RangeFieldType>(quad_points_[0].size()));
    const auto& faces = basis_functions_.triangulation().faces();
    const auto slope_func =
        (slope_type == SlopeLimiterType::minmod) ? XT::Common::minmod<RangeFieldType> : superbee<RangeFieldType>;
    auto& vals_left = reconstructed_values[0];
    auto& vals_right = reconstructed_values[1];
    thread_local LocalVectorType face_flux(0.);
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      face_flux *= 0.;
      const bool positive_dir = v_positive_[jj][dd];
      auto& outside_vals = positive_dir ? vals_right : vals_left;
      const auto& vertices = faces[jj]->vertices();
      const size_t num_quad_points = quad_weights_[jj].size();
      const auto& non_reconstructed_values = (*ansatz_distribution_values[1])[jj];
      // reconstruct densities
      if (slope_type == SlopeLimiterType::no_slope) {
        for (size_t ll = 0; ll < num_quad_points; ++ll)
          outside_vals[ll] = non_reconstructed_values[ll];
      } else {
        const auto sign = positive_dir ? 1. : -1.;
        for (size_t ll = 0; ll < num_quad_points; ++ll) {
          const auto slope = slope_func(non_reconstructed_values[ll] - (*ansatz_distribution_values[0])[jj][ll],
                                        (*ansatz_distribution_values[2])[jj][ll] - non_reconstructed_values[ll]);
          outside_vals[ll] = non_reconstructed_values[ll] + sign * 0.5 * slope;
        } // ll (quad points)
      }
      // calculate fluxes
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        RangeFieldType factor = outside_vals[ll] * quad_weights_[jj][ll] * quad_points_[jj][ll][dd];
        const auto& basis_ll = M_[jj][ll];
        for (size_t ii = 0; ii < 3; ++ii)
          face_flux[ii] += basis_ll[ii] * factor;
      } // ll (quad points)
      auto& val = positive_dir ? right_flux_value : left_flux_value;
      for (size_t ii = 0; ii < 3; ++ii)
        val[vertices[ii]->index()] += face_flux[ii];
    } // jj
  } // void calculate_reconstructed_fluxes(...)


  // ============================================================================================
  // ================================== Helper functions ========================================
  // ============================================================================================


  std::unique_ptr<VectorType> get_isotropic_alpha(const RangeFieldType density) const
  {
    const auto alpha_iso_dynvector = basis_functions_.alpha_iso(density);
    auto ret = std::make_unique<VectorType>(alpha_iso_dynvector.size(), 0., 0);
    for (size_t ii = 0; ii < ret->size(); ++ii)
      (*ret)[ii] = alpha_iso_dynvector[ii];
    return ret;
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const DomainType& u) const
  {
    return get_isotropic_alpha(basis_functions_.density(u));
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const VectorType& u) const
  {
    auto u_domain = std::make_unique<DomainType>();
    for (size_t ii = 0; ii < dimFlux; ++ii)
      (*u_domain)[ii] = u.get_entry(ii);
    return get_isotropic_alpha(*u_domain);
  }

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  static bool is_realizable(const VectorType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  // temporary vectors to store inner products and exponentials
  std::vector<std::vector<RangeFieldType>>& working_storage() const
  {
    thread_local std::vector<std::vector<RangeFieldType>> work_vecs;
    work_vecs.resize(num_faces_);
    for (size_t jj = 0; jj < num_faces_; ++jj)
      work_vecs[jj].resize(quad_points_[jj].size());
    return work_vecs;
  }

  bool all_positive(const QuadratureWeightsType& vals) const
  {
    for (size_t jj = 0; jj < num_faces_; ++jj)
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = vals[jj][ll];
        if (val < 0. || std::isinf(val) || std::isnan(val))
          return false;
      }
    return true;
  }

  void calculate_scalar_products(const VectorType& alpha, QuadratureWeightsType& scalar_products) const
  {
    LocalVectorType local_alpha;
    const auto& faces = basis_functions_.triangulation().faces();
    scalar_products.resize(num_faces_);
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      scalar_products[jj].resize(quad_weights_[jj].size());
      const auto& vertices = faces[jj]->vertices();
      for (size_t ii = 0; ii < 3; ++ii)
        local_alpha[ii] = alpha.get_entry(vertices[ii]->index());
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
        scalar_products[jj][ll] = local_alpha * M_[jj][ll];
    } // jj
  }

  void apply_exponential(QuadratureWeightsType& values) const
  {
    for (size_t jj = 0; jj < num_faces_; ++jj) {
      assert(values[jj].size() < std::numeric_limits<int>::max());
      XT::Common::Mkl::exp(static_cast<int>(values[jj].size()), values[jj].data(), values[jj].data());
    }
  }

  const MomentBasis& basis_functions_;
  QuadraturePointsType quad_points_;
  QuadratureWeightsType quad_weights_;
  std::vector<std::bitset<3>> v_positive_;
  BasisValuesMatrixType M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const size_t num_faces_;
  XT::LA::SparsityPatternDefault pattern_;
};
#endif // ENTROPY_FLUX_USE_3D_HATFUNCTIONS_SPECIALIZATION

#if ENTROPY_FLUX_USE_1D_HATFUNCTIONS_SPECIALIZATION
#  if ENTROPY_FLUX_1D_HATFUNCTIONS_USE_ANALYTICAL_INTEGRALS
/**
 * Specialization of EntropyBasedFluxImplementation for 1D Hatfunctions with MaxwellBoltzmann entropy
 * (no change of basis, analytic integrals + Taylor)
 */
template <class D, class R, size_t dimRange, EntropyType entropy>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 1, R, dimRange, 1, 1, entropy>>
  : public XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>
{
  using BaseType = typename XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using MomentBasis = HatFunctionMomentBasis<D, 1, R, dimRange, 1, 1, entropy>;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using FluxDomainType = FieldVector<DomainFieldType, dimFlux>;
  using VectorType = XT::Common::FieldVector<RangeFieldType, basis_dimRange>;
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;
  static_assert(entropy == EntropyType::MaxwellBoltzmann, "Not implemented for other entropies");

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const bool /*disable_realizability_check*/,
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
    , v_points_(basis_functions_.partitioning())
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


  // ============================================================================================
  // ============================= FunctionInterface methods ====================================
  // ============================================================================================


  int order(const XT::Common::Parameter& /*param*/) const override
  {
    return 1;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
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
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
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


  // ============================================================================================
  // =============================== Kinetic fluxes =============================================
  // ============================================================================================


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, \psi is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType
  evaluate_kinetic_flux(const DomainType& u_i, const DomainType& u_j, const FluxDomainType& n_ij, const size_t dd) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first;
    evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const FluxDomainType& n_ij,
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


  // ============================================================================================
  // ============ Evaluations of ansatz distribution, moments, hessian etc. =====================
  // ============================================================================================


  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();
    // rescale u such that the density <psi> is 1
    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    static const auto alpha_one = basis_functions_.alpha_one();
    VectorType phi = u / density;
    VectorType alpha_initial = alpha_in - alpha_one * std::log(density);
    RangeFieldType tau_prime =
        std::min(tau_ / ((1 + std::sqrt(dimRange) * phi.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_);
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
    for (const auto& rr : r_sequence_) {
      // regularize u
      v = phi;
      if (rr > 0) {
        alpha_k = *get_isotropic_alpha(v);
        VectorType r_times_u_iso(u_iso);
        r_times_u_iso *= rr;
        v *= 1 - rr;
        v += r_times_u_iso;
      }

      // calculate f_0
      RangeFieldType f_k = calculate_f(alpha_k, v);

      int backtracking_failed = 0;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase r if too many iterations are used
        if (kk > k_0_ && rr < r_max)
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
          if (rr < r_max)
            break;
          else
            DUNE_THROW(Dune::MathError, "Failure to converge!");
        }

        const auto& alpha_tilde = alpha_k;
        const auto u_alpha_tilde = g_k + v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;
        const auto alpha_prime = alpha_tilde - alpha_one * std::log(density_tilde);
        const auto u_alpha_prime = calculate_u(alpha_prime);
        auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)) {
          ret->first = alpha_prime + alpha_one * std::log(density);
          ret->second = std::make_pair(v * density, rr);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
            // while (backtracking_failed >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.) {
            // calculate alpha_new = alpha_k + zeta_k d_k
            auto alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;
            // calculate f(alpha_new)
            RangeFieldType f_new = calculate_f(alpha_new, v);
            if (backtracking_failed >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              backtracking_failed = 0.;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
          if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    DUNE_THROW(MathError, "Failed to converge");

    return ret;
  } // ... get_alpha(...)

  DomainType get_u(const DomainType& alpha) const
  {
    return calculate_u(alpha);
  }

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


  // ============================================================================================
  // ================================== Helper functions ========================================
  // ============================================================================================


  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  static bool is_realizable(const DomainType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const RangeFieldType density) const
  {
    return std::make_unique<VectorType>(basis_functions_.alpha_iso(density));
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const DomainType& u) const
  {
    return get_isotropic_alpha(basis_functions_.density(u));
  }

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
template <class D, class R, size_t dimRange, EntropyType entropy>
class EntropyBasedFluxImplementation<HatFunctionMomentBasis<D, 1, R, dimRange, 1, 1, entropy>>
  : public XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>
{
  using BaseType = typename XT::Functions::FunctionInterface<dimRange, 1, dimRange, R>;
  using ThisType = EntropyBasedFluxImplementation;

public:
  using MomentBasis = HatFunctionMomentBasis<D, 1, R, dimRange, 1, 1, entropy>;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = dimRange;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRowDerivativeRangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using FluxDomainType = FieldVector<DomainFieldType, dimFlux>;
  using VectorType = XT::Common::FieldVector<RangeFieldType, basis_dimRange>;
  using AlphaReturnType = std::pair<VectorType, std::pair<DomainType, RangeFieldType>>;
  static const size_t num_intervals = dimRange - 1;
  static const size_t block_size = 2;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;
  using BasisValuesMatrixType = FieldVector<std::vector<LocalVectorType>, num_intervals>;
  using QuadraturePointsType = FieldVector<std::vector<RangeFieldType>, num_intervals>;
  using QuadratureWeightsType = QuadraturePointsType;

  explicit EntropyBasedFluxImplementation(const MomentBasis& basis_functions,
                                          const RangeFieldType tau,
                                          const bool /*disable_realizability_check*/,
                                          const RangeFieldType epsilon_gamma,
                                          const RangeFieldType chi,
                                          const RangeFieldType xi,
                                          const std::vector<RangeFieldType> r_sequence,
                                          const size_t k_0,
                                          const size_t k_max,
                                          const RangeFieldType epsilon)
    : basis_functions_(basis_functions)
    , grid_points_(basis_functions_.partitioning())
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


  // ============================================================================================
  // ============================= FunctionInterface methods ====================================
  // ============================================================================================


  int order(const XT::Common::Parameter& /*param*/) const override
  {
    return 1;
  }

  virtual RangeReturnType evaluate(const DomainType& u,
                                   const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
    return evaluate_with_alpha(alpha);
  }

  virtual RangeReturnType evaluate_with_alpha(const VectorType& alpha) const
  {
    RangeReturnType ret(0.);
    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M_[jj][ll];
        auto factor_ll = eta_ast_prime_vals[jj][ll] * quad_points_[jj][ll] * quad_weights_[jj][ll];
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
    const auto alpha = get_alpha(u, *get_isotropic_alpha(u), true)->first;
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


  // ============================================================================================
  // ============ Evaluations of ansatz distribution, moments, hessian etc. =====================
  // ============================================================================================


  std::unique_ptr<AlphaReturnType> get_alpha(const DomainType& u) const
  {
    return get_alpha(u, *get_isotropic_alpha(u), true);
  }

  // returns (alpha, (actual_u, r)), where r is the regularization parameter and actual_u the regularized u
  std::unique_ptr<AlphaReturnType>
  get_alpha(const DomainType& u, const VectorType& alpha_in, const bool regularize) const
  {
    auto ret = std::make_unique<AlphaReturnType>();

    constexpr bool rescale = (entropy == EntropyType::MaxwellBoltzmann);

    // rescale u such that the density <psi> is 1 if rescale is true
    RangeFieldType density = basis_functions_.density(u);
    if (!(density > 0.) || std::isinf(density))
      DUNE_THROW(Dune::MathError, "Negative, inf or NaN density!");
    static const auto alpha_one = basis_functions_.alpha_one();
    VectorType phi = rescale ? u / density : u;
    VectorType alpha_initial = alpha_in;
    if (rescale)
      alpha_initial -= alpha_one * std::log(density);
    RangeFieldType tau_prime =
        rescale
            ? std::min(tau_ / ((1 + std::sqrt(dimRange) * phi.two_norm()) * density + std::sqrt(dimRange) * tau_), tau_)
            : tau_;
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
    for (const auto& rr : r_sequence_) {
      // regularize u
      v = phi;
      if (rr > 0) {
        alpha_k = *get_isotropic_alpha(v);
        VectorType r_times_u_iso(u_iso);
        r_times_u_iso *= rr;
        v *= 1 - rr;
        v += r_times_u_iso;
      }

      // calculate f_0
      RangeFieldType f_k = calculate_f(alpha_k, v);

      int backtracking_failed = 0;
      VectorType g_k, d_k, minus_g_k, u_alpha_prime;
      for (size_t kk = 0; kk < k_max_; ++kk) {
        // exit inner for loop to increase r if too many iterations are used
        if (kk > k_0_ && rr < r_max)
          break;
        // calculate gradient g
        calculate_gradient(alpha_k, v, g_k);
        // calculate Hessian H
        calculate_hessian(alpha_k, M_, H_diag, H_subdiag, entropy == EntropyType::MaxwellBoltzmann);
        // calculate descent direction d_k;
        minus_g_k = g_k;
        minus_g_k *= -1;
        try {
          d_k = minus_g_k;
          XT::LA::solve_sym_tridiag_posdef(H_diag, H_subdiag, d_k);
        } catch (const Dune::MathError&) {
          if (rr < r_max)
            break;
          else
            DUNE_THROW(Dune::MathError, "Failure to converge!");
        }

        const auto& alpha_tilde = alpha_k;
        const auto u_alpha_tilde = g_k + v;
        auto density_tilde = basis_functions_.density(u_alpha_tilde);
        if (!(density_tilde > 0.) || std::isinf(density_tilde))
          break;
        auto alpha_prime = alpha_tilde;
        if (rescale) {
          alpha_prime -= alpha_one * std::log(density_tilde);
          calculate_u(alpha_prime, u_alpha_prime);
        } else {
          u_alpha_prime = u_alpha_tilde;
        }
        auto u_eps_diff = v - u_alpha_prime * (1 - epsilon_gamma_);
        auto& eta_ast_prime_vals = working_storage();
        // if rescale is true, working storage already contains the eta_ast_prime evaluations due to the call to
        // calculate_u above
        if (!rescale)
          evaluate_eta_ast_prime(alpha_prime, eta_ast_prime_vals);
        // checking realizability is cheap so we do not need the second stopping criterion
        if (g_k.two_norm() < tau_prime && is_realizable(u_eps_diff)
            && (entropy == EntropyType::MaxwellBoltzmann || all_positive(eta_ast_prime_vals))) {
          ret->first = rescale ? alpha_prime + alpha_one * std::log(density) : alpha_prime;
          ret->second = std::make_pair(rescale ? v * density : v, rr);
          return ret;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          while (backtracking_failed >= 2 || zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {
            // calculate alpha_new = alpha_k + zeta_k d_k
            auto alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;
            // calculate f(alpha_new)
            RangeFieldType f_new = calculate_f(alpha_new, v);
            if (backtracking_failed >= 2 || XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              backtracking_failed = 0.;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
          // if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm() * 100.)
          if (zeta_k <= epsilon_ * alpha_k.two_norm() / d_k.two_norm())
            ++backtracking_failed;
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // rr loop (Regularization parameter)
    DUNE_THROW(MathError, "Failed to converge");
    return ret;
  } // ... get_alpha(...)

  // returns density rho = < eta_ast_prime(alpha * b(v)) >
  RangeFieldType get_rho(const DomainType& alpha) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_intervals; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_prime_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  // returns < eta_ast(alpha * b(v)) >
  RangeFieldType get_eta_ast_integrated(const DomainType& alpha) const
  {
    auto& eta_ast_vals = working_storage();
    evaluate_eta_ast(alpha, eta_ast_vals);
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_intervals; ++jj)
      ret += std::inner_product(
          quad_weights_[jj].begin(), quad_weights_[jj].end(), eta_ast_vals[jj].begin(), RangeFieldType(0.));
    return ret;
  }

  DomainType get_u(const DomainType& alpha) const
  {
    DomainType ret;
    calculate_u(alpha, ret);
    return ret;
  }

  DomainType get_u(const QuadratureWeightsType& eta_ast_prime_vals) const
  {
    DomainType ret;
    calculate_u(eta_ast_prime_vals, ret);
    return ret;
  }

  void calculate_u(const DomainType& alpha, DomainType& u) const
  {
    auto& eta_ast_prime_vals = working_storage();
    evaluate_eta_ast_prime(alpha, eta_ast_prime_vals);
    calculate_u(eta_ast_prime_vals, u);
  } // void calculate_u(...)

  void calculate_u(const QuadratureWeightsType& eta_ast_prime_vals, DomainType& u) const
  {
    std::fill(u.begin(), u.end(), 0.);
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        auto factor_ll = eta_ast_prime_vals[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 2; ++ii)
          u[jj + ii] += M_[jj][ll][ii] * factor_ll;
      } // ll (quad points)
    } // jj (intervals)
  }

  RangeFieldType calculate_f(const VectorType& alpha, const VectorType& v) const
  {
    return get_eta_ast_integrated(alpha) - alpha * v;
  } // void calculate_u(...)

  void calculate_gradient(const VectorType& alpha, const VectorType& v, VectorType& g_k) const
  {
    calculate_u(alpha, g_k);
    g_k -= v;
  }

  void calculate_hessian(const QuadratureWeightsType& eta_ast_twoprime_vals,
                         const BasisValuesMatrixType& M,
                         VectorType& H_diag,
                         FieldVector<RangeFieldType, dimRange - 1>& H_subdiag) const
  {
    std::fill(H_diag.begin(), H_diag.end(), 0.);
    std::fill(H_subdiag.begin(), H_subdiag.end(), 0.);
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M[jj][ll];
        const auto factor = eta_ast_twoprime_vals[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 2; ++ii)
          H_diag[jj + ii] += std::pow(basis_ll[ii], 2) * factor;
        H_subdiag[jj] += basis_ll[0] * basis_ll[1] * factor;
      } // ll (quad points)
    } // jj (intervals)
  } // void calculate_hessian(...)

  void calculate_hessian(const DomainType& alpha,
                         const BasisValuesMatrixType& M,
                         VectorType& H_diag,
                         FieldVector<RangeFieldType, dimRange - 1>& H_subdiag,
                         const bool use_working_storage = false) const
  {
    auto& eta_ast_twoprime_vals = working_storage();
    if (!use_working_storage)
      evaluate_eta_ast_twoprime(alpha, eta_ast_twoprime_vals);
    calculate_hessian(eta_ast_twoprime_vals, M, H_diag, H_subdiag);
  } // void calculate_hessian(...)

  // J = df/dalpha is the derivative of the flux with respect to alpha.
  // As F = (f_1, f_2, f_3) is matrix-valued
  // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
  // vector-valued),
  // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
  // this function returns the dd-th matrix df_dd/dalpha of J
  // assumes work_vecs already contains the needed eta_ast_twoprime(alpha * b) values
  void calculate_J(const BasisValuesMatrixType& M,
                   VectorType& J_diag,
                   FieldVector<RangeFieldType, dimRange - 1>& J_subdiag) const
  {
    std::fill(J_diag.begin(), J_diag.end(), 0.);
    std::fill(J_subdiag.begin(), J_subdiag.end(), 0.);
    const auto& eta_ast_twoprime_vals = working_storage();
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& basis_ll = M[jj][ll];
        const auto factor = eta_ast_twoprime_vals[jj][ll] * quad_points_[jj][ll] * quad_weights_[jj][ll];
        for (size_t ii = 0; ii < 2; ++ii)
          J_diag[jj + ii] += std::pow(basis_ll[ii], 2) * factor;
        J_subdiag[jj] += basis_ll[0] * basis_ll[1] * factor;
      } // ll (quad points)
    } // jj (intervals)
  } // void calculate_J(...)

  void
  apply_inverse_hessian(const QuadratureWeightsType& density_evaluations, const DomainType& u, DomainType& Hinv_u) const
  {
    thread_local VectorType H_diag;
    thread_local FieldVector<RangeFieldType, dimRange - 1> H_subdiag;
    calculate_hessian(density_evaluations, M_, H_diag, H_subdiag);
    // factorize H = LDL^T, where L is unit lower bidiagonal and D is diagonal
    // H_diag is overwritten by the diagonal elements of D
    // H_subdiag is overwritten by the subdiagonal elements of L
    XT::LA::tridiagonal_ldlt(H_diag, H_subdiag);
    // Solve H ret = u
    Hinv_u = u;
    XT::LA::solve_tridiagonal_ldlt_factorized(H_diag, H_subdiag, Hinv_u);
  }

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


  // ============================================================================================
  // ============================= Entropy evaluations ==========================================
  // ============================================================================================


  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast(ret);
  }

  // evaluates \eta_{\ast}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_intervals; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] = -std::log(1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_prime(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast_prime(ret);
  }

  // evaluates \eta_{\ast}^{\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already contains
  // exp(alpha^T b(v_i))
  void evaluate_eta_ast_prime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_intervals; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= (1 - ret[jj][ll]);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i
  void evaluate_eta_ast_twoprime(const VectorType& alpha, QuadratureWeightsType& ret) const
  {
    calculate_scalar_products(alpha, ret);
    apply_exponential(ret);
    evaluate_eta_ast_twoprime(ret);
  }

  // evaluates \eta_{\ast}^{\prime\prime}(\alpha^T b(v_i)) for all quadrature points v_i, assumes that ret already
  // contains exp(alpha^T b(v_i))
  void evaluate_eta_ast_twoprime(QuadratureWeightsType& ret) const
  {
    if (entropy == EntropyType::BoseEinstein)
      for (size_t jj = 0; jj < num_intervals; ++jj)
        for (size_t ll = 0; ll < ret[jj].size(); ++ll)
          ret[jj][ll] /= std::pow(1 - ret[jj][ll], 2);
  }

  // stores evaluations of exp(alpha^T b(v_i)) for all quadrature points v_i
  void store_exp_evaluations(QuadratureWeightsType& exp_evaluations, const DomainType& alpha) const
  {
    this->calculate_scalar_products(XT::LA::convert_to<VectorType>(alpha), exp_evaluations);
    this->apply_exponential(exp_evaluations);
  }

  void store_eta_ast_prime_vals(const QuadratureWeightsType& exp_evaluations, QuadratureWeightsType& eta_ast_prime_vals)
  {
    eta_ast_prime_vals = exp_evaluations;
    evaluate_eta_ast_prime(eta_ast_prime_vals);
  }

  void store_eta_ast_twoprime_vals(const QuadratureWeightsType& exp_evaluations,
                                   QuadratureWeightsType& eta_ast_twoprime_vals)
  {
    eta_ast_twoprime_vals = exp_evaluations;
    evaluate_eta_ast_twoprime(eta_ast_twoprime_vals);
  }

  // stores evaluations of a given boundary distribution psi(v) at all quadrature points v_i
  void store_boundary_distribution_evaluations(
      QuadratureWeightsType& boundary_distribution_evaluations,
      const std::function<RangeFieldType(const FluxDomainType&)>& boundary_distribution) const
  {
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      boundary_distribution_evaluations[jj].resize(quad_points_[jj].size());
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll)
        boundary_distribution_evaluations[jj][ll] = boundary_distribution(quad_points_[jj][ll]);
    }
  }


  // ============================================================================================
  // =============================== Kinetic fluxes =============================================
  // ============================================================================================


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, \psi is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  DomainType
  evaluate_kinetic_flux(const DomainType& u_i, const DomainType& u_j, const FluxDomainType& n_ij, const size_t dd) const
  {
    const auto alpha_i = get_alpha(u_i, *get_isotropic_alpha(u_i), true)->first;
    const auto alpha_j = get_alpha(u_j, *get_isotropic_alpha(u_j), true)->first;
    return evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // DomainType evaluate_kinetic_flux(...)

  DomainType evaluate_kinetic_flux_with_alphas(const VectorType& alpha_i,
                                               const VectorType& alpha_j,
                                               const FluxDomainType& n_ij,
                                               const size_t dd) const
  {
    assert(dd == 0);
    thread_local FieldVector<QuadratureWeightsType, 2> eta_ast_prime_vals;
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      eta_ast_prime_vals[0][jj].resize(quad_points_[jj].size());
      eta_ast_prime_vals[1][jj].resize(quad_points_[jj].size());
    }
    evaluate_eta_ast_prime(alpha_i, eta_ast_prime_vals[0]);
    evaluate_eta_ast_prime(alpha_j, eta_ast_prime_vals[1]);
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    DomainType ret(0);
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto position = quad_points_[jj][ll];
        RangeFieldType factor =
            position * n_ij[dd] > 0. ? eta_ast_prime_vals[0][jj][ll] : eta_ast_prime_vals[1][jj][ll];
        factor *= quad_weights_[jj][ll] * position;
        for (size_t ii = 0; ii < 2; ++ii)
          ret[jj + ii] += M_[jj][ll][ii] * factor;
      } // ll (quad points)
    } // jj (faces)
    ret *= n_ij[dd];
    return ret;
  } // DomainType evaluate_kinetic_flux(...)

  // Calculates left and right kinetic flux with reconstructed densities. Ansatz distribution values contains
  // evaluations of the ansatz distribution at each quadrature point for a stencil of three entities. The distributions
  // are reconstructed pointwise for each quadrature point and the resulting (part of) the kinetic flux is <
  // psi_reconstr * b * v>_{+/-}.
  template <SlopeLimiterType slope_type, class FluxesMapType>
  void calculate_reconstructed_fluxes(const FieldVector<const QuadratureWeightsType*, 3>& ansatz_distribution_values,
                                      FluxesMapType& flux_values,
                                      const size_t dd) const
  {
    assert(dd == 0);
    // get flux storage
    BasisDomainType coord(0.5);
    coord[dd] = 0;
    auto& left_flux_value = flux_values[coord];
    coord[dd] = 1;
    auto& right_flux_value = flux_values[coord];
    right_flux_value = left_flux_value = DomainType(0.);
    thread_local XT::Common::FieldVector<std::vector<RangeFieldType>, 2> reconstructed_values(
        std::vector<RangeFieldType>(quad_points_[0].size()));
    const auto slope_func =
        (slope_type == SlopeLimiterType::minmod) ? XT::Common::minmod<RangeFieldType> : superbee<RangeFieldType>;
    auto& vals_left = reconstructed_values[0];
    auto& vals_right = reconstructed_values[1];
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      // reconstruct densities
      if (slope_type == SlopeLimiterType::no_slope) {
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
          vals_left[ll] = vals_right[ll] = (*ansatz_distribution_values[1])[jj][ll];
      } else {
        for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
          const auto slope =
              slope_func((*ansatz_distribution_values[1])[jj][ll] - (*ansatz_distribution_values[0])[jj][ll],
                         (*ansatz_distribution_values[2])[jj][ll] - (*ansatz_distribution_values[1])[jj][ll]);
          vals_left[ll] = (*ansatz_distribution_values[1])[jj][ll] - 0.5 * slope;
          vals_right[ll] = (*ansatz_distribution_values[1])[jj][ll] + 0.5 * slope;
        } // ll (quad points)
      }
      // calculate fluxes
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll) {
        const auto& position = quad_points_[jj][ll];
        RangeFieldType factor = position > 0. ? vals_right[ll] : vals_left[ll];
        factor *= quad_weights_[jj][ll] * position;
        auto& val = position > 0. ? right_flux_value : left_flux_value;
        const auto& basis_ll = M_[jj][ll];
        for (size_t ii = 0; ii < 2; ++ii)
          val[jj + ii] += basis_ll[ii] * factor;
      } // ll (quad points)
    } // jj
  } // void calculate_reconstructed_fluxes(...)


  // ============================================================================================
  // ================================== Helper functions ========================================
  // ============================================================================================


  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

  void calculate_scalar_products(const VectorType& alpha, QuadratureWeightsType& scalar_products) const
  {
    LocalVectorType local_alpha;
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      scalar_products[jj].resize(quad_weights_[jj].size());
      for (size_t ii = 0; ii < 2; ++ii)
        local_alpha[ii] = alpha[jj + ii];
      for (size_t ll = 0; ll < quad_weights_[jj].size(); ++ll)
        scalar_products[jj][ll] = local_alpha * M_[jj][ll];
    } // jj
  }

  void apply_exponential(QuadratureWeightsType& values) const
  {
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      assert(values[jj].size() < std::numeric_limits<int>::max());
      XT::Common::Mkl::exp(static_cast<int>(values[jj].size()), values[jj].data(), values[jj].data());
    }
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const RangeFieldType density) const
  {
    return std::make_unique<VectorType>(basis_functions_.alpha_iso(density));
  }

  std::unique_ptr<VectorType> get_isotropic_alpha(const DomainType& u) const
  {
    return get_isotropic_alpha(basis_functions_.density(u));
  }

  static bool is_realizable(const DomainType& u)
  {
    for (const auto& u_i : u)
      if (!(u_i > 0.) || std::isinf(u_i))
        return false;
    return true;
  }

  QuadratureWeightsType& working_storage() const
  {
    thread_local QuadratureWeightsType work_vec;
    for (size_t jj = 0; jj < num_intervals; ++jj)
      work_vec[jj].resize(quad_points_[jj].size());
    return work_vec;
  }

  bool all_positive(const QuadratureWeightsType& vals) const
  {
    for (size_t jj = 0; jj < num_intervals; ++jj)
      for (size_t ll = 0; ll < quad_points_[jj].size(); ++ll) {
        const auto val = vals[jj][ll];
        if (val < 0. || std::isinf(val) || std::isnan(val))
          return false;
      }
    return true;
  }

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
