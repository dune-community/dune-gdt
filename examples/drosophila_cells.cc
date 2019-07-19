// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/solver/istl/saddlepoint.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/quadratic.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>

#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>


using namespace Dune;
using namespace Dune::GDT;


// some global defines
using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using PGV = XT::Grid::PeriodicGridView<GV>;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using PI = XT::Grid::extract_intersection_t<PGV>;
using MatrixType = XT::LA::IstlRowMajorSparseMatrix<double>;
using VectorType = XT::LA::IstlDenseVector<double>;

template <class DiscreteFunc>
static std::string get_filename(const std::string& prefix, const DiscreteFunc& func, const size_t step)
{
  return prefix + "_" + func.name() + "_" + XT::Common::to_string(step) + ".txt";
}

template <class DiscreteFunc>
static std::string
get_rankfilename(const std::string& prefix, const DiscreteFunc& func, const size_t step, const int rank)
{
  return prefix + "_rank_" + XT::Common::to_string(rank) + "_" + func.name() + "_" + XT::Common::to_string(step)
         + ".txt";
}

template <class DiscreteFunc>
static void write_to_textfile(const DiscreteFunc& func, const std::string& prefix, const size_t step, const double t)
{
  const auto& grid_view = func.space().grid_view();
  const auto local_func = func.local_function();
  // write one file per MPI rank
  std::ofstream rankfile(get_rankfilename(prefix, func, step, grid_view.comm().rank()));
  for (const auto& entity : elements(grid_view, Dune::Partitions::interiorBorder)) {
    local_func->bind(entity);
    const auto entity_center = entity.geometry().center();
    auto position = entity_center;
    assert(position.size() == d);
    const auto val = local_func->evaluate(entity.geometry().local(position), {"t", t});
    for (size_t ii = 0; ii < d; ++ii)
      rankfile << XT::Common::to_string(position[ii], 15) << " ";
    rankfile << val << std::endl;
  }
  rankfile.close();
  // Wait till files on all ranks are written
  grid_view.comm().barrier();
  // Merge files
  if (grid_view.comm().rank() == 0) {
    const std::string merged_file_name = get_filename(prefix, func, step);
    std::remove(merged_file_name.c_str());
    std::ofstream merged_file(merged_file_name, std::ios_base::binary | std::ios_base::app);
    for (int ii = 0; ii < grid_view.comm().size(); ++ii) {
      const std::string rankfile_to_merge_name = get_rankfilename(prefix, func, step, ii);
      std::ifstream rankfile_to_merge(rankfile_to_merge_name, std::ios_base::binary);
      merged_file << rankfile_to_merge.rdbuf();
      rankfile_to_merge.close();
      std::remove(rankfile_to_merge_name.c_str());
    } // ii
    merged_file.close();
  } // if (rank == 0)
} // void write_to_textfile()

template <class DiscreteFunc, class VecDiscreteFunc>
static void write_files(const bool visualize,
                        const VecDiscreteFunc& u,
                        const DiscreteFunc& p,
                        const VecDiscreteFunc& P,
                        const VecDiscreteFunc& Pnat,
                        const DiscreteFunc& phi,
                        const DiscreteFunc& phi_zero,
                        const DiscreteFunc& phinat,
                        const DiscreteFunc& mu,
                        const std::string& prefix,
                        const size_t step,
                        const double t)
{
  std::string postfix = "_" + Dune::XT::Common::to_string(step);
  if (visualize) {
    u.visualize(prefix + "_" + u.name() + postfix);
    p.visualize(prefix + "_" + p.name() + postfix);
    P.visualize(prefix + "_" + P.name() + postfix);
    Pnat.visualize(prefix + "_" + Pnat.name() + postfix);
    phi.visualize(prefix + "_" + phi.name() + postfix);
    phi_zero.visualize(prefix + "_" + phi_zero.name() + postfix);
    phinat.visualize(prefix + "_" + phinat.name() + postfix);
    mu.visualize(prefix + "_" + mu.name() + postfix);
  }
  write_to_textfile(u, prefix, step, t);
  write_to_textfile(p, prefix, step, t);
  write_to_textfile(P, prefix, step, t);
  write_to_textfile(Pnat, prefix, step, t);
  write_to_textfile(phi, prefix, step, t);
  write_to_textfile(phi_zero, prefix, step, t);
  write_to_textfile(phinat, prefix, step, t);
  write_to_textfile(mu, prefix, step, t);
}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_navier_stokes(VecDiscreteFunc& u,
                         DiscreteFunc& p,
                         const VecDiscreteFunc& P,
                         const VecDiscreteFunc& Pnat,
                         const DiscreteFunc& phi,
                         const DiscreteFunc& phinat,
                         const DiscreteFunc& /*mu*/,
                         const double Re,
                         const double Fa,
                         const double xi,
                         const double vol_domain,
                         const XT::Grid::BoundaryInfo<PI>& boundary_info)
{
  if (Re > 1e-10)
    DUNE_THROW(Dune::NotImplemented, "No Navier-Stokes solver implemented yet!");
  const auto& u_space = u.space();
  const auto& p_space = p.space();

  const auto& grid_view = u_space.grid_view();
  // Setup spaces and matrices and vectors
  // Equations are
  // \int \nabla u \nabla v - \int p div v = \int ff v
  // \int (div u) q = \int gg q
  // System is [A B; B^T C] [u; p] = [f; g]
  // Dimensions are: A: n x n, B: n x m, C: m x m, u: n, f: n, p: m, g: m
  const size_t m = u_space.mapper().size();
  const size_t n = p_space.mapper().size();
  auto pattern_A = make_element_sparsity_pattern(u_space, u_space, grid_view);
  auto pattern_B = make_element_sparsity_pattern(u_space, p_space, grid_view);
  auto pattern_C = make_element_sparsity_pattern(p_space, p_space, grid_view);
  MatrixType A(m, m, pattern_A);
  MatrixType B(m, n, pattern_B);
  MatrixOperator<MatrixType, PGV, d> A_operator(grid_view, u_space, u_space, A);
  MatrixOperator<MatrixType, PGV, 1, 1, d> B_operator(grid_view, p_space, u_space, B);
  // calculate A_{ij} as \int \nabla v_i \nabla v_j
  A_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(1.)));
  // calculate B_{ij} as \int \nabla p_i div(v_j)
  B_operator.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
      LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
  // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
  VectorType f_vector(m), g_vector(n, 0.), p_basis_integrated_vector(n);
  auto f_functional = make_vector_functional(u_space, f_vector);
  auto P_local = P.local_function();
  auto Pnat_local = Pnat.local_function();
  auto phi_local = phi.local_function();
  auto phinat_local = phinat.local_function();
  const auto Fa_inv = 1. / Fa;
  using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;
  using RangeFieldType = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
  XT::Functions::GenericGridFunction<E, d> f(
      /*order = */ 3 * u_space.max_polorder(),
      /*post_bind_func*/
      [&P_local, &Pnat_local, &phi_local, &phinat_local](const E& element) {
        P_local->bind(element);
        Pnat_local->bind(element);
        phi_local->bind(element);
        phinat_local->bind(element);
      },
      /*evaluate_func*/
      [&P_local, &Pnat_local, &phi_local, &phinat_local, Fa_inv, xi](const DomainType& x_local,
                                                                     const XT::Common::Parameter& param) {
        // evaluate P, Pnat, phi, phinat, \nabla P, \nabla Pnat, \nabla phi, div P, div Pnat and phi_tilde = (phi + 1)/2
        // return type of the jacobians is a FieldMatrix<r, d>
        const auto P = P_local->evaluate(x_local, param);
        const auto Pnat = Pnat_local->evaluate(x_local, param);
        const auto phi = phi_local->evaluate(x_local, param)[0];
        const auto phinat = phinat_local->evaluate(x_local, param)[0];
        const auto grad_P = P_local->jacobian(x_local, param);
        const auto grad_Pnat = Pnat_local->jacobian(x_local, param);
        RangeFieldType div_P(0.), div_Pnat(0.);
        for (size_t ii = 0; ii < d; ++ii) {
          div_P += grad_P[ii][ii];
          div_Pnat += grad_Pnat[ii][ii];
        }
        const auto grad_phi = phi_local->jacobian(x_local, param)[0];
        const auto phi_tilde = (phi + 1) / 2;

        // evaluate rhs terms
        const auto Fa_inv_times_P_otimes_P_times_grad_phi_tilde = P * ((P * grad_phi) * Fa_inv * 0.5);
        auto grad_P_times_P = P;
        grad_P.mv(P, grad_P_times_P);
        const auto Fa_inv_times_phi_tilde_times_grad_P_times_P = grad_P_times_P * (phi_tilde * Fa_inv);
        const auto Fa_inv_times_phi_tilde_times_P_times_div_P = P * (div_P * phi_tilde * Fa_inv);
        const auto xi_plus = (xi + 1.) / 2.;
        const auto xi_minus = (xi - 1.) / 2.;
        auto grad_Pnat_times_P = P;
        grad_Pnat.mv(P, grad_Pnat_times_P);
        const auto xi_plus_times_grad_Pnat_times_P = grad_Pnat_times_P * xi_plus;
        const auto xi_plus_times_Pnat_times_div_P = Pnat * (div_P * xi_plus);
        auto grad_P_times_Pnat = P;
        grad_P.mv(Pnat, grad_P_times_Pnat);
        const auto xi_minus_times_grad_P_times_Pnat = grad_P_times_Pnat * xi_minus;
        const auto xi_minus_times_P_times_div_Pnat = P * (div_Pnat * xi_minus);
        const auto phinat_grad_phi = grad_phi * phinat;
        auto grad_P_T_times_Pnat = P;
        grad_P.mtv(Pnat, grad_P_T_times_Pnat);
        return Fa_inv_times_P_otimes_P_times_grad_phi_tilde + Fa_inv_times_phi_tilde_times_grad_P_times_P
               + Fa_inv_times_phi_tilde_times_P_times_div_P + xi_plus_times_grad_Pnat_times_P
               + xi_plus_times_Pnat_times_div_P + xi_minus_times_grad_P_times_Pnat + xi_minus_times_P_times_div_Pnat
               + phinat_grad_phi + grad_P_T_times_Pnat;
      });
  f_functional.append(LocalElementIntegralFunctional<E, d>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), f)));
  A_operator.append(f_functional);
  auto p_basis_integrated_functional = make_vector_functional(p_space, p_basis_integrated_vector);
  XT::Functions::ConstantGridFunction<E> one_function(1);
  p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), one_function)));
  B_operator.append(p_basis_integrated_functional);
  // Dirichlet constrainst for u
  auto dirichlet_constraints = make_dirichlet_constraints(u_space, boundary_info);
  A_operator.append(dirichlet_constraints);
  // assemble everything
  A_operator.assemble(true);
  B_operator.assemble(true);
  // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
  // Dirichlet DoF, so we also need to clear the row of B.
  dirichlet_constraints.apply(A, f_vector);
  for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs())
    B.clear_row(DoF);

  // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
  // [A B; B^T 0] to the unit vector.
  MatrixType C(n, n, pattern_C);
  size_t dof_index = 0;
  B.clear_col(dof_index);
  g_vector.set_entry(dof_index, 0.);
  C.set_entry(dof_index, dof_index, 1.);

  // now solve the system
  XT::LA::SaddlePointSolver<double> solver(A, B, B, C);
  // solve by schurcomplement (where the schur complement is inverted by CG and the inner
  // solves with A are using a direct method)
  std::string type = "cg_direct_schurcomplement";
  solver.apply(f_vector, g_vector, u.dofs().vector(), p.dofs().vector(), type);

  // ensure int_\Omega p = 0 (TODO: remove, not necessary as p is not used anywhere)
  auto p_integral = p_basis_integrated_vector * p.dofs().vector();
  auto p_correction = make_discrete_function<VectorType>(p_space, "p_corr");
  XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain);
  default_interpolation(const_p_integral_func, p_correction);
  p -= p_correction;
}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_ofield(const VecDiscreteFunc& u,
                  VecDiscreteFunc& P,
                  VecDiscreteFunc& Pnat,
                  const DiscreteFunc& phi,
                  const double xi,
                  const double kappa,
                  const double c_1,
                  const double Pa,
                  const double beta,
                  const double dt)
{
  const auto& P_space = P.space();
  const auto& Pnat_space = Pnat.space();

  const auto& grid_view = P_space.grid_view();
  // Setup spaces and matrices and vectors
  // System is [S_{00} S_{01}; S_{10} S_{11}] [P; Pnat] = [f_{of}; g_{of}]
  // All matrices have dimension n x n, all vectors have dimension n
  const size_t n = P_space.mapper().size();
  assert(Pnat_space.mapper().size() == n);
  // Use same pattern for all submatrices
  auto submatrix_pattern = make_element_sparsity_pattern(P_space, P_space, grid_view);
  XT::LA::SparsityPatternDefault pattern(2 * n);
  for (size_t ii = 0; ii < n; ++ii)
    for (const auto& jj : submatrix_pattern.inner(ii)) {
      pattern.insert(ii, jj);
      pattern.insert(ii, n + jj);
      pattern.insert(n + ii, jj);
      pattern.insert(n + ii, n + jj);
    }
  pattern.sort();
  // create system matrix and matrix views for the submatrices
  MatrixType S(2 * n, 2 * n, pattern);
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;
  MatrixViewType S_00(S, 0, n, 0, n);
  MatrixViewType S_01(S, 0, n, n, 2 * n);
  MatrixViewType S_10(S, n, 2 * n, 0, n);
  MatrixViewType S_11(S, n, 2 * n, n, 2 * n);
  // assemble matrix S_{00} = M + dt A
  MatrixOperator<MatrixViewType, PGV, d> S_00_operator(grid_view, P_space, Pnat_space, S_00);
  // calculate M_{ij} as \int \nabla \psi_i \nabla phi_j
  S_00_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(1.)));
  // calculate dt A
  auto u_local = u.local_function();
  using DomainType = typename XT::Functions::GenericGridFunction<E, d, d>::DomainType;
  using RangeFieldType = typename XT::Functions::GenericGridFunction<E, d, d>::RangeFieldType;
  // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
  XT::Functions::GenericGridFunction<E, d, d> dt_Omega_minus_xi_D(
      /*order = */ std::max(u.space().max_polorder() - 1, 0),
      /*post_bind_func*/
      [&u_local](const E& element) { u_local->bind(element); },
      /*evaluate_func*/
      [&u_local, xi, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate \nabla u
        auto grad_u = u_local->jacobian(x_local, param);
        auto grad_u_T = grad_u;
        grad_u_T.transpose();
        auto& ret = grad_u_T;
        ret *= dt * (1. - xi) / 2.;
        grad_u *= dt * (1 + xi) / 2.;
        ret -= grad_u;
        return ret;
      });
  S_00_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalElementQuadraticIntegrand<E>(dt_Omega_minus_xi_D)));
  XT::Functions::GenericGridFunction<E, d, 1> dt_u(
      /*order = */ u.space().max_polorder(),
      /*post_bind_func*/
      [&u_local](const E& element) { u_local->bind(element); },
      /*evaluate_func*/
      [&u_local, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
        auto ret = u_local->evaluate(x_local, param);
        ret *= dt;
        return ret;
      });
  S_00_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(dt_u)));

  // calculate S_01 = dt B
  MatrixOperator<MatrixViewType, PGV, d> S_01_operator(grid_view, Pnat_space, Pnat_space, S_01);
  S_01_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(dt / kappa)));

  // calculate S_10 = C
  MatrixOperator<MatrixViewType, PGV, d> S_10_operator(grid_view, P_space, P_space, S_10);
  S_10_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(-1. / Pa)));
  auto P_local = P.local_function();
  auto phi_local = phi.local_function();
  XT::Functions::GenericGridFunction<E, 1, 1> frac_c1_Pa_times_phi_minus_Pn_times_Pn(
      /*order = */ 2 * P.space().max_polorder(),
      /*post_bind_func*/
      [&P_local, &phi_local](const E& element) {
        P_local->bind(element);
        phi_local->bind(element);
      },
      /*evaluate_func*/
      [&P_local, &phi_local, c_1, Pa](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto P_n = P_local->evaluate(x_local, param);
        const auto phi = phi_local->evaluate(x_local, param);
        return c_1 / Pa * (phi - P_n * P_n);
      });
  S_10_operator.append(LocalElementIntegralBilinearForm<E, d>(
      LocalElementProductIntegrand<E, d>(frac_c1_Pa_times_phi_minus_Pn_times_Pn)));
  XT::Functions::GenericGridFunction<E, d, d> minus_two_frac_c1_Pa_Pn_otimes_Pn(
      /*order = */ 2 * P.space().max_polorder(),
      /*post_bind_func*/
      [&P_local](const E& element) { P_local->bind(element); },
      /*evaluate_func*/
      [&P_local, c_1, Pa](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto P_n = P_local->evaluate(x_local, param);
        FieldMatrix<RangeFieldType, d, d> ret;
        for (size_t ii = 0; ii < d; ++ii)
          for (size_t jj = 0; jj < d; ++jj)
            ret[ii][jj] = P_n[ii] * P_n[jj];
        ret *= -2. * c_1 / Pa;
        return ret;
      });
  S_10_operator.append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementQuadraticIntegrand<E>(minus_two_frac_c1_Pa_Pn_otimes_Pn)));

  // calculate S_11 = D
  MatrixOperator<MatrixViewType, PGV, d> S_11_operator(grid_view, Pnat_space, P_space, S_11);
  S_11_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(1.)));

  // create rhs = (f_{of}, g_{of}) = (0, g_{of})
  VectorType rhs_vector(2 * n, 0.);
  XT::LA::VectorView<VectorType> g_vector(rhs_vector, n, 2 * n);
  auto g_functional = make_vector_functional(P_space, g_vector);

  using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;
  using RangeFieldType = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
  XT::Functions::GenericGridFunction<E, d> g(
      /*order = */ 3 * P_space.max_polorder(),
      /*post_bind_func*/
      [&P_local, &phi_local](const E& element) {
        P_local->bind(element);
        phi_local->bind(element);
      },
      /*evaluate_func*/
      [&P_local, &phi_local, beta, Pa, c_1](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate P, phi, \nabla phi
        const auto P_n = P_local->evaluate(x_local, param);
        const auto grad_phi = phi_local->jacobian(x_local, param)[0];

        // evaluate rhs terms
        auto frac_beta_Pa_grad_phi = grad_phi;
        frac_beta_Pa_grad_phi *= beta / Pa;
        auto minus_two_frac_c1_Pa_Pn_times_Pn_Pn = P_n;
        minus_two_frac_c1_Pa_Pn_times_Pn_Pn *= -2. * c_1 / Pa * (P_n * P_n);
        return frac_beta_Pa_grad_phi + minus_two_frac_c1_Pa_Pn_times_Pn_Pn;
      });
  g_functional.append(LocalElementIntegralFunctional<E, d>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), g)));
  S_00_operator.append(g_functional);
  // assemble everything
  S_00_operator.assemble(true);
  S_01_operator.assemble(true);
  S_10_operator.assemble(true);
  S_11_operator.assemble(true);

  // solve system
  auto ret = XT::LA::solve(S, rhs_vector);

  // copy to vectors
  for (size_t ii = 0; ii < n; ++ii)
    P.dofs().vector()[ii] = ret[ii];
  for (size_t ii = 0; ii < n; ++ii)
    Pnat.dofs().vector()[ii] = ret[n + ii];
}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_pfield(const VecDiscreteFunc& u,
                  const VecDiscreteFunc& P,
                  DiscreteFunc& phi,
                  DiscreteFunc& phi_zero,
                  DiscreteFunc& phinat,
                  DiscreteFunc& mu,
                  const DiscreteFunc& g_D,
                  const double gamma,
                  const double c_1,
                  const double Pa,
                  const double Be,
                  const double Ca,
                  const double beta,
                  const double epsilon,
                  const double dt,
                  const XT::Grid::BoundaryInfo<PI>& boundary_info)
{
  const auto& grid_view = phi.space().grid_view();
  // Setup spaces and matrices and vectors
  // System is [S_{00} 0 S_{02}; S_{10} S_{11} 0; S_{20} S_{21} S_{22}] [phi; phinat; mu] = [f_{pf}; g_{pf}; h_{pf}]
  // All matrices have dimension n x n, all vectors have dimension n
  const size_t n = phi.space().mapper().size();
  assert(phinat.space().mapper().size() == n);
  assert(mu.space().mapper().size() == n);
  // Use same pattern for all submatrices
  auto submatrix_pattern = make_element_sparsity_pattern(phi.space(), phi.space(), grid_view);
  XT::LA::SparsityPatternDefault pattern(3 * n);
  for (size_t ii = 0; ii < n; ++ii)
    for (const auto& jj : submatrix_pattern.inner(ii)) {
      pattern.insert(ii, jj); // S_{00}
      pattern.insert(ii, 2 * n + jj); // S_{02}
      pattern.insert(n + ii, jj); // S_{10}
      pattern.insert(n + ii, n + jj); // S_{11}
      pattern.insert(2 * n + ii, jj); // S_{20}
      pattern.insert(2 * n + ii, n + jj); // S_{21}
      pattern.insert(2 * n + ii, 2 * n + jj); // S_{22}
    }
  pattern.sort();
  // create system matrix and matrix views for the submatrices
  MatrixType S(3 * n, 3 * n, pattern);
  MatrixType M(n, n, submatrix_pattern);
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;
  MatrixViewType S_00(S, 0, n, 0, n);
  MatrixViewType S_02(S, 0, n, 2 * n, 3 * n);
  MatrixViewType S_10(S, n, 2 * n, 0, n);
  MatrixViewType S_11(S, n, 2 * n, n, 2 * n);
  MatrixViewType S_20(S, 2 * n, 3 * n, 0, n);
  MatrixViewType S_21(S, 2 * n, 3 * n, n, 2 * n);
  MatrixViewType S_22(S, 2 * n, 3 * n, 2 * n, 3 * n);
  // assemble matrix S_{00} = A
  MatrixOperator<MatrixViewType, PGV, 1> S_00_operator(grid_view, phi.space(), phi.space(), S_00);
  const auto phi_zero_local = phi_zero.local_function();
  const auto g_D_local = g_D.local_function();
  using DomainType = typename XT::Functions::GenericGridFunction<E, 1, 1>::DomainType;
  XT::Functions::GenericGridFunction<E, 1, 1> A_prefactor(
      /*order = */ 2 * phi.space().max_polorder(),
      /*post_bind_func*/
      [&phi_zero_local, &g_D_local](const E& element) {
        phi_zero_local->bind(element);
        g_D_local->bind(element);
      },
      /*evaluate_func*/
      [&phi_zero_local, &g_D_local, epsilon](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
        const auto g_D = g_D_local->evaluate(x_local, param);
        return 1. / epsilon * (3. * phi_zero_n * phi_zero_n + 6. * g_D * phi_zero_n + 3. * g_D * g_D - 1);
      });
  S_00_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(A_prefactor)));
  S_00_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalEllipticIntegrand<E, 1>(epsilon)));
  // assemble matrix S_{02} = C
  MatrixOperator<MatrixViewType, PGV, 1> S_02_operator(grid_view, phinat.space(), phi.space(), S_02);
  S_02_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
  // assemble matrix S_{10} = M + dt D
  const auto u_local = u.local_function();
  MatrixOperator<MatrixViewType, PGV, 1> S_10_operator(grid_view, phinat.space(), phi.space(), S_10);
  XT::Functions::GenericGridFunction<E, d, 1> minus_u_dt(
      /*order = */ u.space().max_polorder(),
      /*post_bind_func*/
      [&u_local](const E& element) { u_local->bind(element); },
      /*evaluate_func*/
      [&u_local, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
        auto ret = u_local->evaluate(x_local, param);
        ret *= -dt;
        return ret;
      });
  S_10_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1>(minus_u_dt)));
  S_10_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
  // assemble matrix S_{11} = dt E
  MatrixOperator<MatrixViewType, PGV, 1> S_11_operator(grid_view, phinat.space(), phi.space(), S_11);
  S_11_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalEllipticIntegrand<E, 1>(dt * gamma)));
  // assemble matrix S_{20} = G
  const auto mu_local = mu.local_function();
  MatrixOperator<MatrixViewType, PGV, 1> S_20_operator(grid_view, phinat.space(), phi.space(), S_20);
  XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
      /*order = */ 2 * phi.space().max_polorder(),
      /*post_bind_func*/
      [&phi_zero_local, &g_D_local, &mu_local](const E& element) {
        phi_zero_local->bind(element);
        g_D_local->bind(element);
        mu_local->bind(element);
      },
      /*evaluate_func*/
      [&phi_zero_local, &g_D_local, &mu_local, epsilon, Be](const DomainType& x_local,
                                                            const XT::Common::Parameter& param) {
        const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
        const auto g_D = g_D_local->evaluate(x_local, param);
        const auto mu_n = mu_local->evaluate(x_local, param);
        return 6. / (Be * std::pow(epsilon, 2)) * (phi_zero_n * mu_n + g_D * mu_n);
      });
  S_20_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(G_prefactor)));
  // assemble matrix S_{21} = H
  MatrixOperator<MatrixViewType, PGV, 1> S_21_operator(grid_view, phinat.space(), phi.space(), S_21);
  S_21_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
  // assemble matrix S_{22} = J
  MatrixOperator<MatrixViewType, PGV, 1> S_22_operator(grid_view, phinat.space(), phi.space(), S_22);
  XT::Functions::GenericGridFunction<E, 1, 1> J_prefactor(
      /*order = */ 2 * phi.space().max_polorder(),
      /*post_bind_func*/
      [&phi_zero_local, &g_D_local](const E& element) {
        phi_zero_local->bind(element);
        g_D_local->bind(element);
      },
      /*evaluate_func*/
      [&phi_zero_local, &g_D_local, Ca, epsilon, Be](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
        const auto g_D = g_D_local->evaluate(x_local, param);
        return 1. / Ca
               + 1. / (Be * std::pow(epsilon, 2))
                     * (3. * g_D * g_D - 1. + 3. * std::pow(phi_zero_n, 2) + 6. * g_D * phi_zero_n);
      });
  S_22_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(J_prefactor)));
  S_22_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalEllipticIntegrand<E, 1>(1. / Be)));

  // create rhs = (f_{pf}, g_{pf}, h_{pf})
  VectorType rhs_vector(3 * n, 0.);
  XT::LA::VectorView<VectorType> f_vector(rhs_vector, 0, n);
  XT::LA::VectorView<VectorType> g_vector(rhs_vector, n, 2 * n);
  XT::LA::VectorView<VectorType> h_vector(rhs_vector, 2 * n, 3 * n);
  auto f_functional = make_vector_functional(phi_zero.space(), f_vector);
  auto g_functional = make_vector_functional(phi.space(), g_vector);
  auto h_functional = make_vector_functional(phi.space(), h_vector);

  using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;
  XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
      /*order = */ 3 * phi.space().max_polorder(),
      /*post_bind_func*/
      [&phi_zero_local, &g_D_local](const E& element) {
        phi_zero_local->bind(element);
        g_D_local->bind(element);
      },
      /*evaluate_func*/
      [&phi_zero_local, &g_D_local, epsilon](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate phi_zero, g_D
        const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
        const auto g_D = g_D_local->evaluate(x_local, param);
        return 1. / epsilon
               * (2. * std::pow(phi_zero_n, 3) + 3. * g_D * std::pow(phi_zero_n, 2) + g_D - std::pow(g_D, 3));
      });
  f_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), f_pf)));
  f_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalEllipticIntegrand<E, 1>(-epsilon), g_D)));
  S_00_operator.append(f_functional);

  using RangeFieldType = typename XT::Functions::GenericGridFunction<E, 1>::RangeFieldType;
  g_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementGradientValueIntegrand<E>(u), g_D)));
  S_11_operator.append(g_functional);
  MatrixOperator<MatrixType, PGV, 1> M_operator(grid_view, phi.space(), phi.space(), M);
  M_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));

  const auto P_local = P.local_function();
  XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
      /*order = */ 3 * phi.space().max_polorder(),
      /*post_bind_func*/
      [&phi_zero_local, &g_D_local, &mu_local, &P_local](const E& element) {
        phi_zero_local->bind(element);
        g_D_local->bind(element);
        mu_local->bind(element);
        P_local->bind(element);
      },
      /*evaluate_func*/
      [&phi_zero_local, &g_D_local, &mu_local, &P_local, Be, epsilon, c_1, Pa, beta](
          const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate phi_zero, g_D, mu, P, divP
        const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
        const auto g_D = g_D_local->evaluate(x_local, param);
        const auto mu_n = mu_local->evaluate(x_local, param);
        const auto Pn = P_local->evaluate(x_local, param);
        const auto grad_P = P_local->jacobian(x_local, param);
        RangeFieldType div_P(0.);
        for (size_t ii = 0; ii < d; ++ii)
          div_P += grad_P[ii][ii];
        return 6. / (Be * std::pow(epsilon, 2)) * (std::pow(phi_zero_n, 2) * mu_n + g_D * phi_zero_n * mu_n)
               - c_1 / (2. * Pa) * Pn * Pn - beta / Pa * div_P;
      });
  h_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), h_pf)));
  S_11_operator.append(h_functional);

  // append Dirichlet constraints for phi
  auto dirichlet_constraints = make_dirichlet_constraints(phi.space(), boundary_info);
  S_00_operator.append(dirichlet_constraints);

  // assemble everything
  S_00_operator.assemble(true);
  S_02_operator.assemble(true);
  S_10_operator.assemble(true);
  S_11_operator.assemble(true);
  S_20_operator.assemble(true);
  S_21_operator.assemble(true);
  S_22_operator.assemble(true);
  M_operator.assemble(true);

  // apply dirichlet constraints for phi
  dirichlet_constraints.apply(S_00, f_vector);
  for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs()) {
    S_02.clear_row(DoF);
    S_10.clear_col(DoF);
    S_20.clear_col(DoF);
  }

  // solve system
  g_vector *= dt;
  g_vector += M * phi_zero.dofs().vector();
  auto ret = XT::LA::solve(S, rhs_vector);

  // copy to vectors
  for (size_t ii = 0; ii < n; ++ii)
    phi_zero.dofs().vector()[ii] = ret[ii];
  for (size_t ii = 0; ii < n; ++ii)
    phinat.dofs().vector()[ii] = ret[n + ii];
  for (size_t ii = 0; ii < n; ++ii)
    mu.dofs().vector()[ii] = ret[2 * n + ii];
  phi.dofs().vector() = phi_zero.dofs().vector() + g_D.dofs().vector();
}

int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
#if HAVE_TBB
    DXTC_CONFIG.set("threading.partition_factor", 1, true);
    XT::Common::threadManager().set_max_threads(1);
#endif

    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // read configuration
    XT::Common::Configuration config("activepolargels.ini");

    // grid config
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(16));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(4));
    std::string periodic_dirs = config.get("grid.periodic_dirs", "01");

    // timestepping
    double t_end = config.template get<double>("fem.t_end", 340.);
    double dt = config.template get<double>("fem.dt", 0.005);

    // problem parameters
    double L = config.template get<double>("problem.L", 1e-6);
    double U = config.template get<double>("problem.U", 1e-6);
    double rho = config.template get<double>("problem.rho", 1.e3);
    double eta = config.template get<double>("problem.eta", 2.e3);
    double sigma = config.template get<double>("problem.sigma", 0.0188);
    double b_N = config.template get<double>("problem.b_N", 1.26e-14);
    double k = config.template get<double>("problem.k", 2.e-9);
    double xi = config.template get<double>("problem.xi", 1.1);
    double eta_rot = config.template get<double>("problem.eta_rot", 3.3e3);
    double zeta = config.template get<double>("problem.zeta", 2.e3);
    double epsilon = config.template get<double>("problem.epsilon", 0.21);
    double gamma = config.template get<double>("problem.gamma", 0.025);
    double c_1 = config.template get<double>("problem.c_1", 5.);
    double beta = config.template get<double>("problem.beta", 0.);
    double Re = rho * U * L / eta;
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;

    // output
    std::string filename = config.get("output.filename", "drosophila");
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);

    // create grid for [0, 160] x [0, 40] with periodic boundaries in x-direction
    FieldVector<double, d> lower_left{{0., 0.}};
    FieldVector<double, d> upper_right{{160., 40.}};
    auto grid = XT::Grid::make_cube_grid<G>(lower_left, upper_right, /*num_elements=*/{num_elements_x, num_elements_y});
    const double vol_domain = (upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]);
    auto nonperiodic_grid_view = grid.leaf_view();
    std::bitset<d> periodic_directions(periodic_dirs);
    PGV grid_view(nonperiodic_grid_view, periodic_directions);

    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    XT::Grid::AllDirichletBoundaryInfo<PI> dirichlet_boundary_info;
    XT::Grid::AllNeumannBoundaryInfo<PI> neumann_boundary_info;

    const XT::Functions::ConstantFunction<d> lambda(1);

    // create spaces for the variables
    // Taylor-Hood P2-P1 space for the Stokes variables (u, p) \in S_{stokes} = (H^1_0)^d x L^2\R, test functions (v, q)
    // \in S_{stokes}
    auto u_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    auto p_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/1);

    // P2 spaces for the orientation field variables P, P^{\natural} \in S_{ofield} = (H^1)^d x (L^2)^d, test functions
    // (Q, Q^\natural) \in S_{ofield}
    auto P_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    auto Pnat_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);

    // P2 spaces for the phase field variables \phi, \phi^\natural, \mu \in S_{pfield} = H^1_0 x H^1 x H^1, test
    // functions (\chi, \chi^\natural, \nu) \in S_{pfield}
    auto phi_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);
    auto phinat_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);
    auto mu_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);

    // create discrete functions for the variables
    auto u = make_discrete_function<VectorType>(u_space, "u");
    auto p = make_discrete_function<VectorType>(p_space, "p");
    auto P = make_discrete_function<VectorType>(P_space, "P");
    auto Pnat = make_discrete_function<VectorType>(Pnat_space, "Pnat");
    auto phi = make_discrete_function<VectorType>(phi_space, "phi");
    auto phi_zero = make_discrete_function<VectorType>(phi_space, "phizero");
    auto phinat = make_discrete_function<VectorType>(phinat_space, "phinat");
    auto mu = make_discrete_function<VectorType>(mu_space, "mu");
    auto g_D_phi = make_discrete_function<VectorType>(phi_space, "g_D_phi");

    // create and project initial values
    // we only need initial values for P and phi
    // Initially, cell is circular with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 160] \times [0, 40].
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
    FieldVector<double, d> center{upper_right[0] / 2, upper_right[1] / 2};
    auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
    const XT::Functions::GenericFunction<d> phi_initial(
        50,
        /*evaluate=*/
        [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2) * epsilon)); },
        /*name=*/"phi_initial");

    // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
    // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
    std::srand(1); // set seed for std::rand to 1
    const XT::Functions::GenericFunction<d, d> P_initial(3,
                                                         /*evaluate=*/
                                                         [phi_initial](const auto& x, const auto& param) {
                                                           auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           auto ret = FieldVector<double, d>({1. + rand1, 0. + rand2});
                                                           ret *= (phi_initial.evaluate(x, param) + 1.) / 2.;
                                                           return ret;
                                                         },
                                                         /*name=*/"P_initial");

    const XT::Functions::ConstantFunction<d> minus_one(-1.);

    // interpolate initial and boundary values
    default_interpolation(phi_initial, phi);
    default_interpolation(P_initial, P);
    XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
    boundary_interpolation(minus_one, g_D_phi, all_dirichlet_boundary_info, XT::Grid::DirichletBoundary{});
    phi_zero.dofs().vector() = phi.dofs().vector();
    phi_zero -= g_D_phi;

    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    write_files(true, u, p, P, Pnat, phi, phi_zero, phinat, mu, filename, 0, t);

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t << std::endl;
      solve_navier_stokes(u, p, P, Pnat, phi, phinat, mu, Re, Fa, xi, vol_domain, dirichlet_boundary_info);
      std::cout << "Stokes done" << std::endl;
      solve_ofield(u, P, Pnat, phi, xi, kappa, c_1, Pa, beta, dt);
      std::cout << "Ofield done" << std::endl;
      solve_pfield(
          u, P, phi, phi_zero, phinat, mu, g_D_phi, gamma, c_1, Pa, Be, Ca, beta, epsilon, dt, dirichlet_boundary_info);
      std::cout << "Pfield done" << std::endl;

      t += actual_dt;

      // check if data should be written in this timestep (and write)
      if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t, next_save_time)) {
        write_files(true, u, p, P, Pnat, phi, phi_zero, phinat, mu, filename, save_step_counter, t);
        next_save_time += write_step;
        ++save_step_counter;
      }
    } // while (t < t_end)
  } catch (Exception& e) {
    std::cerr << "\nDUNE reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "\nstl reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occured!" << std::endl;
    return EXIT_FAILURE;
  } // try
  return EXIT_SUCCESS;
} // ... main(...)
