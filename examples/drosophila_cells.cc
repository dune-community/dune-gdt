// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <chrono>
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
#include <dune/gdt/local/integrands/symmetric_elliptic.hh>
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
using G = ALU_2D_SIMPLEX_CONFORMING;
// using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using PGV = XT::Grid::PeriodicGridView<GV>;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using PI = XT::Grid::extract_intersection_t<PGV>;
using MatrixType = XT::LA::EigenRowMajorSparseMatrix<double>;
using MatrixViewType = XT::LA::MatrixView<MatrixType>;
using VectorType = XT::LA::EigenDenseVector<double>;
using R = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
using DiscreteFunctionType = DiscreteFunction<VectorType, PGV, 1, 1, R>;
using VectorDiscreteFunctionType = DiscreteFunction<VectorType, PGV, d, 1, R>;
using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;

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

static void write_files(const bool visualize,
                        const VectorDiscreteFunctionType& u,
                        const DiscreteFunctionType& p,
                        const VectorDiscreteFunctionType& P,
                        const VectorDiscreteFunctionType& Pnat,
                        const DiscreteFunctionType& phi,
                        const DiscreteFunctionType& phi_zero,
                        const DiscreteFunctionType& phinat,
                        const DiscreteFunctionType& mu,
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

enum class StokesSolverType
{
  schurcomplement_cg_direct,
  eigen_sparse_lu
};

struct StokesSolver
{
#if HAVE_EIGEN
  using ColMajorBackendType = ::Eigen::SparseMatrix<R, ::Eigen::ColMajor>;
  using SolverType = ::Eigen::SparseLU<ColMajorBackendType>;
#endif // HAVE_EIGEN

  StokesSolver(VectorDiscreteFunctionType& u,
               DiscreteFunctionType& p,
               const VectorDiscreteFunctionType& P,
               const VectorDiscreteFunctionType& Pnat,
               const DiscreteFunctionType& phi,
               const DiscreteFunctionType& phinat,
               const double Re,
               const double Fa,
               const double xi,
               const double vol_domain,
               const XT::Grid::BoundaryInfo<PI>& boundary_info,
               const StokesSolverType solver_type)
    : u_(u)
    , p_(p)
    , P_(P)
    , Pnat_(Pnat)
    , phi_(phi)
    , phinat_(phinat)
    , Re_(Re)
    , Fa_inv_(1. / Fa)
    , xi_(xi)
    , vol_domain_(vol_domain)
    , boundary_info_(boundary_info)
    , solver_type_(solver_type)
    , grid_view_(u_.space().grid_view())
    , m_(u_.space().mapper().size())
    , n_(p_.space().mapper().size())
    , pattern_A_(make_element_sparsity_pattern(u_.space(), u_.space(), grid_view_))
    , pattern_B_(make_element_sparsity_pattern(u_.space(), p_.space(), grid_view_))
    , pattern_C_(n_)
    , S_(m_ + n_, m_ + n_, create_pattern(m_, n_, pattern_A_, pattern_B_))
    , A_(m_, m_, pattern_A_)
    , B_(m_, n_, pattern_B_)
    , rhs_vector_(m_ + n_, 0.)
    , f_vector_(m_, 0.)
    , g_vector_(n_, 0.)
    , p_basis_integrated_vector_(n_)
    , dirichlet_constraints_(make_dirichlet_constraints(u_.space(), boundary_info_))
    , A_operator_(grid_view_, u_.space(), u_.space(), A_)
  {
    if (Re_ > 1e-2)
      DUNE_THROW(Dune::NotImplemented, "No Navier-Stokes solver implemented yet!");
    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
    // [A B; B^T 0] to the unit vector.
    const size_t dof_index = 0;
    pattern_C_.insert(dof_index, dof_index);
    C_ = MatrixType(n_, n_, pattern_C_);
    C_.set_entry(dof_index, dof_index, 1.);
    B_.clear_col(dof_index);
    g_vector_.set_entry(dof_index, 0.);

    MatrixOperator<MatrixType, PGV, 1, 1, d> B_operator(grid_view_, p_.space(), u_.space(), B_);
    // calculate A_{ij} as \int \nabla v_i \nabla v_j
    A_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalSymmetricEllipticIntegrand<E, d>(1.)));
    // calculate B_{ij} as \int \nabla p_i div(v_j)
    B_operator.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
        LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));

    auto p_basis_integrated_functional = make_vector_functional(p_.space(), p_basis_integrated_vector_);
    const XT::Functions::ConstantGridFunction<E> one_function(1);
    p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), one_function)));
    B_operator.append(p_basis_integrated_functional);

    // Dirichlet constrainst for u
    A_operator_.append(dirichlet_constraints_);
    // assemble everything
    A_operator_.assemble(false);
    B_operator.assemble(false);

    dirichlet_constraints_.apply(A_);
    for (const auto& DoF : dirichlet_constraints_.dirichlet_DoFs())
      B_.clear_row(DoF);

    if (solver_type_ == StokesSolverType::eigen_sparse_lu) {
#if HAVE_EIGEN
      for (size_t ii = 0; ii < m_; ++ii)
        for (const auto& jj : pattern_A_.inner(ii))
          S_.set_entry(ii, jj, A_.get_entry(ii, jj));

      for (size_t ii = 0; ii < m_; ++ii)
        for (const auto& jj : pattern_B_.inner(ii)) {
          S_.set_entry(ii, m_ + jj, B_.get_entry(ii, jj));
          S_.set_entry(m_ + jj, ii, -B_.get_entry(ii, jj));
        }

      for (size_t ii = 0; ii < n_; ++ii)
        for (const auto& jj : pattern_C_.inner(ii))
          S_.set_entry(m_ + ii, m_ + jj, C_.get_entry(ii, jj));

      S_colmajor_ = S_.backend();
      solver_.analyzePattern(S_colmajor_);
      solver_.factorize(S_colmajor_);
#else // HAVE_EIGEN
      DUNE_THROW(Dune::NotImplemented, "You are missing Eigen!");
#endif // HAVE_EIGEN
    }
  }

  static XT::LA::SparsityPatternDefault create_pattern(const size_t m,
                                                       const size_t n,
                                                       const XT::LA::SparsityPatternDefault& pattern_A,
                                                       const XT::LA::SparsityPatternDefault& pattern_B)
  {
    XT::LA::SparsityPatternDefault pattern(m + n);
    for (size_t ii = 0; ii < m; ++ii) {
      for (const auto& jj : pattern_A.inner(ii))
        pattern.insert(ii, jj);
      for (const auto& jj : pattern_B.inner(ii)) {
        pattern.insert(ii, m + jj);
        pattern.insert(m + jj, ii);
      }
    }
    pattern.insert(m, m);
    pattern.sort();
    return pattern;
  }

  void apply()
  {
    auto begin = std::chrono::steady_clock::now();

    auto f_functional = make_vector_functional(u_.space(), f_vector_);
    auto P_local = P_.local_function();
    auto Pnat_local = Pnat_.local_function();
    auto phi_local = phi_.local_function();
    auto phinat_local = phinat_.local_function();

    // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
    const R Fa_inv = Fa_inv_;
    const R xi = xi_;
    XT::Functions::GenericGridFunction<E, d> f(
        /*order = */ 3 * u_.space().max_polorder(),
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
          // evaluate P, Pnat, phi, phinat, \nabla P, \nabla Pnat, \nabla phi, div P, div Pnat and phi_tilde = (phi +
          // 1)/2 return type of the jacobians is a FieldMatrix<r, d>
          const auto P = P_local->evaluate(x_local, param);
          const auto Pnat = Pnat_local->evaluate(x_local, param);
          const auto phi = phi_local->evaluate(x_local, param)[0];
          const auto phinat = phinat_local->evaluate(x_local, param)[0];
          const auto grad_P = P_local->jacobian(x_local, param);
          const auto grad_Pnat = Pnat_local->jacobian(x_local, param);
          const auto grad_phi = phi_local->jacobian(x_local, param)[0];
          R div_P(0.), div_Pnat(0.);
          for (size_t ii = 0; ii < d; ++ii) {
            div_P += grad_P[ii][ii];
            div_Pnat += grad_Pnat[ii][ii];
          }
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
    A_operator_.clear();
    A_operator_.append(f_functional);
    f_vector_ *= 0.;
    A_operator_.assemble(false);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling Stokes took: " << time.count() << " s!" << std::endl;

    // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
    // Dirichlet DoF, so we also need to clear the row of B.
    dirichlet_constraints_.apply(f_vector_);

    // now solve the system
    begin = std::chrono::steady_clock::now();
    if (solver_type_ == StokesSolverType::schurcomplement_cg_direct) {
      XT::LA::SaddlePointSolver<VectorType, MatrixType> solver(A_, B_, B_, C_);
      // solve by schurcomplement (where the schur complement is inverted by CG and the inner
      // solves with A are using a direct method)
      static std::string type = "cg_direct_schurcomplement";
      solver.apply(f_vector_, g_vector_, u_.dofs().vector(), p_.dofs().vector(), type);
#if HAVE_EIGEN
    } else if (solver_type_ == StokesSolverType::eigen_sparse_lu) {
      for (size_t ii = 0; ii < m_; ++ii)
        rhs_vector_[ii] = f_vector_[ii];
      for (size_t ii = 0; ii < n_; ++ii)
        rhs_vector_[m_ + ii] = g_vector_[ii];
      VectorType ret(m_ + n_);
      ret.backend() = solver_.solve(rhs_vector_.backend());
      for (size_t ii = 0; ii < m_; ++ii)
        u_.dofs().vector()[ii] = ret[ii];
      for (size_t ii = 0; ii < n_; ++ii)
        p_.dofs().vector()[ii] = ret[m_ + ii];
#endif
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown solver type");
    }
    time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Stokes took: " << time.count() << " s!" << std::endl;

    // ensure int_\Omega p = 0 (TODO: remove, not necessary as p is not used anywhere)
    auto p_integral = p_basis_integrated_vector_ * p_.dofs().vector();
    auto p_correction = make_discrete_function<VectorType>(p_.space(), "p_corr");
    XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain_);
    default_interpolation(const_p_integral_func, p_correction);
    p_ -= p_correction;
  }

  VectorDiscreteFunctionType& u_;
  DiscreteFunctionType& p_;
  const VectorDiscreteFunctionType& P_;
  const VectorDiscreteFunctionType& Pnat_;
  const DiscreteFunctionType& phi_;
  const DiscreteFunctionType& phinat_;
  const double Re_;
  const double Fa_inv_;
  const double xi_;
  const double vol_domain_;
  const XT::Grid::BoundaryInfo<PI>& boundary_info_;
  const StokesSolverType solver_type_;
  const PGV& grid_view_;
  const size_t m_;
  const size_t n_;
  XT::LA::SparsityPatternDefault pattern_A_;
  XT::LA::SparsityPatternDefault pattern_B_;
  XT::LA::SparsityPatternDefault pattern_C_;
  MatrixType S_;
#if HAVE_EIGEN
  ColMajorBackendType S_colmajor_;
  SolverType solver_;
#endif // HAVE_EIGEN
  MatrixType A_;
  MatrixType B_;
  MatrixType C_;
  VectorType rhs_vector_;
  VectorType f_vector_;
  VectorType g_vector_;
  VectorType p_basis_integrated_vector_;
  DirichletConstraints<PI, SpaceInterface<PGV, d, 1, R>> dirichlet_constraints_;
  MatrixOperator<MatrixType, PGV, d> A_operator_;
};

struct OfieldSolver
{
  // Setup spaces and matrices and vectors
  // System is [S_{00} S_{01}; S_{10} S_{11}] [P; Pnat] = [f_{of}; g_{of}]
  // All matrices have dimension n x n, all vectors have dimension n
  // Use same pattern for all submatrices
  OfieldSolver(const VectorDiscreteFunctionType& u,
               VectorDiscreteFunctionType& P,
               VectorDiscreteFunctionType& Pnat,
               const DiscreteFunctionType& phi,
               const double xi,
               const double kappa,
               const double c_1,
               const double Pa,
               const double beta)
    : u_(u)
    , P_(P)
    , Pnat_(Pnat)
    , phi_(phi)
    , xi_(xi)
    , kappa_(kappa)
    , c_1_(c_1)
    , Pa_inv_(1. / Pa)
    , beta_(beta)
    , grid_view_(P_.space().grid_view())
    , n_(P_.space().mapper().size())
    , submatrix_pattern_(make_element_sparsity_pattern(P_.space(), P_.space(), grid_view_))
    , pattern_(create_pattern(n_, submatrix_pattern_))
    , S_(2 * n_, 2 * n_, pattern_)
    , M_(n_, n_, submatrix_pattern_)
    , S_00_(S_, 0, n_, 0, n_)
    , S_01_(S_, 0, n_, n_, 2 * n_)
    , S_10_(S_, n_, 2 * n_, 0, n_)
    , S_11_(S_, n_, 2 * n_, n_, 2 * n_)
    , S_00_operator_(grid_view_, P_.space(), Pnat_.space(), S_00_)
    , S_10_operator_(grid_view_, P_.space(), P_.space(), S_10_)
    , M_operator_(grid_view_, P_.space(), P_.space(), M_)
    , rhs_vector_(2 * n_, 0.)
    , f_vector_(rhs_vector_, 0, n_)
    , g_vector_(rhs_vector_, n_, 2 * n_)
  {
    // calculate M_{ij} as \int \psi_i phi_j
    M_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(1.)));
    M_operator_.assemble(false);
    // set S_11 = D = M
    S_11_ = M_;
  }

  static XT::LA::SparsityPatternDefault create_pattern(const size_t n,
                                                       const XT::LA::SparsityPatternDefault& submatrix_pattern)
  {
    XT::LA::SparsityPatternDefault pattern(2 * n);
    for (size_t ii = 0; ii < n; ++ii)
      for (const auto& jj : submatrix_pattern.inner(ii)) {
        pattern.insert(ii, jj);
        pattern.insert(ii, n + jj);
        pattern.insert(n + ii, jj);
        pattern.insert(n + ii, n + jj);
      }
    pattern.sort();
    return pattern;
  }

  void apply(const double dt)
  {
    auto begin = std::chrono::steady_clock::now();
    // calculate S_01 = dt B
    S_01_ = M_;
    S_01_ *= dt / kappa_;
    // assemble matrix S_{00} = M + dt A
    S_00_operator_.clear();
    S_00_ = M_;
    // calculate dt A
    auto u_local = u_.local_function();
    // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
    const R xi = xi_;
    XT::Functions::GenericGridFunction<E, d, d> dt_Omega_minus_xi_D(
        /*order = */ std::max(u_.space().max_polorder() - 1, 0),
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
    S_00_operator_.append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementQuadraticIntegrand<E>(dt_Omega_minus_xi_D)));
    XT::Functions::GenericGridFunction<E, d, 1> dt_u(
        /*order = */ u_.space().max_polorder(),
        /*post_bind_func*/
        [&u_local](const E& element) { u_local->bind(element); },
        /*evaluate_func*/
        [&u_local, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = u_local->evaluate(x_local, param);
          ret *= dt;
          return ret;
        });
    S_00_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(dt_u)));

    // calculate S_10 = C
    S_10_operator_.clear();
    S_10_ *= 0.;
    S_10_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(-Pa_inv_)));
    auto P_local = P_.local_function();
    auto phi_local = phi_.local_function();
    const R Pa_inv = Pa_inv_;
    const R c_1 = c_1_;
    XT::Functions::GenericGridFunction<E, 1, 1> frac_c1_Pa_times_phi_minus_Pn_times_Pn(
        /*order = */ 2 * P_.space().max_polorder(),
        /*post_bind_func*/
        [&P_local, &phi_local](const E& element) {
          P_local->bind(element);
          phi_local->bind(element);
        },
        /*evaluate_func*/
        [&P_local, &phi_local, c_1, Pa_inv](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = P_local->evaluate(x_local, param);
          const auto phi = phi_local->evaluate(x_local, param);
          return c_1 * Pa_inv * (phi - P_n * P_n);
        });
    S_10_operator_.append(LocalElementIntegralBilinearForm<E, d>(
        LocalElementProductIntegrand<E, d>(frac_c1_Pa_times_phi_minus_Pn_times_Pn)));
    XT::Functions::GenericGridFunction<E, d, d> minus_two_frac_c1_Pa_Pn_otimes_Pn(
        /*order = */ 2 * P_.space().max_polorder(),
        /*post_bind_func*/
        [&P_local](const E& element) { P_local->bind(element); },
        /*evaluate_func*/
        [&P_local, c_1, Pa_inv](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = P_local->evaluate(x_local, param);
          FieldMatrix<R, d, d> ret;
          for (size_t ii = 0; ii < d; ++ii)
            for (size_t jj = 0; jj < d; ++jj)
              ret[ii][jj] = P_n[ii] * P_n[jj];
          ret *= -2. * c_1 * Pa_inv;
          return ret;
        });
    S_10_operator_.append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementQuadraticIntegrand<E>(minus_two_frac_c1_Pa_Pn_otimes_Pn)));

    auto g_functional = make_vector_functional(P_.space(), g_vector_);
    g_vector_ *= 0.;
    const auto beta = beta_;
    XT::Functions::GenericGridFunction<E, d> g(
        /*order = */ 3 * P_.space().max_polorder(),
        /*post_bind_func*/
        [&P_local, &phi_local](const E& element) {
          P_local->bind(element);
          phi_local->bind(element);
        },
        /*evaluate_func*/
        [&P_local, &phi_local, beta, Pa_inv, c_1](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, phi, \nabla phi
          const auto P_n = P_local->evaluate(x_local, param);
          const auto grad_phi = phi_local->jacobian(x_local, param)[0];

          // evaluate rhs terms
          auto frac_beta_Pa_grad_phi = grad_phi;
          frac_beta_Pa_grad_phi *= beta * Pa_inv;
          auto minus_two_frac_c1_Pa_Pn_times_Pn_Pn = P_n;
          minus_two_frac_c1_Pa_Pn_times_Pn_Pn *= -2. * c_1 * Pa_inv * (P_n * P_n);
          return frac_beta_Pa_grad_phi + minus_two_frac_c1_Pa_Pn_times_Pn_Pn;
        });
    g_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), g)));
    S_10_operator_.append(g_functional);
    // assemble everything
    S_00_operator_.assemble(false);
    S_10_operator_.assemble(false);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    M_.mv(P_.dofs().vector(), f_vector_);
    std::cout << "Assembling Ofield took: " << time.count() << " s!" << std::endl;

    // solve system
    begin = std::chrono::steady_clock::now();
    const auto ret = XT::LA::solve(S_, rhs_vector_);
    time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Ofield took: " << time.count() << " s!" << std::endl;

    // copy to vectors
    for (size_t ii = 0; ii < n_; ++ii)
      P_.dofs().vector()[ii] = ret[ii];
    for (size_t ii = 0; ii < n_; ++ii)
      Pnat_.dofs().vector()[ii] = ret[n_ + ii];
  }

  const VectorDiscreteFunctionType& u_;
  VectorDiscreteFunctionType& P_;
  VectorDiscreteFunctionType& Pnat_;
  const DiscreteFunctionType& phi_;
  const R xi_;
  const R kappa_;
  const R c_1_;
  const R Pa_inv_;
  const R beta_;
  const PGV& grid_view_;
  const size_t n_;
  const XT::LA::SparsityPatternDefault submatrix_pattern_;
  const XT::LA::SparsityPatternDefault pattern_;
  MatrixType S_;
  MatrixType M_;
  MatrixViewType S_00_;
  MatrixViewType S_01_;
  MatrixViewType S_10_;
  MatrixViewType S_11_;
  MatrixOperator<MatrixViewType, PGV, d> S_00_operator_;
  MatrixOperator<MatrixViewType, PGV, d> S_10_operator_;
  MatrixOperator<MatrixType, PGV, d> M_operator_;
  VectorType rhs_vector_;
  XT::LA::VectorView<VectorType> f_vector_;
  XT::LA::VectorView<VectorType> g_vector_;
};

struct PfieldSolver
{
  // Setup spaces and matrices and vectors
  // System is [S_{00} 0 S_{02}; S_{10} S_{11} 0; S_{20} S_{21} S_{22}] [phi; phinat; mu] = [f_{pf}; g_{pf}; h_{pf}]
  // All matrices have dimension n x n, all vectors have dimension n
  PfieldSolver(const VectorDiscreteFunctionType& u,
               const VectorDiscreteFunctionType& P,
               DiscreteFunctionType& phi,
               DiscreteFunctionType& phi_zero,
               DiscreteFunctionType& phinat,
               DiscreteFunctionType& mu,
               const DiscreteFunctionType& g_D,
               const double gamma,
               const double c_1,
               const double Pa,
               const double Be,
               const double Ca,
               const double beta,
               const double epsilon,
               const XT::Grid::BoundaryInfo<PI>& boundary_info)
    : u_(u)
    , P_(P)
    , phi_(phi)
    , phi_zero_(phi_zero)
    , phinat_(phinat)
    , mu_(mu)
    , g_D_(g_D)
    , gamma_(gamma)
    , c_1_(c_1)
    , Pa_(Pa)
    , Be_(Be)
    , Ca_(Ca)
    , beta_(beta)
    , epsilon_(epsilon)
    , boundary_info_(boundary_info)
    , grid_view_(phi_.space().grid_view())
    , n_(phi_.space().mapper().size())
    , submatrix_pattern_(make_element_sparsity_pattern(phi_.space(), phi_.space(), grid_view_))
    , pattern_(create_pattern(n_, submatrix_pattern_))
    , S_(3 * n_, 3 * n_, pattern_)
    , M_(n_, n_, submatrix_pattern_)
    , E_(n_, n_, submatrix_pattern_)
    , S_00_(S_, 0, n_, 0, n_)
    , S_02_(S_, 0, n_, 2 * n_, 3 * n_)
    , S_10_(S_, n_, 2 * n_, 0, n_)
    , S_11_(S_, n_, 2 * n_, n_, 2 * n_)
    , S_20_(S_, 2 * n_, 3 * n_, 0, n_)
    , S_21_(S_, 2 * n_, 3 * n_, n_, 2 * n_)
    , S_22_(S_, 2 * n_, 3 * n_, 2 * n_, 3 * n_)
    , S_00_operator_(grid_view_, phi_.space(), phi_.space(), S_00_)
    , S_10_operator_(grid_view_, phinat_.space(), phi_.space(), S_10_)
    , S_20_operator_(grid_view_, mu_.space(), phi_.space(), S_20_)
    , S_22_operator_(grid_view_, mu_.space(), mu_.space(), S_22_)
    , rhs_vector_(3 * n_, 0.)
    , f_vector_(rhs_vector_, 0, n_)
    , g_vector_(rhs_vector_, n_, 2 * n_)
    , h_vector_(rhs_vector_, 2 * n_, 3 * n_)
    , dirichlet_constraints_(make_dirichlet_constraints(phi_.space(), boundary_info_))
  {
    assert(phinat_.space().mapper().size() == n_);
    assert(mu_.space().mapper().size() == n_);
    MatrixOperator<MatrixType, PGV, 1> M_operator(grid_view_, phi_.space(), phi_.space(), M_);
    M_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
    MatrixOperator<MatrixType, PGV, 1> E_operator(grid_view_, phinat_.space(), phinat_.space(), E_);
    E_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalEllipticIntegrand<E, 1>(1.)));
    M_operator.append(dirichlet_constraints_);
    M_operator.assemble(false);
    E_operator.assemble(false);
    // Set matrix S_{02} = C = M
    S_02_ = M_;
    // Set matrix S_{21} = H = M
    S_21_ = M_;
  }

  static XT::LA::SparsityPatternDefault create_pattern(const size_t n,
                                                       const XT::LA::SparsityPatternDefault& submatrix_pattern)
  {
    // Use same pattern for all submatrices
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
    return pattern;
  }

  void apply(const double dt)
  {
    auto begin = std::chrono::steady_clock::now();
    S_00_operator_.clear();
    S_10_operator_.clear();
    S_20_operator_.clear();
    S_22_operator_.clear();
    // assemble matrix S_{00} = A
    S_00_ = E_;
    S_00_ *= epsilon_;
    const auto phi_zero_local = phi_zero_.local_function();
    const auto g_D_local = g_D_.local_function();
    const R epsilon_inv = 1. / epsilon_;
    XT::Functions::GenericGridFunction<E, 1, 1> A_prefactor(
        /*order = */ 2 * phi_.space().max_polorder(),
        /*post_bind_func*/
        [&phi_zero_local, &g_D_local](const E& element) {
          phi_zero_local->bind(element);
          g_D_local->bind(element);
        },
        /*evaluate_func*/
        [&phi_zero_local, &g_D_local, epsilon_inv](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
          const auto g_D = g_D_local->evaluate(x_local, param);
          return epsilon_inv * (3. * std::pow(phi_zero_n, 2) + 6. * g_D * phi_zero_n + 3. * std::pow(g_D, 2) - 1);
        });
    S_00_operator_.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(A_prefactor)));
    // assemble matrix S_{10} = M + dt D
    S_10_ = M_;
    const auto u_local = u_.local_function();
    XT::Functions::GenericGridFunction<E, d, 1> minus_u_dt(
        /*order = */ u_.space().max_polorder(),
        /*post_bind_func*/
        [&u_local](const E& element) { u_local->bind(element); },
        /*evaluate_func*/
        [&u_local, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = u_local->evaluate(x_local, param);
          ret *= -dt;
          return ret;
        });
    S_10_operator_.append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u_dt)));
    // assemble matrix S_{11} = dt E
    S_11_ = E_;
    S_11_ *= gamma_ * dt;
    // assemble matrix S_{20} = G
    S_20_ *= 0.;
    const auto mu_local = mu_.local_function();
    const auto six_inv_Be_eps2 = 6. / (Be_ * std::pow(epsilon_, 2));
    XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
        /*order = */ 2 * phi_.space().max_polorder(),
        /*post_bind_func*/
        [&phi_zero_local, &g_D_local, &mu_local](const E& element) {
          phi_zero_local->bind(element);
          g_D_local->bind(element);
          mu_local->bind(element);
        },
        /*evaluate_func*/
        [&phi_zero_local, &g_D_local, &mu_local, six_inv_Be_eps2](const DomainType& x_local,
                                                                  const XT::Common::Parameter& param) {
          const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
          const auto g_D = g_D_local->evaluate(x_local, param);
          const auto mu_n = mu_local->evaluate(x_local, param);
          return six_inv_Be_eps2 * (phi_zero_n * mu_n + g_D * mu_n);
        });
    S_20_operator_.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(G_prefactor)));
    // assemble matrix S_{22} = J
    S_22_ = E_;
    S_22_ *= 1. / Be_;
    const R Ca_inv = 1. / Ca_;
    const R Be_eps2_inv = 1. / (Be_ * std::pow(epsilon_, 2));
    XT::Functions::GenericGridFunction<E, 1, 1> J_prefactor(
        /*order = */ 2 * phi_.space().max_polorder(),
        /*post_bind_func*/
        [&phi_zero_local, &g_D_local](const E& element) {
          phi_zero_local->bind(element);
          g_D_local->bind(element);
        },
        /*evaluate_func*/
        [&phi_zero_local, &g_D_local, Ca_inv, Be_eps2_inv](const DomainType& x_local,
                                                           const XT::Common::Parameter& param) {
          const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
          const auto g_D = g_D_local->evaluate(x_local, param);
          return Ca_inv
                 + Be_eps2_inv * (3. * std::pow(g_D, 2) - 1. + 3. * std::pow(phi_zero_n, 2) + 6. * g_D * phi_zero_n);
        });
    S_22_operator_.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(J_prefactor)));

    // create rhs = (f_{pf}, g_{pf}, h_{pf})
    auto f_functional = make_vector_functional(phi_zero_.space(), f_vector_);
    auto g_functional = make_vector_functional(phinat_.space(), g_vector_);
    auto h_functional = make_vector_functional(mu_.space(), h_vector_);

    XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
        /*order = */ 3 * phi_.space().max_polorder(),
        /*post_bind_func*/
        [&phi_zero_local, &g_D_local](const E& element) {
          phi_zero_local->bind(element);
          g_D_local->bind(element);
        },
        /*evaluate_func*/
        [&phi_zero_local, &g_D_local, epsilon_inv](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate phi_zero, g_D
          const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
          const auto g_D = g_D_local->evaluate(x_local, param);
          return epsilon_inv
                 * (2. * std::pow(phi_zero_n, 3) + 3. * g_D * std::pow(phi_zero_n, 2) + g_D - std::pow(g_D, 3));
        });
    f_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), f_pf)));
    f_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalEllipticIntegrand<E, 1>(-epsilon_), g_D_)));
    S_00_operator_.append(f_functional);

    XT::Functions::GenericGridFunction<E, d, 1> dt_u(
        /*order = */ u_.space().max_polorder(),
        /*post_bind_func*/
        [&u_local](const E& element) { u_local->bind(element); },
        /*evaluate_func*/
        [&u_local, dt](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = u_local->evaluate(x_local, param);
          ret *= dt;
          return ret;
        });
    g_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementGradientValueIntegrand<E>(dt_u), g_D_)));
    S_10_operator_.append(g_functional);

    const auto P_local = P_.local_function();
    const R c_1 = c_1_;
    const R Pa = Pa_;
    const R beta = beta_;
    XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
        /*order = */ 3 * phi_.space().max_polorder(),
        /*post_bind_func*/
        [&phi_zero_local, &g_D_local, &mu_local, &P_local](const E& element) {
          phi_zero_local->bind(element);
          g_D_local->bind(element);
          mu_local->bind(element);
          P_local->bind(element);
        },
        /*evaluate_func*/
        [&phi_zero_local, &g_D_local, &mu_local, &P_local, six_inv_Be_eps2, c_1, Pa, beta](
            const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate phi_zero, g_D, mu, P, divP
          const auto phi_zero_n = phi_zero_local->evaluate(x_local, param);
          const auto g_D = g_D_local->evaluate(x_local, param);
          const auto mu_n = mu_local->evaluate(x_local, param);
          const auto Pn = P_local->evaluate(x_local, param);
          const auto grad_P = P_local->jacobian(x_local, param);
          R div_P(0.);
          for (size_t ii = 0; ii < d; ++ii)
            div_P += grad_P[ii][ii];
          return six_inv_Be_eps2 * (std::pow(phi_zero_n, 2) * mu_n + g_D * phi_zero_n * mu_n)
                 - c_1 / (2. * Pa) * (Pn * Pn) - beta / Pa * div_P;
        });
    h_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), h_pf)));
    S_20_operator_.append(h_functional);

    // assemble everything
    f_vector_ *= 0.;
    h_vector_ *= 0.;
    M_.mv(phi_zero_.dofs().vector(), g_vector_);
    S_00_operator_.assemble(false);
    S_10_operator_.assemble(false);
    S_20_operator_.assemble(false);
    S_22_operator_.assemble(false);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling Pfield took: " << time.count() << " s!" << std::endl;

    // apply dirichlet constraints for phi
    dirichlet_constraints_.apply(S_00_, f_vector_);
    for (const auto& DoF : dirichlet_constraints_.dirichlet_DoFs()) {
      S_.unit_row(DoF);
      S_.unit_col(DoF);
    }

    // solve system
    begin = std::chrono::steady_clock::now();
    const auto ret = XT::LA::solve(S_, rhs_vector_);
    time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Pfield took: " << time.count() << " s!" << std::endl;

    // copy to vectors
    for (size_t ii = 0; ii < n_; ++ii)
      phi_zero_.dofs().vector()[ii] = ret[ii];
    for (size_t ii = 0; ii < n_; ++ii)
      phinat_.dofs().vector()[ii] = ret[n_ + ii];
    for (size_t ii = 0; ii < n_; ++ii)
      mu_.dofs().vector()[ii] = ret[2 * n_ + ii];
    phi_.dofs().vector() = phi_zero_.dofs().vector() + g_D_.dofs().vector();
  }

  const VectorDiscreteFunctionType& u_;
  const VectorDiscreteFunctionType& P_;
  DiscreteFunctionType& phi_;
  DiscreteFunctionType& phi_zero_;
  DiscreteFunctionType& phinat_;
  DiscreteFunctionType& mu_;
  const DiscreteFunctionType& g_D_;
  const double gamma_;
  const double c_1_;
  const double Pa_;
  const double Be_;
  const double Ca_;
  const double beta_;
  const double epsilon_;
  const XT::Grid::BoundaryInfo<PI>& boundary_info_;
  const PGV& grid_view_;
  const size_t n_;
  const XT::LA::SparsityPatternDefault submatrix_pattern_;
  const XT::LA::SparsityPatternDefault pattern_;
  MatrixType S_;
  MatrixType M_;
  MatrixType E_;
  MatrixViewType S_00_;
  MatrixViewType S_02_;
  MatrixViewType S_10_;
  MatrixViewType S_11_;
  MatrixViewType S_20_;
  MatrixViewType S_21_;
  MatrixViewType S_22_;
  MatrixOperator<MatrixViewType, PGV, 1> S_00_operator_;
  MatrixOperator<MatrixViewType, PGV, 1> S_10_operator_;
  MatrixOperator<MatrixViewType, PGV, 1> S_20_operator_;
  MatrixOperator<MatrixViewType, PGV, 1> S_22_operator_;
  VectorType rhs_vector_;
  XT::LA::VectorView<VectorType> f_vector_;
  XT::LA::VectorView<VectorType> g_vector_;
  XT::LA::VectorView<VectorType> h_vector_;
  DirichletConstraints<PI, SpaceInterface<PGV, 1, 1, R>> dirichlet_constraints_;
};

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
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

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
    FieldVector<double, d> center{upper_right[0] / 2., upper_right[1] / 2.};
    auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
    const XT::Functions::GenericFunction<d> phi_initial(
        50,
        /*evaluate=*/
        [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2) * epsilon)); },
        /*name=*/"phi_initial");

    // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
    // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
    std::srand(1); // set seed for std::rand to 1
    const XT::Functions::GenericFunction<d, d> P_initial(50,
                                                         /*evaluate=*/
                                                         [phi_initial](const auto& x, const auto& param) {
                                                           //                                                           auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           //                                                           auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           //                                                           auto ret = FieldVector<double, d>({1. + rand1, 0. + rand2});
                                                           auto ret = FieldVector<double, d>({1., 0.});
                                                           ret *= (phi_initial.evaluate(x, param) + 1.) / 2.;
                                                           return ret;
                                                         },
                                                         /*name=*/"P_initial");

    const XT::Functions::ConstantFunction<d> minus_one(-1.);
    const XT::Functions::ConstantFunction<d, d> u_initial(0.);

    // interpolate initial and boundary values
    //    default_interpolation(phi_initial, phi);

    const size_t m = phi.space().mapper().size();
    const auto M_pattern = make_element_sparsity_pattern(phi.space(), phi.space(), grid_view);
    MatrixType M(phi.space().mapper().size(), phi.space().mapper().size(), M_pattern);
    MatrixOperator<MatrixType, PGV, 1> M_operator(grid_view, phi.space(), phi.space(), M);
    VectorType rhs_vector(m, 0.);
    M_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
    auto f_functional = make_vector_functional(phi.space(), rhs_vector);
    f_functional.append(LocalElementIntegralFunctional<E, 1>(local_binary_to_unary_element_integrand(
        LocalElementProductIntegrand<E, 1>(), phi_initial.as_grid_function(grid_view))));
    M_operator.append(f_functional);
    M_operator.assemble(false);
    phi.dofs().vector() = XT::LA::solve(M, rhs_vector);


    //    default_interpolation(P_initial, P);

    const size_t m2 = P.space().mapper().size();
    const auto M2_pattern = make_element_sparsity_pattern(P.space(), P.space(), grid_view);
    MatrixType M2(m2, m2, M2_pattern);
    MatrixOperator<MatrixType, PGV, d> M2_operator(grid_view, P.space(), P.space(), M2);
    VectorType rhs_vector2(m2, 0.);
    M2_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(1.)));
    auto f2_functional = make_vector_functional(P.space(), rhs_vector2);
    f2_functional.append(LocalElementIntegralFunctional<E, d>(local_binary_to_unary_element_integrand(
        LocalElementProductIntegrand<E, d>(), P_initial.as_grid_function(grid_view))));
    M2_operator.append(f2_functional);
    M2_operator.assemble(false);
    P.dofs().vector() = XT::LA::solve(M2, rhs_vector2);


    default_interpolation(u_initial, u);
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

    OfieldSolver ofield_solver(u, P, Pnat, phi, xi, kappa, c_1, Pa, beta);
    PfieldSolver pfield_solver(
        u, P, phi, phi_zero, phinat, mu, g_D_phi, gamma, c_1, Pa, Be, Ca, beta, epsilon, dirichlet_boundary_info);
    StokesSolver stokes_solver(
        u, p, P, Pnat, phi, phinat, Re, Fa, xi, vol_domain, dirichlet_boundary_info, StokesSolverType::eigen_sparse_lu);

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t << std::endl;
      stokes_solver.apply();
      std::cout << "Stokes done" << std::endl;
      pfield_solver.apply(actual_dt);
      std::cout << "Pfield done" << std::endl;
      ofield_solver.apply(actual_dt);
      std::cout << "Ofield done" << std::endl;

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
