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
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/symmetric_elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/gradient_value.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/norms.hh>

#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>

#include <Eigen/IterativeLinearSolvers>


using namespace Dune;
using namespace Dune::GDT;


// some global defines
// using G = ALU_2D_SIMPLEX_CONFORMING;
using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using PGV = XT::Grid::PeriodicGridView<GV>;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using PI = XT::Grid::extract_intersection_t<PGV>;
using MatrixType = XT::LA::EigenRowMajorSparseMatrix<double>;
using VectorType = XT::LA::EigenDenseVector<double>;
// using MatrixType = XT::LA::IstlRowMajorSparseMatrix<double>;
// using VectorType = XT::LA::IstlDenseVector<double>;
using MatrixViewType = XT::LA::MatrixView<MatrixType>;
using VectorViewType = XT::LA::VectorView<VectorType>;
using R = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
using DiscreteFunctionType = DiscreteFunction<VectorType, PGV, 1, 1, R>;
using LocalDiscreteFunctionType = typename DiscreteFunctionType::LocalFunctionType;
using VectorDiscreteFunctionType = DiscreteFunction<VectorType, PGV, d, 1, R>;
using VectorLocalDiscreteFunctionType = typename VectorDiscreteFunctionType::LocalFunctionType;
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
                        const std::vector<VectorDiscreteFunctionType>& P,
                        const std::vector<VectorDiscreteFunctionType>& Pnat,
                        const std::vector<DiscreteFunctionType>& phi,
                        const std::vector<DiscreteFunctionType>& phinat,
                        const std::vector<DiscreteFunctionType>& mu,
                        const std::string& prefix,
                        const size_t step,
                        const double t,
                        const bool subsampling)
{
  auto vtk_writer = u.create_vtkwriter(u.space().grid_view(), subsampling);
  const size_t num_cells = phi.size();
  if (visualize) {
    u.add_to_vtkwriter(*vtk_writer);
    p.add_to_vtkwriter(*vtk_writer);
    for (size_t kk = 0; kk < num_cells; ++kk) {
      P[kk].add_to_vtkwriter(*vtk_writer);
      Pnat[kk].add_to_vtkwriter(*vtk_writer);
      phi[kk].add_to_vtkwriter(*vtk_writer);
      phinat[kk].add_to_vtkwriter(*vtk_writer);
      mu[kk].add_to_vtkwriter(*vtk_writer);
      phi[kk].add_gradient_to_vtkwriter(*vtk_writer);
      phinat[kk].add_gradient_to_vtkwriter(*vtk_writer);
      mu[kk].add_gradient_to_vtkwriter(*vtk_writer);
    }
  }
  std::string postfix = "_" + Dune::XT::Common::to_string(step);
  u.write_visualization(*vtk_writer, prefix + postfix);
  write_to_textfile(u, prefix, step, t);
  write_to_textfile(p, prefix, step, t);
  for (size_t kk = 0; kk < num_cells; ++kk) {
    write_to_textfile(P[kk], prefix, step, t);
    write_to_textfile(Pnat[kk], prefix, step, t);
    write_to_textfile(phi[kk], prefix, step, t);
    write_to_textfile(phinat[kk], prefix, step, t);
    write_to_textfile(mu[kk], prefix, step, t);
  }
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
               const std::vector<VectorDiscreteFunctionType>& P,
               const std::vector<VectorDiscreteFunctionType>& Pnat,
               const std::vector<DiscreteFunctionType>& phi,
               const std::vector<DiscreteFunctionType>& phinat,
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
    , C_(n_, n_, pattern_C_)
    , rhs_vector_(m_ + n_, 0.)
    , f_vector_(m_, 0., 100)
    , g_vector_(n_, 0., 0)
    , p_basis_integrated_vector_(n_)
    , dirichlet_constraints_(make_dirichlet_constraints(u_.space(), boundary_info_))
    , A_operator_(grid_view_, u_.space(), u_.space(), A_)
  {
    if (Re_ > 1e-2)
      DUNE_THROW(Dune::NotImplemented, "No Navier-Stokes solver implemented yet!");

    MatrixOperator<MatrixType, PGV, 1, 1, d> B_operator(grid_view_, p_.space(), u_.space(), B_);
    // calculate A_{ij} as \int \nabla v_i \nabla v_j
    A_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalSymmetricEllipticIntegrand<E>(1.)));
    //    A_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(1.)));
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
    A_operator_.assemble(true);
    B_operator.assemble(true);

    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
    // [A B; B^T 0] to the unit vector.
    const size_t dof_index = 0;
    pattern_C_.insert(dof_index, dof_index);
    C_ = MatrixType(n_, n_, pattern_C_);
    C_.set_entry(dof_index, dof_index, 1.);
    B_.clear_col(dof_index);
    g_vector_.set_entry(dof_index, 0.);

    dirichlet_constraints_.apply(A_, false, true);
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

    // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
    const R Fa_inv = Fa_inv_;
    const R xi = xi_;

    const size_t num_cells = phi_.size();
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells);
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phinat_local_(num_cells);
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> P_local_(num_cells);
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> Pnat_local_(num_cells);

    f_functional.append(LocalElementIntegralFunctional<E, d>(
        /*order*/ [& P_space_ = P_[0].space()](
                      const auto& test_basis,
                      const auto& param) { return 3 * P_space_.max_polorder() + test_basis.order(param); },
        /*evaluate_func*/
        [Fa_inv, xi](const auto& test_basis,
                     const DomainType& x_local,
                     DynamicVector<R>& result,
                     const XT::Common::Parameter& param) {
          const size_t sz = test_basis.size(param);
          if (result.size() < sz)
            result.resize(sz);
          std::fill(result.begin(), result.end(), 0.);
          using TestBasisType = typename LocalElementIntegralFunctional<E, d>::GenericIntegrand::LocalBasisType;
          thread_local std::vector<typename TestBasisType::RangeType> test_basis_values_;
          thread_local std::vector<typename TestBasisType::DerivativeRangeType> test_basis_grads_;
          test_basis.evaluate(x_local, test_basis_values_, param);
          test_basis.jacobians(x_local, test_basis_grads_, param);

          // evaluate P, Pnat, phi, phinat, \nabla P, \nabla Pnat, \nabla phi, div P, div Pnat and phi_tilde = (phi +
          // 1)/2 return type of the jacobians is a FieldMatrix<r, d>
          const size_t num_cells = P_local_.size();
          for (size_t kk = 0; kk < num_cells; ++kk) {
            const auto P = P_local_[kk]->evaluate(x_local, param);
            const auto Pnat = Pnat_local_[kk]->evaluate(x_local, param);
            const auto phi = phi_local_[kk]->evaluate(x_local, param)[0];
            const auto phinat = phinat_local_[kk]->evaluate(x_local, param)[0];
            const auto grad_P = P_local_[kk]->jacobian(x_local, param);
            const auto grad_phi = phi_local_[kk]->jacobian(x_local, param)[0];
            const auto phi_tilde = (phi + 1.) / 2.;

            // evaluate rhs terms
            const auto phinat_grad_phi = grad_phi * phinat;
            auto grad_P_T_times_Pnat = P;
            grad_P.mtv(Pnat, grad_P_T_times_Pnat);
            for (size_t ii = 0; ii < sz; ++ii) {
              for (size_t mm = 0; mm < d; ++mm) {
                result[ii] += (phinat_grad_phi[mm] + grad_P_T_times_Pnat[mm]) * test_basis_values_[ii][mm];
                for (size_t nn = 0; nn < d; ++nn)
                  result[ii] += (-Fa_inv * phi_tilde * P[mm] * P[nn] - 0.5 * (xi + 1) * Pnat[mm] * P[nn]
                                 - 0.5 * (xi - 1) * P[mm] * Pnat[nn])
                                * test_basis_grads_[ii][mm][nn];
              } // mm
            } // ii
          } // kk
        },
        /*post_bind_func*/
        [& phi_ = phi_, &phinat_ = phinat_, &P_ = P_, &Pnat_ = Pnat_](const E& element) {
          for (size_t kk = 0; kk < P_.size(); ++kk) {
            if (!phi_local_[kk])
              phi_local_[kk] = phi_[kk].local_function();
            if (!phinat_local_[kk])
              phinat_local_[kk] = phinat_[kk].local_function();
            if (!P_local_[kk])
              P_local_[kk] = P_[kk].local_function();
            if (!Pnat_local_[kk])
              Pnat_local_[kk] = Pnat_[kk].local_function();
            P_local_[kk]->bind(element);
            Pnat_local_[kk]->bind(element);
            phi_local_[kk]->bind(element);
            phinat_local_[kk]->bind(element);
          } // kk
        }));
    A_operator_.clear();
    A_operator_.append(f_functional);
    f_vector_ *= 0.;
    A_operator_.assemble(true);
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
  const std::vector<VectorDiscreteFunctionType>& P_;
  const std::vector<VectorDiscreteFunctionType>& Pnat_;
  const std::vector<DiscreteFunctionType>& phi_;
  const std::vector<DiscreteFunctionType>& phinat_;
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
#if HAVE_EIGEN
  using RowMajorBackendType = typename MatrixType::BackendType;
  using SolverType = ::Eigen::BiCGSTAB<RowMajorBackendType>;
  // using SolverType = ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::IncompleteLUT<R>>;
#endif // HAVE_EIGEN

  // Setup spaces and matrices and vectors
  // System is [S_{00} S_{01}; S_{10} S_{11}] [P; Pnat] = [f_{of}; g_{of}]
  // All matrices have dimension n x n, all vectors have dimension n
  // Use same pattern for all submatrices
  OfieldSolver(const VectorDiscreteFunctionType& u,
               std::vector<VectorDiscreteFunctionType>& P,
               std::vector<VectorDiscreteFunctionType>& Pnat,
               const std::vector<DiscreteFunctionType>& phi,
               const double xi,
               const double kappa,
               const double c_1,
               const double Pa,
               const double beta,
               const bool linearize)
    : u_(u)
    , P_(P)
    , Pnat_(Pnat)
    , phi_(phi)
    , xi_(xi)
    , kappa_(kappa)
    , c_1_(c_1)
    , Pa_inv_(1. / Pa)
    , beta_(beta)
    , num_cells_(P_.size())
    , P_space_(P_[0].space())
    , grid_view_(P_space_.grid_view())
    , n_(P_space_.mapper().size())
    , submatrix_pattern_(make_element_sparsity_pattern(P_space_, P_space_, grid_view_))
    , pattern_(create_pattern(n_, submatrix_pattern_))
    , S_(2 * n_, 2 * n_, pattern_, 100)
    , M_(n_, n_, submatrix_pattern_)
    , C_elliptic_part_(n_, n_, submatrix_pattern_)
    , C_linear_part_(n_, n_, submatrix_pattern_)
    , S_00_(S_, 0, n_, 0, n_)
    , S_01_(S_, 0, n_, n_, 2 * n_)
    , S_10_(S_, n_, 2 * n_, 0, n_)
    , S_11_(S_, n_, 2 * n_, n_, 2 * n_)
    , S_00_operator_(grid_view_, P_space_, P_space_, S_00_)
    , S_10_operator_(grid_view_, P_space_, P_space_, S_10_)
    , M_operator_(grid_view_, P_space_, P_space_, M_)
    , rhs_vector_(2 * n_, 0., 100)
    , old_result_(num_cells_, VectorType(2 * n_, 0.))
    , f_vector_(rhs_vector_, 0, n_)
    , g_vector_(rhs_vector_, n_, 2 * n_)
    , linearize_(linearize)
    , residual_(2 * n_, 0.)
    , x_n_(2 * n_, 0.)
    , update_(2 * n_, 0.)
    , candidate_(2 * n_, 0.)
  {
    // calculate M_{ij} as \int \psi_i phi_j
    M_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(1.)));
    M_operator_.assemble(true);
    // set S_11 = D = M
    S_11_ = M_;
    // calculate S_01 = B
    S_01_ = M_;
    S_01_ *= 1. / kappa_;

    MatrixOperator<MatrixType, PGV, d> elliptic_operator(grid_view_, P_space_, P_space_, C_elliptic_part_);
    elliptic_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-Pa_inv_)));
    elliptic_operator.assemble(true);

    solver_.analyzePattern(S_.backend());
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

  double compute_residual(const VectorType& x_n, const size_t ll, double l2_ref_P = 1., double l2_ref_Pnat = 1.)
  {
    // linear part
    S_.mv(x_n, residual_);
    // subtract rhs
    residual_ -= rhs_vector_;
    if (linearize_)
      return -1.;

    // nonlinear part
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> P_local_(num_cells_);
    VectorViewType res0_vec(residual_, 0, n_);
    VectorViewType res1_vec(residual_, n_, 2 * n_);
    const auto res0 = make_discrete_function(P_space_, res0_vec);
    const auto res1 = make_discrete_function(P_space_, res1_vec);
    auto nonlinear_res_functional = make_vector_functional(P_space_, res1_vec);
    XT::Functions::GenericGridFunction<E, d, 1> nonlinear_res_pf(
        /*order = */ 3 * P_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &P_ = P_](const E& element) {
          if (!P_local_[ll])
            P_local_[ll] = P_[ll].local_function();
          P_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, factor = -c_1_ * Pa_inv_](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto P_n = P_local_[ll]->evaluate(x_local, param);
          return P_n * (factor * (P_n * P_n));
        });
    nonlinear_res_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), nonlinear_res_pf)));
    S_00_operator_.clear();
    S_00_operator_.append(nonlinear_res_functional);
    S_00_operator_.assemble(true);
    // relative error if l2_norm is > 1, else absolute error
    l2_ref_P = l2_ref_P < 1. ? 1. : l2_ref_P;
    l2_ref_Pnat = l2_ref_Pnat < 1. ? 1. : l2_ref_Pnat;
    return l2_norm(grid_view_, res0) / l2_ref_P + l2_norm(grid_view_, res1) / l2_ref_Pnat;
  }

  void assemble_rhs(const double dt, const size_t ll)
  {
    M_.mv(P_[ll].dofs().vector(), f_vector_);
    f_vector_ /= dt;

    auto g_functional = make_vector_functional(P_space_, g_vector_);
    g_vector_ *= 0.;
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells_);
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> P_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, d> g(
        /*order = */ 3 * P_space_.max_polorder(),
        /*post_bind_func*/
        [ll, linearize_ = linearize_, &P_ = P_, &phi_ = phi_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
          if (linearize_) {
            if (!P_local_[ll])
              P_local_[ll] = P_[ll].local_function();
            P_local_[ll]->bind(element);
          }
        },
        /*evaluate_func*/
        [ll, linearize_ = linearize_, factor1 = beta_ * Pa_inv_, factor2 = -2. * c_1_ * Pa_inv_](
            const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate rhs terms
          const auto grad_phi = phi_local_[ll]->jacobian(x_local, param)[0];
          auto ret = grad_phi;
          ret *= factor1;
          if (linearize_) {
            const auto P_n = P_local_[ll]->evaluate(x_local, param);
            auto ret2 = P_n;
            ret2 *= factor2 * (P_n * P_n);
            ret += ret2;
          }
          return ret;
        });
    g_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), g)));
    S_10_operator_.clear();
    S_10_operator_.append(g_functional);
    S_10_operator_.assemble(true);
  }

  void assemble_linear_jacobian(const double dt, const size_t ll)
  {
    // assemble matrix S_{00} = M/dt + A
    S_00_operator_.clear();
    S_00_ = M_;
    S_00_ *= 1 / dt;
    // calculate A
    // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
    thread_local std::unique_ptr<VectorLocalDiscreteFunctionType> u_local_;
    const R xi = xi_;
    XT::Functions::GenericGridFunction<E, d, d> Omega_minus_xi_D_transposed(
        /*order = */ std::max(u_.space().max_polorder() - 1, 0),
        /*post_bind_func*/
        [& u_ = u_](const E& element) {
          if (!u_local_)
            u_local_ = u_.local_function();
          u_local_->bind(element);
        },
        /*evaluate_func*/
        [xi](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate \nabla u
          auto grad_u = u_local_->jacobian(x_local, param);
          auto grad_u_T = grad_u;
          grad_u_T.transpose();
          auto& ret = grad_u;
          ret *= (1. - xi) / 2.;
          grad_u_T *= (1. + xi) / 2.;
          ret -= grad_u_T;
          return ret;
        });
    S_00_operator_.append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(Omega_minus_xi_D_transposed)));
    S_00_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(u_)));
    S_00_operator_.assemble(true);

    // calculate linear part S_10 = C
    S_10_operator_.clear();
    S_10_ = C_elliptic_part_;
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_inv_phi(
        /*order = */ P_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &phi_ = phi_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, factor = c_1_ * Pa_inv_](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi = phi_local_[ll]->evaluate(x_local, param);
          return factor * phi;
        });
    S_10_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(c1_Pa_inv_phi)));
    S_10_operator_.assemble(true);
    C_linear_part_ = S_10_;

    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_nonlinear_jacobian(ll);
  }

  void assemble_nonlinear_jacobian(const size_t ll)
  {
    S_10_operator_.clear();
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> P_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_P2(
        /*order = */ 2. * P_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &P_ = P_](const E& element) {
          if (!P_local_[ll])
            P_local_[ll] = P_[ll].local_function();
          P_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, factor = -c_1_ * Pa_inv_](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = P_local_[ll]->evaluate(x_local, param);
          return factor * P_n.two_norm2();
        });
    S_10_operator_.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(c1_Pa_P2)));
    XT::Functions::GenericGridFunction<E, d, d> minus_two_frac_c1_Pa_Pn_otimes_Pn(
        /*order = */ 2 * P_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &P_ = P_](const E& element) {
          if (!P_local_[ll])
            P_local_[ll] = P_[ll].local_function();
          P_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, factor = -2. * c_1_ * Pa_inv_](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = P_local_[ll]->evaluate(x_local, param);
          FieldMatrix<R, d, d> ret;
          for (size_t ii = 0; ii < d; ++ii)
            for (size_t jj = 0; jj < d; ++jj)
              ret[ii][jj] = P_n[ii] * P_n[jj];
          ret *= factor;
          return ret;
        });
    // Pn_otimes_Pn is symmetric, so no need to transpose
    S_10_operator_.append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(minus_two_frac_c1_Pa_Pn_otimes_Pn)));
    S_10_operator_.assemble(true);
  }

  void copy_P_Pnat_to(VectorType& vec, const size_t ll)
  {
    for (size_t ii = 0; ii < n_; ++ii) {
      vec[ii] = P_[ll].dofs().vector()[ii];
      vec[n_ + ii] = Pnat_[ll].dofs().vector()[ii];
    }
  }

  void fill_P_Pnat_from(const VectorType& vec, const size_t ll)
  {
    for (size_t ii = 0; ii < n_; ++ii) {
      P_[ll].dofs().vector()[ii] = vec[ii];
      Pnat_[ll].dofs().vector()[ii] = vec[n_ + ii];
    }
  }

  void revert_jacobian_to_linear()
  {
    S_10_ = C_linear_part_;
  }

  void solve_linear_system(const size_t iter, const size_t ll)
  {
    const auto begin = std::chrono::steady_clock::now();
    //    std::ofstream S_file("S_" + XT::Common::to_string(dt) + ".txt");
    //    S_file << S_ << std::endl;
    //    S_file.close();
    //    DUNE_THROW(Dune::NotImplemented, "");
    //    const auto ret = XT::LA::solve(S_, rhs_vector_, XT::LA::SolverOptions<MatrixType>::options("lu.umfpack"));
    residual_ *= -1.;
    //      update_ = XT::LA::solve(S_, residual_);
    solver_.compute(S_.backend());
    update_.backend() = solver_.solveWithGuess(residual_.backend(), old_result_[ll].backend());
    old_result_[ll] = update_;
    const std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Ofield in iteration " << iter << " took: " << time.count() << " s!" << std::endl;
  }

  void apply(const double dt, const size_t ll)
  {
    auto begin = std::chrono::steady_clock::now();
    // *********** assemble linear part of jacobian **********
    begin = std::chrono::steady_clock::now();
    assemble_linear_jacobian(dt, ll);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling linear part of jacobian took: " << time.count() << " s!" << std::endl;

    // ************ create rhs = (f_{pf}, g_{pf}, h_{pf}) *************
    begin = std::chrono::steady_clock::now();
    assemble_rhs(dt, ll);
    time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling rhs took: " << time.count() << " s!" << std::endl;

    // fill x_n
    copy_P_Pnat_to(x_n_, ll);

    if (linearize_) {
      compute_residual(x_n_, ll);
      solve_linear_system(0, ll);
      x_n_ += update_;
      fill_P_Pnat_from(x_n_, ll);
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 1000;
      const auto max_dampening_iter = 1000;

      const auto l2_ref_P = l2_norm(grid_view_, P_[ll]);
      const auto l2_ref_Pnat = l2_norm(grid_view_, Pnat_[ll]);

      size_t iter = 0;
      while (true) {

        // ********* compute residual *********
        begin = std::chrono::steady_clock::now();
        const auto res = compute_residual(x_n_, ll, l2_ref_P, l2_ref_Pnat);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

        if (res < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_nonlinear_jacobian(ll);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        solve_linear_system(iter, ll);

        DUNE_THROW_IF(iter >= max_iter, Exceptions::operator_error, "max iterations reached!\n|residual|_l2 = " << res);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        revert_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res << "\niter = " << iter << "\n");
          candidate_ = x_n_ + update_ * lambda;
          fill_P_Pnat_from(candidate_, ll);
          candidate_res = compute_residual(candidate_, ll, l2_ref_P, l2_ref_Pnat);
          std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        std::cout << "Current res: " << candidate_res << std::endl;
        x_n_ = candidate_;
        iter += 1;
      } // while (true)
    }
  }

  const VectorDiscreteFunctionType& u_;
  std::vector<VectorDiscreteFunctionType>& P_;
  std::vector<VectorDiscreteFunctionType>& Pnat_;
  const std::vector<DiscreteFunctionType>& phi_;
  const R xi_;
  const R kappa_;
  const R c_1_;
  const R Pa_inv_;
  const R beta_;
  const size_t num_cells_;
  const SpaceInterface<PGV, d, 1, R>& P_space_;
  const PGV& grid_view_;
  const size_t n_;
  const XT::LA::SparsityPatternDefault submatrix_pattern_;
  const XT::LA::SparsityPatternDefault pattern_;
  MatrixType S_;
  MatrixType M_;
  MatrixType C_elliptic_part_;
  MatrixType C_linear_part_;
  MatrixViewType S_00_;
  MatrixViewType S_01_;
  MatrixViewType S_10_;
  MatrixViewType S_11_;
  MatrixOperator<MatrixViewType, PGV, d> S_00_operator_;
  MatrixOperator<MatrixViewType, PGV, d> S_10_operator_;
  MatrixOperator<MatrixType, PGV, d> M_operator_;
  VectorType rhs_vector_;
  std::vector<VectorType> old_result_;
  XT::LA::VectorView<VectorType> f_vector_;
  XT::LA::VectorView<VectorType> g_vector_;
  SolverType solver_;
  const bool linearize_;
  VectorType residual_;
  VectorType x_n_;
  VectorType update_;
  VectorType candidate_;
};

struct PfieldSolver
{
#if HAVE_EIGEN
  using RowMajorBackendType = typename MatrixType::BackendType;
  using SolverType = ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::IncompleteLUT<R>>;
#endif // HAVE_EIGEN

  // Setup spaces and matrices and vectors
  // System is [S_{00} 0 S_{02}; S_{10} S_{11} 0; S_{20} S_{21} S_{22}] [phi; phinat; mu] = [f_{pf}; g_{pf}; h_{pf}]
  // All matrices have dimension n x n, all vectors have dimension n
  PfieldSolver(const VectorDiscreteFunctionType& u,
               const std::vector<VectorDiscreteFunctionType>& P,
               std::vector<DiscreteFunctionType>& phi,
               std::vector<DiscreteFunctionType>& phinat,
               std::vector<DiscreteFunctionType>& mu,
               const double gamma,
               const double c_1,
               const double Pa,
               const double Be,
               const double Ca,
               const double beta,
               const double epsilon,
               const double In,
               const XT::Grid::BoundaryInfo<PI>& boundary_info,
               const bool linearize)
    : u_(u)
    , P_(P)
    , phi_(phi)
    , phinat_(phinat)
    , mu_(mu)
    , gamma_(gamma)
    , c_1_(c_1)
    , Pa_(Pa)
    , Be_(Be)
    , Ca_(Ca)
    , beta_(beta)
    , epsilon_(epsilon)
    , In_(In)
    , boundary_info_(boundary_info)
    , linearize_(linearize)
    , num_cells_(phi_.size())
    , phi_space_(phi_[0].space())
    , grid_view_(phi_space_.grid_view())
    , n_(phi_space_.mapper().size())
    , submatrix_pattern_(make_element_sparsity_pattern(phi_space_, phi_space_, grid_view_))
    , pattern_(create_pattern(n_, submatrix_pattern_))
    , S_(3 * n_, 3 * n_, pattern_, 100)
    , M_(n_, n_, submatrix_pattern_)
    , elliptic_matrix_(n_, n_, submatrix_pattern_)
    , A_linear_part_(n_, n_, submatrix_pattern_)
    , A_nonlinear_part_(n_, n_, submatrix_pattern_)
    , J_linear_part_(n_, n_, submatrix_pattern_)
    , S_00_(S_, 0, n_, 0, n_)
    , S_01_(S_, 0, n_, n_, 2 * n_)
    , S_10_(S_, n_, 2 * n_, 0, n_)
    , S_11_(S_, n_, 2 * n_, n_, 2 * n_)
    , S_12_(S_, n_, 2 * n_, 2 * n_, 3 * n_)
    , S_20_(S_, 2 * n_, 3 * n_, 0, n_)
    , S_22_(S_, 2 * n_, 3 * n_, 2 * n_, 3 * n_)
    , S_00_operator_(grid_view_, phi_space_, phi_space_, S_00_)
    , S_10_operator_(grid_view_, phi_space_, phi_space_, S_10_)
    , A_nonlinear_part_operator_(grid_view_, phi_space_, phi_space_, A_nonlinear_part_)
    , rhs_vector_(3 * n_, 0., 100)
    , old_result_(num_cells_, VectorType(3 * n_, 0.))
    , f_vector_(rhs_vector_, 2 * n_, 3 * n_)
    , g_vector_(rhs_vector_, 0, n_)
    , h_vector_(rhs_vector_, n_, 2 * n_)
    , dirichlet_constraints_(make_dirichlet_constraints(phi_space_, boundary_info_))
    , residual_(3 * n_, 0.)
    , x_n_(3 * n_, 0.)
    , update_(3 * n_, 0.)
    , candidate_(3 * n_, 0.)
  {
    assert(phinat_[0].space().mapper().size() == n_);
    assert(mu_[0].space().mapper().size() == n_);
    MatrixOperator<MatrixType, PGV, 1> M_operator(grid_view_, phi_space_, phi_space_, M_);
    M_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
    MatrixOperator<MatrixType, PGV, 1> elliptic_operator(grid_view_, phi_space_, phi_space_, elliptic_matrix_);
    elliptic_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
    M_operator.append(dirichlet_constraints_);
    M_operator.assemble(true);
    elliptic_operator.assemble(true);
    A_linear_part_ = elliptic_matrix_ * epsilon_;
    J_linear_part_ = elliptic_matrix_ * 1. / Be_;
    J_linear_part_ += M_ * 1. / Ca_;
    // Set matrix S_{22} = C = M
    S_22_ = M_;
    // Set matrix S_{11} = H = M
    S_11_ = M_;
    // Set matrix S_{01} = E
    S_01_ = elliptic_matrix_;
    S_01_ *= gamma_;

    // apply Dirichlet constraints to linear part
    for (const auto& DoF : dirichlet_constraints_.dirichlet_DoFs()) {
      A_linear_part_.clear_row(DoF);
      S_22_.clear_row(DoF);
    }

    solver_.analyzePattern(S_.backend());
  }

  static XT::LA::SparsityPatternDefault create_pattern(const size_t n,
                                                       const XT::LA::SparsityPatternDefault& submatrix_pattern)
  {
    // Use same pattern for all submatrices
    XT::LA::SparsityPatternDefault pattern(3 * n);
    for (size_t ii = 0; ii < n; ++ii)
      for (const auto& jj : submatrix_pattern.inner(ii)) {
        pattern.insert(ii, jj); // S_{00}
        pattern.insert(ii, n + jj); // S_{01}
        pattern.insert(n + ii, jj); // S_{10}
        pattern.insert(n + ii, n + jj); // S_{11}
        pattern.insert(n + ii, 2 * n + jj); // S_{12}
        pattern.insert(2 * n + ii, jj); // S_{20}
        pattern.insert(2 * n + ii, 2 * n + jj); // S_{22}
      }
    pattern.sort();
    return pattern;
  }

  double compute_residual(
      const VectorType& x_n, const size_t ll, double l2_ref_phi = 1., double l2_ref_phinat = 1., double l2_ref_mu = 1.)
  {
    VectorViewType res2_vec(residual_, 2 * n_, 3 * n_);
    // linear part
    S_.mv(x_n, residual_);
    // subtract rhs
    residual_ -= rhs_vector_;
    if (linearize_) {
      dirichlet_constraints_.apply(res2_vec);
      return 10.;
    }

    // nonlinear part
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells_);
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> mu_local_(num_cells_);
    VectorViewType res0_vec(residual_, 0, n_);
    VectorViewType res1_vec(residual_, n_, 2 * n_);
    const auto res0 = make_discrete_function(phi_space_, res0_vec);
    const auto res1 = make_discrete_function(phi_space_, res1_vec);
    const auto res2 = make_discrete_function(phi_space_, res2_vec);
    auto nonlinear_res1_functional = make_vector_functional(phi_space_, res1_vec);
    auto nonlinear_res2_functional = make_vector_functional(phi_space_, res2_vec);
    const auto Bfunc =
        [epsilon_inv = 1. / epsilon_](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi_n = phi_local_[kk]->evaluate(x_local, param)[0];
          return epsilon_inv * std::pow(std::pow(phi_n, 2) - 1, 2);
        };
    const auto wfunc = [](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
      const auto phi_n = phi_local_[kk]->evaluate(x_local, param)[0];
      if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
        return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
      else
        return 0.;
    };
    XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf1(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, num_cells = num_cells_, &phi_ = phi_, &mu_ = mu_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
          if (!mu_local_[ll])
            mu_local_[ll] = mu_[ll].local_function();
          mu_local_[ll]->bind(element);
          if (num_cells > 1) {
            for (size_t kk = 0; kk < num_cells; ++kk) {
              if (!phi_local_[kk])
                phi_local_[kk] = phi_[kk].local_function();
              phi_local_[kk]->bind(element);
            }
          }
        },
        /*evaluate_func*/
        [wfunc,
         Bfunc,
         ll,
         In_inv = 1. / In_,
         eps_inv = 1. / epsilon_,
         num_cells = num_cells_,
         inv_Be_eps2 = 1. / (Be_ * std::pow(epsilon_, 2))](const DomainType& x_local,
                                                           const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto phi_n = phi_local_[ll]->evaluate(x_local, param);
          const auto mu_n = mu_local_[ll]->evaluate(x_local, param);
          auto ret = inv_Be_eps2 * (3. * phi_n * phi_n - 1) * mu_n;
          if (num_cells > 1) {
            R wsum = 0.;
            R Bsum = 0.;
            for (size_t kk = 0; kk < num_cells; ++kk) {
              if (kk != ll) {
                wsum += wfunc(kk, x_local, param);
                Bsum += Bfunc(kk, x_local, param);
              }
            } // kk
            ret += In_inv * 4 * eps_inv * (std::pow(phi_n, 3) - phi_n) * wsum;
            auto w_prime = 0;
            if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.)) {
              const auto ln = std::log((1 + phi_n) / (1 - phi_n));
              const auto ln2 = std::pow(ln, 2);
              w_prime = 2 * std::exp(-0.5 * ln2) * ln / (std::pow(phi_n, 2) - 1);
            }
            ret += In_inv * w_prime * Bsum;
          } // num_cells > 1
          return ret;
        });
    XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf2(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &phi_ = phi_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, inv_eps = 1. / epsilon_](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto phi_n = phi_local_[ll]->evaluate(x_local, param);
          return inv_eps * (phi_n * phi_n - 1) * phi_n;
        });
    nonlinear_res1_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), nonlinear_res_pf1)));
    nonlinear_res2_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), nonlinear_res_pf2)));
    S_00_operator_.clear();
    S_00_operator_.append(nonlinear_res1_functional);
    S_00_operator_.append(nonlinear_res2_functional);
    S_00_operator_.assemble(true);
    dirichlet_constraints_.apply(res2_vec);
    // relative error if l2_norm is > 1, else absolute error
    l2_ref_phi = l2_ref_phi < 1. ? 1. : l2_ref_phi;
    l2_ref_phinat = l2_ref_phinat < 1. ? 1. : l2_ref_phinat;
    l2_ref_mu = l2_ref_mu < 1. ? 1. : l2_ref_mu;
    return l2_norm(grid_view_, res0) / l2_ref_phi + l2_norm(grid_view_, res1) / l2_ref_phinat
           + l2_norm(grid_view_, res2) / l2_ref_mu;
  }

  void assemble_rhs(const double dt, const size_t ll)
  {
    S_00_operator_.clear();
    auto f_functional = make_vector_functional(phi_space_, f_vector_);
    auto h_functional = make_vector_functional(phi_space_, h_vector_);

    // calculate f
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells_);
    if (linearize_)
      f_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &phi_ = phi_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll, two_epsilon_inv = 2. / epsilon_](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate phi_
          const auto phi_n = phi_local_[ll]->evaluate(x_local, param);
          return two_epsilon_inv * std::pow(phi_n, 3);
        });
    f_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), f_pf)));
    if (linearize_)
      S_00_operator_.append(f_functional);

    // calculate g
    M_.mv(phi_[ll].dofs().vector(), g_vector_);
    g_vector_ /= dt;

    // calculate h
    h_vector_ *= 0.;
    thread_local std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>> P_local_(num_cells_);
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> mu_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, linearize_ = linearize_, &phi_ = phi_, &mu_ = mu_, &P_ = P_](const E& element) {
          if (!P_local_[ll])
            P_local_[ll] = P_[ll].local_function();
          P_local_[ll]->bind(element);
          if (linearize_) {
            if (!phi_local_[ll])
              phi_local_[ll] = phi_[ll].local_function();
            phi_local_[ll]->bind(element);
            if (!mu_local_[ll])
              mu_local_[ll] = mu_[ll].local_function();
            mu_local_[ll]->bind(element);
          }
        },
        /*evaluate_func*/
        [ll,
         linearize_ = linearize_,
         factor0 = 6. / (Be_ * std::pow(epsilon_, 2)),
         factor1 = -c_1_ / (2. * Pa_),
         factor2 = -beta_ / Pa_](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto Pn = P_local_[ll]->evaluate(x_local, param);
          const auto grad_P = P_local_[ll]->jacobian(x_local, param);
          R div_P(0.);
          for (size_t ii = 0; ii < d; ++ii)
            div_P += grad_P[ii][ii];
          auto ret = factor1 * (Pn * Pn) + factor2 * div_P;
          if (linearize_) {
            const auto phi_n = phi_local_[ll]->evaluate(x_local, param);
            const auto mu_n = mu_local_[ll]->evaluate(x_local, param);
            ret += factor0 * std::pow(phi_n, 2) * mu_n;
          }
          return ret;
        });
    h_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(1.), h_pf)));
    S_00_operator_.append(h_functional);

    // assemble rhs
    S_00_operator_.assemble(true);
  }

  void assemble_linear_jacobian(const double dt, const size_t ll)
  {
    // assemble matrix S_{00} = M/dt + D
    S_00_operator_.clear();
    S_00_ = M_;
    S_00_ *= 1. / dt;
    thread_local std::unique_ptr<VectorLocalDiscreteFunctionType> u_local_;
    XT::Functions::GenericGridFunction<E, d, 1> minus_u(
        /*order = */ u_.space().max_polorder(),
        /*post_bind_func*/
        [& u_ = u_](const E& element) {
          if (!u_local_)
            u_local_ = u_.local_function();
          u_local_->bind(element);
        },
        /*evaluate_func*/
        [](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = u_local_->evaluate(x_local, param);
          ret *= -1.;
          return ret;
        });
    S_00_operator_.append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u)));
    S_00_operator_.assemble(true);
    // linear part of matrix S_{12} = J
    S_12_ = J_linear_part_;
    // linear part of matrix S_{20} = A
    S_20_ = A_linear_part_;

    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_nonlinear_jacobian(ll);
  }

  void assemble_nonlinear_jacobian(const size_t ll)
  {
    A_nonlinear_part_operator_.clear();
    A_nonlinear_part_ *= 0.;
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> phi_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, 1, 1> A_nonlinear_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, &phi_ = phi_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
        },
        /*evaluate_func*/
        [ll](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi_n = phi_local_[ll]->evaluate(x_local, param);
          return (3. * phi_n * phi_n - 1.);
        });
    A_nonlinear_part_operator_.append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(A_nonlinear_prefactor)));
    A_nonlinear_part_operator_.assemble(true);
    A_nonlinear_part_ *= 1. / epsilon_;
    S_20_ += A_nonlinear_part_;
    A_nonlinear_part_ *= 1. / (Be_ * epsilon_);
    S_12_ += A_nonlinear_part_;

    // assemble matrix S_{10} = G
    S_10_operator_.clear();
    S_10_ *= 0.;
    const auto Bfunc =
        [epsilon_inv = 1. / epsilon_](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi_n = phi_local_[kk]->evaluate(x_local, param)[0];
          return epsilon_inv * std::pow(std::pow(phi_n, 2) - 1, 2);
        };
    const auto wfunc = [](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
      const auto phi_n = phi_local_[kk]->evaluate(x_local, param)[0];
      if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
        return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
      else
        return 0.;
    };
    thread_local std::vector<std::unique_ptr<LocalDiscreteFunctionType>> mu_local_(num_cells_);
    XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, num_cells = num_cells_, &phi_ = phi_, &mu_ = mu_](const E& element) {
          if (!phi_local_[ll])
            phi_local_[ll] = phi_[ll].local_function();
          phi_local_[ll]->bind(element);
          if (!mu_local_[ll])
            mu_local_[ll] = mu_[ll].local_function();
          mu_local_[ll]->bind(element);
          if (num_cells > 1) {
            for (size_t kk = 0; kk < num_cells; ++kk) {
              if (!phi_local_[kk])
                phi_local_[kk] = phi_[kk].local_function();
              phi_local_[kk]->bind(element);
            }
          }
        },
        /*evaluate_func*/
        [ll,
         In_inv = 1. / In_,
         num_cells = num_cells_,
         eps_inv = 1. / epsilon_,
         six_inv_Be_eps2 = 6. / (Be_ * std::pow(epsilon_, 2)),
         &Bfunc,
         &wfunc](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi_n = phi_local_[ll]->evaluate(x_local, param)[0];
          const auto mu_n = mu_local_[ll]->evaluate(x_local, param)[0];
          auto ret = six_inv_Be_eps2 * phi_n * mu_n;
          if (num_cells > 1) {
            R wsum = 0.;
            R Bsum = 0.;
            for (size_t kk = 0; kk < num_cells; ++kk) {
              if (kk != ll) {
                wsum += wfunc(kk, x_local, param);
                Bsum += Bfunc(kk, x_local, param);
              }
            } // kk
            ret += In_inv * 4 * eps_inv * (3. * std::pow(phi_n, 2) - 1) * wsum;
            auto w_twoprime = 0;
            if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.)) {
              const auto ln = std::log((1 + phi_n) / (1 - phi_n));
              const auto ln2 = std::pow(ln, 2);
              w_twoprime = 4 * std::exp(-0.5 * ln2) * (ln2 - phi_n * ln - 1) / (std::pow(std::pow(phi_n, 2) - 1, 2));
            }
            ret += In_inv * w_twoprime * Bsum;
          } // num_cells > 1
          return ret;
        });
    S_10_operator_.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(G_prefactor)));
    S_10_operator_.assemble(true);

    // clear row, for computation of residual we do not want the unit row
    for (const auto& DoF : dirichlet_constraints_.dirichlet_DoFs())
      S_20_.unit_row(DoF);
  }

  void copy_phi_phinat_mu_to(VectorType& vec, const size_t ll)
  {
    for (size_t ii = 0; ii < n_; ++ii) {
      vec[ii] = phi_[ll].dofs().vector()[ii];
      vec[n_ + ii] = phinat_[ll].dofs().vector()[ii];
      vec[2 * n_ + ii] = mu_[ll].dofs().vector()[ii];
    }
  }

  void fill_phi_phinat_mu_from(const VectorType& vec, const size_t ll)
  {
    for (size_t ii = 0; ii < n_; ++ii) {
      phi_[ll].dofs().vector()[ii] = vec[ii];
      phinat_[ll].dofs().vector()[ii] = vec[n_ + ii];
      mu_[ll].dofs().vector()[ii] = vec[2 * n_ + ii];
    }
  }

  void revert_jacobian_to_linear()
  {
    // clear S_{10} = G
    S_10_ *= 0.;
    // linear part of matrix S_{12} = J
    S_12_ = J_linear_part_;
    // linear part of matrix S_{20} = A
    S_20_ = A_linear_part_;
  }

  void solve_linear_system(const size_t iter, const size_t ll)
  {
    const auto begin = std::chrono::steady_clock::now();
    //    std::ofstream S_file("S_" + XT::Common::to_string(dt) + ".txt");
    //    S_file << S_ << std::endl;
    //    S_file.close();
    //    DUNE_THROW(Dune::NotImplemented, "");
    //    const auto ret = XT::LA::solve(S_, rhs_vector_, XT::LA::SolverOptions<MatrixType>::options("lu.umfpack"));
    residual_ *= -1.;
    //      update_ = XT::LA::solve(S_, residual_);
    solver_.compute(S_.backend());
    update_.backend() = solver_.solveWithGuess(residual_.backend(), old_result_[ll].backend());
    old_result_[ll] = update_;
    const std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Pfield in iteration " << iter << " took: " << time.count() << " s!" << std::endl;
  }

  void apply(const double dt, const size_t ll)
  {
    auto begin = std::chrono::steady_clock::now();
    // *********** assemble linear part of jacobian **********
    begin = std::chrono::steady_clock::now();
    assemble_linear_jacobian(dt, ll);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling linear part of jacobian took: " << time.count() << " s!" << std::endl;

    // ************ create rhs = (f_{pf}, g_{pf}, h_{pf}) *************
    begin = std::chrono::steady_clock::now();
    assemble_rhs(dt, ll);
    time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling rhs took: " << time.count() << " s!" << std::endl;

    // fill x_n
    copy_phi_phinat_mu_to(x_n_, ll);

    if (linearize_) {
      compute_residual(x_n_, ll);
      solve_linear_system(0, ll);
      x_n_ += update_;
      fill_phi_phinat_mu_from(x_n_, ll);
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 1000;
      const auto max_dampening_iter = 1000;

      const auto l2_ref_phi = l2_norm(grid_view_, phi_[ll]);
      const auto l2_ref_phinat = l2_norm(grid_view_, phinat_[ll]);
      const auto l2_ref_mu = l2_norm(grid_view_, mu_[ll]);

      size_t iter = 0;
      while (true) {

        // ********* compute residual *********
        begin = std::chrono::steady_clock::now();
        const auto res = compute_residual(x_n_, ll, l2_ref_phi, l2_ref_phinat, l2_ref_mu);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

        if (res < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_nonlinear_jacobian(ll);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        solve_linear_system(iter, ll);

        DUNE_THROW_IF(iter >= max_iter, Exceptions::operator_error, "max iterations reached!\n|residual|_l2 = " << res);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        revert_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res << "\nl = " << iter << "\n");
          candidate_ = x_n_ + update_ * lambda;
          fill_phi_phinat_mu_from(candidate_, ll);
          candidate_res = compute_residual(candidate_, ll, l2_ref_phi, l2_ref_phinat, l2_ref_mu);
          std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        std::cout << "Current res: " << candidate_res << std::endl;
        x_n_ = candidate_;
        iter += 1;
      } // while (true)
    }
  }

  const VectorDiscreteFunctionType& u_;
  const std::vector<VectorDiscreteFunctionType>& P_;
  std::vector<DiscreteFunctionType>& phi_;
  std::vector<DiscreteFunctionType>& phinat_;
  std::vector<DiscreteFunctionType>& mu_;
  const double gamma_;
  const double c_1_;
  const double Pa_;
  const double Be_;
  const double Ca_;
  const double beta_;
  const double epsilon_;
  const double In_;
  const XT::Grid::BoundaryInfo<PI>& boundary_info_;
  const bool linearize_;
  const size_t num_cells_;
  const SpaceInterface<PGV, 1, 1, R>& phi_space_;
  const PGV& grid_view_;
  const size_t n_;
  const XT::LA::SparsityPatternDefault submatrix_pattern_;
  const XT::LA::SparsityPatternDefault pattern_;
  MatrixType S_;
  MatrixType M_;
  MatrixType elliptic_matrix_;
  MatrixType A_linear_part_;
  MatrixType A_nonlinear_part_;
  MatrixType J_linear_part_;
  MatrixViewType S_00_;
  MatrixViewType S_01_;
  MatrixViewType S_10_;
  MatrixViewType S_11_;
  MatrixViewType S_12_;
  MatrixViewType S_20_;
  MatrixViewType S_22_;
  MatrixOperator<MatrixViewType, PGV, 1> S_00_operator_;
  MatrixOperator<MatrixViewType, PGV, 1> S_10_operator_;
  MatrixOperator<MatrixType, PGV, 1> A_nonlinear_part_operator_;
  VectorType rhs_vector_;
  std::vector<VectorType> old_result_;
  XT::LA::VectorView<VectorType> f_vector_;
  XT::LA::VectorView<VectorType> g_vector_;
  XT::LA::VectorView<VectorType> h_vector_;
  DirichletConstraints<PI, SpaceInterface<PGV, 1, 1, R>> dirichlet_constraints_;
  SolverType solver_;
  VectorType residual_;
  VectorType x_n_;
  VectorType update_;
  VectorType candidate_;
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


    // get testcase
    auto testcase = config.template get<std::string>("problem.testcase");

    // grid config
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(16));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(4));

    // timestepping
    double t_end = config.template get<double>("fem.t_end", 340.);
    double dt = config.template get<double>("fem.dt", 0.005);
    const bool linearize = config.template get<bool>("problem.linearize", false);
    std::cout << "linearize: " << linearize << std::endl;

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
    double In = config.template get<double>("problem.In", 1.);
    double Re = rho * U * L / eta;
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

    // output
    std::string filename = config.get("output.filename", "drosophila") + (linearize ? "_linearized" : "");
    bool subsampling = config.get<bool>("output.subsampling", true);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);

    // create grid for [0, 160] x [0, 40] with periodic boundaries in x-direction
    FieldVector<double, d> lower_left, upper_right;
    std::string periodic_dirs;
    size_t num_cells;
    if (testcase == "single_cell") {
      lower_left = {{0., 0.}};
      upper_right = {{160., 40.}};
      periodic_dirs = "01";
      num_cells = 1;
    } else if (testcase == "two_cells") {
      lower_left = {{0., 0.}};
      upper_right = {{50., 50.}};
      // periodic_dirs = "11";
      periodic_dirs = "00";
      num_cells = 2;
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    }
    auto grid = XT::Grid::make_cube_grid<G>(lower_left, upper_right, /*num_elements=*/{num_elements_x, num_elements_y});
    grid.grid().globalRefine(1);
    const double vol_domain = (upper_right[0] - lower_left[0]) * (upper_right[1] - lower_left[1]);
    auto nonperiodic_grid_view = grid.leaf_view();
    std::bitset<d> periodic_directions(periodic_dirs);
    PGV grid_view(nonperiodic_grid_view, periodic_directions);

    // create spaces for the variables
    // Taylor-Hood P2-P1 space for the Stokes variables (u, p) \in S_{stokes} = (H^1_0)^d x L^2\R, test functions (v, q)
    // \in S_{stokes}
    auto u_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    auto p_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/1);

    // P2 spaces for the orientation field variables P, P^{\natural} \in S_{ofield} = (H^1)^d x (L^2)^d, test functions
    // (Q, Q^\natural) \in S_{ofield}
    auto P_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);

    // P2 spaces for the phase field variables \phi, \phi^\natural, \mu \in S_{pfield} = H^1_0 x H^1 x H^1, test
    // functions (\chi, \chi^\natural, \nu) \in S_{pfield}
    auto phi_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);

    // create discrete functions for the variables
    auto u = make_discrete_function<VectorType>(u_space, "u");
    auto p = make_discrete_function<VectorType>(p_space, "p");
    std::vector<VectorDiscreteFunctionType> P;
    std::vector<VectorDiscreteFunctionType> Pnat;
    std::vector<DiscreteFunctionType> phi;
    std::vector<DiscreteFunctionType> phinat;
    std::vector<DiscreteFunctionType> mu;
    for (size_t ii = 0; ii < num_cells; ++ii) {
      const auto ii_str = XT::Common::to_string(ii);
      P.emplace_back(make_discrete_function<VectorType>(P_space, "P_" + ii_str));
      Pnat.emplace_back(make_discrete_function<VectorType>(P_space, "Pnat_" + ii_str));
      phi.emplace_back(make_discrete_function<VectorType>(phi_space, "phi_" + ii_str));
      phinat.emplace_back(make_discrete_function<VectorType>(phi_space, "phinat_" + ii_str));
      mu.emplace_back(make_discrete_function<VectorType>(phi_space, "mu_" + ii_str));
    }

    if (testcase == "single_cell") {
      // create and project initial values
      // we only need initial values for P and phi
      // mu_initial is only needed if linearization is used
      // Initially, cell is circular with Radius R=5 and placed in the center of the domain
      // \Omega = [0, 160] \times [0, 40].
      // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
      // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
      FieldVector<double, d> center{upper_right[0] / 2., upper_right[1] / 2.};
      auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
      const XT::Functions::GenericFunction<d> phi_initial(
          50,
          /*evaluate=*/
          [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2.) * epsilon)); },
          /*name=*/"phi_initial");
      const XT::Functions::GenericFunction<d> mu_initial(50,
                                                         /*evaluate=*/
                                                         [phi_initial, epsilon](const auto& x, const auto& param) {
                                                           // TODO: add approximation of laplacian term
                                                           const auto phi = phi_initial.evaluate(x, param);
                                                           return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                         },
                                                         /*name=*/"mu_initial");

      // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
      // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
      std::srand(1); // set seed for std::rand to 1
      const XT::Functions::GenericFunction<d, d> P_initial(50,
                                                           /*evaluate=*/
                                                           [phi_initial](const auto& x, const auto& param) {
                                                             // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                             // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                             // auto ret = FieldVector<double, d>({1. + rand1, 0. +
                                                             // rand2});
                                                             auto ret = FieldVector<double, d>({1., 0.});
                                                             ret *= (phi_initial.evaluate(x, param) + 1.) / 2.;
                                                             return ret;
                                                           },
                                                           /*name=*/"P_initial");
      const XT::Functions::ConstantFunction<d, d> u_initial(0.);

      // interpolate initial and boundary values
      default_interpolation(phi_initial, phi[0]);
      default_interpolation(mu_initial, mu[0]);
      default_interpolation(P_initial, P[0]);
      default_interpolation(u_initial, u);
    } else if (testcase == "two_cells") {
      // create and project initial values
      // we only need initial values for P_i and phi_i
      // mu_initial is only needed if linearization is used
      // Initially, cell is circular with Radius R=5 and placed in the center of the domain
      // \Omega = [0, 160] \times [0, 40].
      // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
      // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
      FieldVector<double, d> center1{15, 15};
      FieldVector<double, d> center2{35, 35};
      auto r1 = [center1](const auto& xr) { return 4.0 - (center1 - xr).two_norm(); };
      auto r2 = [center2](const auto& xr) { return 4.0 - (center2 - xr).two_norm(); };
      const XT::Functions::GenericFunction<d> phi1_initial(50,
                                                           /*evaluate=*/
                                                           [r = r1, epsilon](const auto& x, const auto& /*param*/) {
                                                             return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
                                                           },
                                                           /*name=*/"phi1_initial");
      const XT::Functions::GenericFunction<d> mu1_initial(50,
                                                          /*evaluate=*/
                                                          [phi1_initial, epsilon](const auto& x, const auto& param) {
                                                            // TODO: add approximation of laplacian term
                                                            const auto phi = phi1_initial.evaluate(x, param);
                                                            return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                          },
                                                          /*name=*/"mu1_initial");
      const XT::Functions::GenericFunction<d> phi2_initial(50,
                                                           /*evaluate=*/
                                                           [r = r2, epsilon](const auto& x, const auto& /*param*/) {
                                                             return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
                                                           },
                                                           /*name=*/"phi1_initial");
      const XT::Functions::GenericFunction<d> mu2_initial(50,
                                                          /*evaluate=*/
                                                          [phi2_initial, epsilon](const auto& x, const auto& param) {
                                                            // TODO: add approximation of laplacian term
                                                            const auto phi = phi2_initial.evaluate(x, param);
                                                            return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                          },
                                                          /*name=*/"mu1_initial");

      // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
      // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
      std::srand(1); // set seed for std::rand to 1
      const XT::Functions::GenericFunction<d, d> P1_initial(50,
                                                            /*evaluate=*/
                                                            [phi1_initial](const auto& x, const auto& param) {
                                                              // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto ret = FieldVector<double, d>({1. + rand1, 0. +
                                                              // rand2});
                                                              auto ret = FieldVector<double, d>({1., 0.});
                                                              ret *= (phi1_initial.evaluate(x, param) + 1.) / 2.;
                                                              return ret;
                                                            },
                                                            /*name=*/"P_initial");
      const XT::Functions::GenericFunction<d, d> P2_initial(50,
                                                            /*evaluate=*/
                                                            [phi2_initial](const auto& x, const auto& param) {
                                                              // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto ret = FieldVector<double, d>({1. + rand1, 0. +
                                                              // rand2});
                                                              auto ret = FieldVector<double, d>({1., 0.});
                                                              ret *= (phi2_initial.evaluate(x, param) + 1.) / 2.;
                                                              return ret;
                                                            },
                                                            /*name=*/"P_initial");

      const XT::Functions::ConstantFunction<d, d> u_initial(0.);

      // interpolate initial and boundary values
      default_interpolation(phi1_initial, phi[0]);
      default_interpolation(mu1_initial, mu[0]);
      default_interpolation(phi2_initial, phi[1]);
      default_interpolation(mu2_initial, mu[1]);
      default_interpolation(P1_initial, P[0]);
      default_interpolation(P2_initial, P[1]);
      default_interpolation(u_initial, u);
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    }

    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
    const XT::Functions::ConstantFunction<d> minus_one(-1.);
    for (size_t kk = 0; kk < num_cells; ++kk)
      boundary_interpolation(minus_one, phi[kk], all_dirichlet_boundary_info, XT::Grid::DirichletBoundary{});


    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    write_files(true, u, p, P, Pnat, phi, phinat, mu, filename, 0, t, subsampling);

    OfieldSolver ofield_solver(u, P, Pnat, phi, xi, kappa, c_1, Pa, beta, linearize);
    PfieldSolver pfield_solver(
        u, P, phi, phinat, mu, gamma, c_1, Pa, Be, Ca, beta, epsilon, In, all_dirichlet_boundary_info, linearize);
    StokesSolver stokes_solver(u,
                               p,
                               P,
                               Pnat,
                               phi,
                               phinat,
                               Re,
                               Fa,
                               xi,
                               vol_domain,
                               all_dirichlet_boundary_info,
                               StokesSolverType::eigen_sparse_lu);

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t << std::endl;
      for (size_t kk = 0; kk < num_cells; ++kk) {
        pfield_solver.apply(actual_dt, kk);
        std::cout << "Pfield " << kk << " done" << std::endl;
        ofield_solver.apply(actual_dt, kk);
        std::cout << "Ofield " << kk << " done" << std::endl;
      }
      stokes_solver.apply();
      std::cout << "Stokes done" << std::endl;

      t += actual_dt;

      // check if data should be written in this timestep (and write)
      if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t, next_save_time)) {
        write_files(true, u, p, P, Pnat, phi, phinat, mu, filename, save_step_counter, t, subsampling);
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
