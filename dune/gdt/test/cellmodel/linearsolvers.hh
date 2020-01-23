// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_CELLMODEL_LINEARSOLVERS_HH
#define DUNE_GDT_TEST_CELLMODEL_LINEARSOLVERS_HH

#include <vector>
#include <memory>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "linear-solver-types.hh"

namespace Dune {


// forward declarations
template <class X, class Y, class F>
class FGMResSolver;

template <class X, class Y, class F>
class RestartedGMResSolver;

template <class X>
class BiCGSTABSolver;

template <class X, class Y>
class LinearOperator;

template <class X>
class LinearOperatorWrapper;

template <class X>
class ScalarProduct;

template <class X, class Y>
class IterativeSolver;

template <class X, class Y>
class Preconditioner;

template <class X, class M, class T>
class PfieldPhiMatrixOperator;

template <class X, class M>
class MatrixToLinearOperator;

template <class X>
class IdentityPreconditioner;

template <class X>
class IterativeSolverPreconditioner;

template <class X, class M>
class MassMatrixScalarProduct;

template <class X, class M>
class PfieldScalarProduct;

class PfieldLinearSolver
{
  using ThisType = PfieldLinearSolver;

public:
  using R = double;
  using MatrixType = XT::LA::EigenRowMajorSparseMatrix<R>;
  using VectorType = XT::LA::CommonDenseVector<R>;
  using EigenVectorType = XT::LA::EigenDenseVector<R>;
  using ColMajorBackendType = ::Eigen::SparseMatrix<R, ::Eigen::ColMajor>;
  using RowMajorBackendType = typename MatrixType::BackendType;
  using LUSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using CGSolverType = ::Eigen::ConjugateGradient<RowMajorBackendType, Eigen::Lower | Eigen::Upper>;
  using CGIncompleteCholeskySolverType =
      ::Eigen::ConjugateGradient<RowMajorBackendType,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IncompleteCholesky<R, Eigen::Lower | Eigen::Upper>>;
  using FGMResSolverType = Dune::FGMResSolver<EigenVectorType, EigenVectorType, EigenVectorType>;
  using RestartedGMResSolverType = Dune::RestartedGMResSolver<EigenVectorType>;
  using BiCGSTABSolverType = Dune::BiCGSTABSolver<EigenVectorType>;
  using DuneLinearOperatorType = Dune::LinearOperator<EigenVectorType, EigenVectorType>;
  using LinearOperatorType = LinearOperatorWrapper<EigenVectorType>;
  using ScalarProductType = Dune::ScalarProduct<EigenVectorType>;
  using IterativeSolverType = Dune::IterativeSolver<EigenVectorType, EigenVectorType>;
  using SchurMatrixLinearOperatorType = PfieldPhiMatrixOperator<EigenVectorType, MatrixType, ThisType>;
  using SystemMatrixLinearOperatorType = MatrixToLinearOperator<EigenVectorType, MatrixType>;
  using PreconditionerType = Dune::Preconditioner<EigenVectorType, EigenVectorType>;
  using IdentityPreconditionerType = IdentityPreconditioner<EigenVectorType>;
  using IterativeSolverPreconditionerType = IterativeSolverPreconditioner<EigenVectorType>;
  using PhiScalarProductType = MassMatrixScalarProduct<EigenVectorType, MatrixType>;
  using PfieldScalarProductType = PfieldScalarProduct<EigenVectorType, MatrixType>;
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

  PfieldLinearSolver(const double gamma,
                     const double epsilon,
                     const double Be,
                     const double Ca,
                     const MatrixType& M,
                     const MatrixType& M_ell,
                     const MatrixType& D,
                     const MatrixType& G,
                     const MatrixType& M_nonlin,
                     const MatrixType& A_boundary,
                     const PfieldLinearSolverType solver_type,
                     const PfieldMassMatrixSolverType mass_matrix_solver_type,
                     const std::set<size_t>& phi_dirichlet_dofs,
                     const XT::LA::SparsityPatternDefault& submatrix_pattern,
                     const size_t num_cells,
                     const double outer_reduction = 1e-10,
                     const int outer_restart = 100,
                     const int outer_verbose = 0,
                     const double inner_reduction = 1e-3,
                     const int inner_maxit = 10,
                     const int inner_verbose = 0);

  // Has to be called after mass matrix is assembled
  void setup();

  void set_params(const XT::Common::Parameter& param);

  void prepare(const double dt, const size_t cell, const bool restricted = false);

  VectorType apply(const VectorType& rhs, const size_t cell) const;

  // creates sparsity pattern of phasefield system matrix
  static XT::LA::SparsityPatternDefault system_matrix_pattern(const size_t n,
                                                              const XT::LA::SparsityPatternDefault& submatrix_pattern);

  // Calling this method will result in ret = M^{-1} rhs.
  // Note: Writing rhs_mu.backend() = solver->solve(rhs_mu.backend()) gives wrong results for CG solver, ret may not be
  // the same vector as rhs!
  void apply_inverse_mass_matrix(const EigenVectorType& rhs,
                                 EigenVectorType& ret,
                                 const EigenVectorType* initial_guess = nullptr) const;

  const std::set<size_t>& dirichlet_dofs() const;

private:
  void fill_S() const;

  VectorType apply_system_matrix_solver(const VectorType& rhs, const size_t /*cell*/) const;

  VectorType apply_schur_solver(const VectorType& rhs, const size_t cell) const;

  std::shared_ptr<LinearOperatorType> create_linear_operator();

  std::shared_ptr<Dune::ScalarProduct<EigenVectorType>> create_scalar_product();

  std::shared_ptr<IterativeSolverType> create_preconditioner_solver(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                                                    ScalarProductType& scalar_product,
                                                                    const R inner_reduction,
                                                                    const int inner_verbose,
                                                                    const int inner_maxit);

  std::shared_ptr<IterativeSolverType> create_iterative_solver(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                                               ScalarProductType& scalar_product,
                                                               const R outer_reduction,
                                                               const size_t outer_restart,
                                                               const int outer_verbose);

  R gamma_;
  R epsilon_;
  R Be_;
  R Ca_;
  R dt_;
  const MatrixType& M_;
  const MatrixType& M_ell_;
  const MatrixType& D_;
  const MatrixType& G_;
  const MatrixType& M_nonlin_;
  const MatrixType& A_boundary_;
  const PfieldLinearSolverType solver_type_;
  const PfieldMassMatrixSolverType mass_matrix_solver_type_;
  const size_t size_phi_;
  const std::set<size_t>& phi_dirichlet_dofs_;
  const bool is_schur_solver_;
  std::shared_ptr<MatrixType> S_;
  std::shared_ptr<ColMajorBackendType> S_colmajor_;
  mutable MatrixViewType S_00_;
  mutable MatrixViewType S_01_;
  mutable MatrixViewType S_10_;
  mutable MatrixViewType S_11_;
  mutable MatrixViewType S_12_;
  mutable MatrixViewType S_20_;
  mutable MatrixViewType S_22_;
  std::shared_ptr<LUSolverType> mass_matrix_lu_solver_;
  std::shared_ptr<CGSolverType> mass_matrix_cg_solver_;
  std::shared_ptr<CGIncompleteCholeskySolverType> mass_matrix_cg_incomplete_cholesky_solver_;
  std::shared_ptr<LUSolverType> direct_solver_;
  std::shared_ptr<LinearOperatorType> linear_operator_;
  std::shared_ptr<ScalarProductType> scalar_product_;
  std::shared_ptr<IdentityPreconditionerType> identity_preconditioner_;
  std::shared_ptr<IterativeSolverType> preconditioner_solver_;
  std::shared_ptr<IterativeSolverPreconditionerType> preconditioner_;
  std::shared_ptr<IterativeSolverType> outer_solver_;
  mutable EigenVectorType phi_tmp_eigen_;
  mutable EigenVectorType phi_tmp_eigen2_;
  mutable EigenVectorType phi_tmp_eigen3_;
  mutable EigenVectorType phi_tmp_eigen4_;
  mutable EigenVectorType last_pfield_update_;
  mutable std::vector<EigenVectorType> last_phi_update_;
}; // class PfieldLinearSolver<...>


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_LINEARSOLVERS_HH
