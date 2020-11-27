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

#include <dune/istl/paamg/fastamg.hh>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-view.hh>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/UmfPackSupport>

#include "linear-solver-types.hh"

namespace Eigen {


class PfieldIncompleteLUTPreconditioner
{
  typedef double Scalar;
  typedef Matrix<Scalar, Dynamic, 1> Vector;
  typedef IncompleteLUT<double> IncompleteLUTSolverType;

public:
  typedef typename Vector::StorageIndex StorageIndex;
  enum
  {
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic
  };

  PfieldIncompleteLUTPreconditioner()
    : incomplete_lut_solver_(nullptr)
  {}

  explicit PfieldIncompleteLUTPreconditioner(const IncompleteLUTSolverType* ilut_solver)
    : incomplete_lut_solver_(ilut_solver)
  {}

  Index rows() const
  {
    return incomplete_lut_solver_->rows();
  }
  Index cols() const
  {
    return incomplete_lut_solver_->cols();
  }

  template <typename MatType>
  PfieldIncompleteLUTPreconditioner& analyzePattern(const MatType& /*mat*/)
  {
    return *this;
  }

  template <typename MatType>
  PfieldIncompleteLUTPreconditioner& factorize(const MatType& /*mat*/)
  {
    return *this;
  }

  template <typename MatType>
  PfieldIncompleteLUTPreconditioner& compute(const MatType& mat)
  {
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs& b, Dest& x) const
  {
    x = incomplete_lut_solver_->solve(b);
  }

  void set_lut_solver_(const IncompleteLUTSolverType* ilut_solver)
  {
    incomplete_lut_solver_ = ilut_solver;
  }

  template <typename Rhs>
  inline const Solve<PfieldIncompleteLUTPreconditioner, Rhs> solve(const MatrixBase<Rhs>& b) const
  {
    return Solve<PfieldIncompleteLUTPreconditioner, Rhs>(*this, b.derived());
  }

  ComputationInfo info()
  {
    return Success;
  }

protected:
  const IncompleteLUTSolverType* incomplete_lut_solver_;
};

} // namespace Eigen

namespace Dune {


// forward declarations
template <class X, class Y, class F>
class FGMResSolver;

template <class X, class Y, class F>
class GMResSolver;

template <class X>
class BiCGSTABSolver;

template <class X, class Y>
class LinearOperator;

template <class M, class X>
class LinearOperatorWrapper;

template <class X, class M, bool use_incomplete_lut>
class Matrix2InverseOperator;

template <class X>
class ScalarProduct;

template <class X, class Y>
class IterativeSolver;

template <class X, class Y>
class Preconditioner;

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

template <class X, class M>
class OfieldScalarProduct;


class CellModelLinearSolverWrapper
{
public:
  using R = double;
  using MatrixType = XT::LA::EigenRowMajorSparseMatrix<R>;
  using VectorType = XT::LA::CommonDenseVector<R>;
  using EigenVectorType = XT::LA::EigenDenseVector<R>;
  using ColMajorBackendType = ::Eigen::SparseMatrix<R, ::Eigen::ColMajor>;
  using RowMajorBackendType = typename MatrixType::BackendType;
  // using LUSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using LUSolverType = ::Eigen::UmfPackLU<ColMajorBackendType>;
  using LDLTSolverType = ::Eigen::SimplicialLDLT<ColMajorBackendType>;
  using CGSolverType = ::Eigen::ConjugateGradient<RowMajorBackendType, Eigen::Lower | Eigen::Upper>;
  using CGIncompleteCholeskySolverType =
      ::Eigen::ConjugateGradient<RowMajorBackendType,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IncompleteCholesky<R, Eigen::Lower | Eigen::Upper>>;
  using FGMResSolverType = Dune::FGMResSolver<EigenVectorType, EigenVectorType, EigenVectorType>;
  using GMResSolverType = Dune::GMResSolver<EigenVectorType, EigenVectorType, EigenVectorType>;
  using BiCGSTABSolverType = Dune::BiCGSTABSolver<EigenVectorType>;
  using DuneLinearOperatorType = Dune::LinearOperator<EigenVectorType, EigenVectorType>;
  using LinearOperatorType = LinearOperatorWrapper<MatrixType, EigenVectorType>;
  using ScalarProductType = Dune::ScalarProduct<EigenVectorType>;
  using IterativeSolverType = Dune::IterativeSolver<EigenVectorType, EigenVectorType>;
  using InverseOperatorType = Dune::InverseOperator<EigenVectorType, EigenVectorType>;
  using SystemMatrixLinearOperatorType = MatrixToLinearOperator<EigenVectorType, MatrixType>;
  using PreconditionerType = Dune::Preconditioner<EigenVectorType, EigenVectorType>;
  // using AMGPreconditionerType = Dune::Amg::FastAMG<SystemMatrixLinearOperatorType,EigenVectorType>;
  using IdentityPreconditionerType = IdentityPreconditioner<EigenVectorType>;
  using IterativeSolverPreconditionerType = IterativeSolverPreconditioner<EigenVectorType>;
  // using Matrix2InverseOperatorType = Matrix2InverseOperator<EigenVectorType, MatrixType, false>;
  using Matrix2InverseOperatorType = Matrix2InverseOperator<EigenVectorType, MatrixType, true>;

  CellModelLinearSolverWrapper(std::shared_ptr<LinearOperatorType>,
                               std::shared_ptr<ScalarProductType>,
                               const MatrixType& M,
                               MatrixType& S_,
                               std::shared_ptr<MatrixType>& S_preconditioner,
                               const CellModelLinearSolverType solver_type,
                               const CellModelMassMatrixSolverType mass_matrix_solver_type,
                               const size_t num_cells,
                               const double outer_reduction = 1e-13,
                               const int outer_restart = 100,
                               const int outer_verbose = 0,
                               const double inner_reduction = 1e-3,
                               const int inner_maxit = 10,
                               const int inner_verbose = 0);

  static std::shared_ptr<MatrixType>
  create_system_matrix(const bool is_schur_solver, const size_t size, const XT::LA::SparsityPatternDefault& pattern);

  static std::shared_ptr<MatrixType> create_preconditioner_matrix(const CellModelLinearSolverType solver_type,
                                                                  const size_t size);

  // Has to be called after mass matrix is assembled.
  void setup();

  void prepare(const size_t cell, const bool restricted);

  // Calling this method will result in ret = M^{-1} rhs.
  // Note: Writing rhs_mu.backend() = solver->solve(rhs_mu.backend()) gives wrong results for CG solver, ret may not be
  // the same vector as rhs!
  void apply_inverse_mass_matrix(const EigenVectorType& rhs,
                                 EigenVectorType& ret,
                                 const EigenVectorType* initial_guess = nullptr) const;

  VectorType apply_system_matrix_solver(const VectorType& rhs, const size_t /*cell*/) const;

  void apply_outer_solver(EigenVectorType& ret, EigenVectorType& rhs) const;

  std::shared_ptr<InverseOperatorType>
  create_preconditioner_inverse_op(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                   ScalarProductType& scalar_product,
                                   std::shared_ptr<MatrixType>& S_preconditioner,
                                   const R inner_reduction,
                                   const int inner_verbose,
                                   const int inner_maxit,
                                   const size_t vector_size);

  std::shared_ptr<PreconditionerType> create_preconditioner(std::shared_ptr<DuneLinearOperatorType> linear_op);

  std::shared_ptr<IterativeSolverType> create_iterative_solver(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                                               ScalarProductType& scalar_product,
                                                               const R outer_reduction,
                                                               const size_t outer_restart,
                                                               const int outer_verbose,
                                                               const size_t vector_size);


  std::shared_ptr<LinearOperatorType> linear_operator_;
  std::shared_ptr<ScalarProductType> scalar_product_;
  const MatrixType& M_;
  MatrixType& S_;
  std::shared_ptr<MatrixType>& S_preconditioner_;
  const CellModelLinearSolverType solver_type_;
  const CellModelMassMatrixSolverType mass_matrix_solver_type_;
  const bool is_schur_solver_;
  std::shared_ptr<ColMajorBackendType> S_colmajor_;
  std::shared_ptr<LUSolverType> mass_matrix_lu_solver_;
  std::shared_ptr<LDLTSolverType> mass_matrix_ldlt_solver_;
  std::shared_ptr<CGSolverType> mass_matrix_cg_solver_;
  std::shared_ptr<CGIncompleteCholeskySolverType> mass_matrix_cg_incomplete_cholesky_solver_;
  std::shared_ptr<LUSolverType> direct_solver_;
  std::shared_ptr<IdentityPreconditionerType> identity_preconditioner_;
  std::shared_ptr<InverseOperatorType> preconditioner_inverse_op_;
  std::shared_ptr<PreconditionerType> preconditioner_;
  std::shared_ptr<IterativeSolverType> outer_solver_;
  mutable std::vector<EigenVectorType> previous_update_;
  // ::Eigen::IncompleteLUT<double> eigen_lut_solver_;
  // ::Eigen::PfieldIncompleteLUTPreconditioner pfield_lut_preconditioner_;
  // mutable ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::PfieldIncompleteLUTPreconditioner> bicgsolver_;
  mutable Dune::InverseOperatorResult statistics_;
}; // class CellModelLinearSolverWrapper<...>


class PfieldLinearSolver
{
  using ThisType = PfieldLinearSolver;

public:
  using R = typename CellModelLinearSolverWrapper::R;
  using MatrixType = typename CellModelLinearSolverWrapper::MatrixType;
  using VectorType = typename CellModelLinearSolverWrapper::VectorType;
  using EigenVectorType = typename CellModelLinearSolverWrapper::EigenVectorType;
  using LinearOperatorType = typename CellModelLinearSolverWrapper::LinearOperatorType;
  using ScalarProductType = typename CellModelLinearSolverWrapper::ScalarProductType;
  using SystemMatrixLinearOperatorType = MatrixToLinearOperator<EigenVectorType, MatrixType>;
  using PhiScalarProductType = MassMatrixScalarProduct<EigenVectorType, MatrixType>;
  using PfieldScalarProductType = PfieldScalarProduct<EigenVectorType, MatrixType>;
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;

  PfieldLinearSolver(const double dt,
                     const double gamma,
                     const double epsilon,
                     const double Be,
                     const double Ca,
                     const MatrixType& M,
                     const MatrixType& E,
                     const MatrixType& B,
                     const MatrixType& Dphi_f_incl_coeffs_and_sign,
                     const MatrixType& Dmu_f,
                     const CellModelLinearSolverType solver_type,
                     const CellModelMassMatrixSolverType mass_matrix_solver_type,
                     const XT::LA::SparsityPatternDefault& submatrix_pattern,
                     const size_t num_cells,
                     const double outer_reduction = 1e-10,
                     const int outer_restart = 100,
                     const int outer_verbose = 0,
                     const double inner_reduction = 1e-3,
                     const int inner_maxit = 10,
                     const int inner_verbose = 0);

  // Has to be called after mass matrix is assembled.
  void setup();

  void set_params(const double gamma, const double epsilon, const double Be, const double Ca);

  void prepare(const size_t cell, const bool restricted = false);

  VectorType apply(const VectorType& rhs, const size_t cell) const;

  bool is_schur_solver() const;

  std::shared_ptr<MatrixType> create_pfield_preconditioner_matrix(const CellModelLinearSolverType solver_type,
                                                                  const XT::LA::SparsityPatternDefault& pattern);

  // creates sparsity pattern of phasefield system matrix
  static XT::LA::SparsityPatternDefault system_matrix_pattern(const XT::LA::SparsityPatternDefault& submatrix_pattern);

  // creates sparsity pattern of phasefield preconditioner matrix
  static XT::LA::SparsityPatternDefault
  preconditioner_matrix_pattern(const XT::LA::SparsityPatternDefault& submatrix_pattern);

  // Calling this method will result in ret = M^{-1} rhs.
  // Note: Writing rhs_mu.backend() = solver->solve(rhs_mu.backend()) gives wrong results for CG solver, ret may not be
  // the same vector as rhs!
  void apply_inverse_mass_matrix(const EigenVectorType& rhs,
                                 EigenVectorType& ret,
                                 const EigenVectorType* initial_guess = nullptr) const;

  void set_nonlinear_part_of_S() const;

  VectorType apply_schur_solver(const VectorType& rhs, const size_t cell) const;

  std::shared_ptr<LinearOperatorType> create_linear_operator();

  std::shared_ptr<Dune::ScalarProduct<EigenVectorType>> create_scalar_product();

  R dt_;
  R gamma_;
  R epsilon_;
  R Be_;
  R Ca_;
  const MatrixType& M_;
  const MatrixType& E_;
  const MatrixType& B_;
  const MatrixType& Dphi_f_incl_coeffs_and_sign_;
  const MatrixType& Dmu_f_;
  CellModelLinearSolverType solver_type_;
  const bool is_schur_solver_;
  const XT::LA::SparsityPatternDefault& submatrix_pattern_;
  const size_t size_phi_;
  std::shared_ptr<MatrixType> S_;
  std::shared_ptr<MatrixType> S_preconditioner_;
  CellModelLinearSolverWrapper wrapper_;
  mutable MatrixViewType S_00_;
  mutable MatrixViewType S_01_;
  mutable MatrixViewType S_10_;
  mutable MatrixViewType S_11_;
  mutable MatrixViewType S_12_;
  mutable MatrixViewType S_20_;
  mutable MatrixViewType S_22_;
  mutable EigenVectorType phi_tmp_eigen_;
  mutable EigenVectorType phi_tmp_eigen2_;
  mutable EigenVectorType phi_tmp_eigen3_;
  mutable EigenVectorType phi_tmp_eigen4_;
  mutable std::vector<EigenVectorType> previous_phi_update_;
}; // class PfieldLinearSolver<...>


class OfieldLinearSolver
{
  using ThisType = OfieldLinearSolver;

public:
  using R = typename CellModelLinearSolverWrapper::R;
  using MatrixType = typename CellModelLinearSolverWrapper::MatrixType;
  using VectorType = typename CellModelLinearSolverWrapper::VectorType;
  using EigenVectorType = typename CellModelLinearSolverWrapper::EigenVectorType;
  using LinearOperatorType = typename CellModelLinearSolverWrapper::LinearOperatorType;
  using ScalarProductType = typename CellModelLinearSolverWrapper::ScalarProductType;
  using SchurMatrixLinearOperatorType = MatrixToLinearOperator<EigenVectorType, MatrixType>;
  using SystemMatrixLinearOperatorType = SchurMatrixLinearOperatorType;
  using PScalarProductType = MassMatrixScalarProduct<EigenVectorType, MatrixType>;
  using OfieldScalarProductType = OfieldScalarProduct<EigenVectorType, MatrixType>;
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;

  OfieldLinearSolver(const double dt,
                     const double kappa,
                     const double Pa,
                     const MatrixType& M,
                     const MatrixType& E_,
                     const MatrixType& B,
                     const MatrixType& C_incl_coeffs_and_sign,
                     const MatrixType& Dd_f_incl_coeffs_and_sign,
                     const CellModelLinearSolverType solver_type,
                     const CellModelMassMatrixSolverType mass_matrix_solver_type,
                     const XT::LA::SparsityPatternDefault& submatrix_pattern,
                     const size_t num_cells,
                     const double outer_reduction = 1e-10,
                     const int outer_restart = 100,
                     const int outer_verbose = 0,
                     const double inner_reduction = 1e-3,
                     const int inner_maxit = 10,
                     const int inner_verbose = 0);

  // Has to be called after mass matrix is assembled.
  void setup();

  void set_params(const double dt, const double kappa, const double Pa);

  void prepare(const size_t cell, const bool restricted = false);

  VectorType apply(const VectorType& rhs, const size_t cell) const;

  bool is_schur_solver() const;

  // creates sparsity pattern of phasefield system matrix
  static XT::LA::SparsityPatternDefault system_matrix_pattern(const XT::LA::SparsityPatternDefault& submatrix_pattern);

  // Calling this method will result in ret = M^{-1} rhs.
  // Note: Writing rhs_mu.backend() = solver->solve(rhs_mu.backend()) gives wrong results for CG solver, ret may not be
  // the same vector as rhs!
  void apply_inverse_mass_matrix(const EigenVectorType& rhs,
                                 EigenVectorType& ret,
                                 const EigenVectorType* initial_guess = nullptr) const;

  MatrixType& schur_matrix();

  void set_nonlinear_part_of_S() const;

  VectorType apply_schur_solver(const VectorType& rhs, const size_t cell) const;

  std::shared_ptr<LinearOperatorType> create_linear_operator();

  std::shared_ptr<Dune::ScalarProduct<EigenVectorType>> create_scalar_product();

  R dt_;
  R kappa_;
  R Pa_;
  const MatrixType& M_;
  const MatrixType& E_;
  const MatrixType& B_;
  const MatrixType& C_incl_coeffs_and_sign_;
  const MatrixType& Dd_f_incl_coeffs_and_sign_;
  const bool is_schur_solver_;
  const size_t size_P_;
  std::shared_ptr<MatrixType> S_;
  std::shared_ptr<MatrixType> S_preconditioner_;
  std::shared_ptr<MatrixType> S_schur_;
  CellModelLinearSolverWrapper wrapper_;
  mutable MatrixViewType S_00_;
  mutable MatrixViewType S_01_;
  mutable MatrixViewType S_10_;
  mutable MatrixViewType S_11_;
  mutable EigenVectorType P_tmp_eigen_;
  mutable EigenVectorType P_tmp_eigen2_;
  mutable std::vector<EigenVectorType> previous_P_update_;
}; // class OfieldLinearSolver<...>


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_LINEARSOLVERS_HH
