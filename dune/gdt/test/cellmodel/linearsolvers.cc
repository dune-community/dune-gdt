// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include "config.h"

#include <memory>
#include <set>
#include <vector>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver/istl/preconditioners.hh>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "fgmres.hh"
#include "preconditioners.hh"
#include "linearoperators.hh"
#include "scalarproducts.hh"
#include "linearsolvers.hh"

namespace Dune {


PfieldLinearSolver::PfieldLinearSolver(const double gamma,
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
                                       const double outer_reduction,
                                       const int outer_restart,
                                       const int outer_verbose,
                                       const double inner_reduction,
                                       const int inner_maxit,
                                       const int inner_verbose)
  : gamma_(gamma)
  , epsilon_(epsilon)
  , Be_(Be)
  , Ca_(Ca)
  , M_(M)
  , M_ell_(M_ell)
  , D_(D)
  , G_(G)
  , M_nonlin_(M_nonlin)
  , A_boundary_(A_boundary)
  , solver_type_(solver_type)
  , mass_matrix_solver_type_(mass_matrix_solver_type)
  , size_phi_(M_.rows())
  , phi_dirichlet_dofs_(phi_dirichlet_dofs)
  , is_schur_solver_(solver_type_ == PfieldLinearSolverType::schur_gmres
                     || solver_type_ == PfieldLinearSolverType::schur_fgmres_gmres
                     || solver_type_ == PfieldLinearSolverType::schur_fgmres_bicgstab)
  , S_(is_schur_solver_ ? std::make_shared<MatrixType>(0, 0)
                        : std::make_shared<MatrixType>(
                              3 * size_phi_, 3 * size_phi_, system_matrix_pattern(size_phi_, submatrix_pattern), 0))
  , S_colmajor_(solver_type_ == PfieldLinearSolverType::direct ? std::make_shared<ColMajorBackendType>(S_->backend())
                                                               : nullptr)
  , S_00_(*S_, 0, size_phi_, 0, size_phi_)
  , S_01_(*S_, 0, size_phi_, size_phi_, 2 * size_phi_)
  , S_10_(*S_, size_phi_, 2 * size_phi_, 0, size_phi_)
  , S_11_(*S_, size_phi_, 2 * size_phi_, size_phi_, 2 * size_phi_)
  , S_12_(*S_, size_phi_, 2 * size_phi_, 2 * size_phi_, 3 * size_phi_)
  , S_20_(*S_, 2 * size_phi_, 3 * size_phi_, 0, size_phi_)
  , S_22_(*S_, 2 * size_phi_, 3 * size_phi_, 2 * size_phi_, 3 * size_phi_)
  , mass_matrix_lu_solver_(
        mass_matrix_solver_type_ == PfieldMassMatrixSolverType::sparse_lu ? std::make_shared<LUSolverType>() : nullptr)
  , mass_matrix_cg_solver_(mass_matrix_solver_type_ == PfieldMassMatrixSolverType::cg ? std::make_shared<CGSolverType>()
                                                                                      : nullptr)
  , mass_matrix_cg_incomplete_cholesky_solver_(mass_matrix_solver_type_
                                                       == PfieldMassMatrixSolverType::cg_incomplete_cholesky
                                                   ? std::make_shared<CGIncompleteCholeskySolverType>()
                                                   : nullptr)
  , direct_solver_(solver_type_ == PfieldLinearSolverType::direct ? std::make_shared<LUSolverType>() : nullptr)
  , linear_operator_(create_linear_operator())
  , scalar_product_(create_scalar_product())
  , identity_preconditioner_(std::make_shared<IdentityPreconditionerType>(SolverCategory::Category::sequential))
  , preconditioner_solver_(
        create_preconditioner_solver(linear_operator_, *scalar_product_, inner_reduction, inner_verbose, inner_maxit))
  , preconditioner_(std::make_shared<IterativeSolverPreconditionerType>(preconditioner_solver_,
                                                                        SolverCategory::Category::sequential))
  , outer_solver_(
        create_iterative_solver(linear_operator_, *scalar_product_, outer_reduction, outer_restart, outer_verbose))
  , phi_tmp_eigen_(size_phi_, 0., 0)
  , phi_tmp_eigen2_(size_phi_, 0., 0)
  , phi_tmp_eigen3_(size_phi_, 0., 0)
  , phi_tmp_eigen4_(size_phi_, 0., 0)
  , last_pfield_update_(3 * size_phi_, 0., 0)
  , last_phi_update_(num_cells, EigenVectorType(size_phi_, 0., 0))
{
  if (direct_solver_)
    direct_solver_->analyzePattern(*S_colmajor_);
  if (!is_schur_solver_) {
    S_11_ = M_;
    S_22_ = M_;
  }
}

// Has to be called after mass matrix is assembled
void PfieldLinearSolver::setup()
{
  if (mass_matrix_lu_solver_)
    mass_matrix_lu_solver_->compute(M_.backend());
  if (mass_matrix_cg_solver_)
    mass_matrix_cg_solver_->compute(M_.backend());
  if (mass_matrix_cg_incomplete_cholesky_solver_)
    mass_matrix_cg_incomplete_cholesky_solver_->compute(M_.backend());
}

void PfieldLinearSolver::set_params(const XT::Common::Parameter& param)
{
  gamma_ = param.get("gamma")[0];
  epsilon_ = param.get("epsilon")[0];
  Be_ = param.get("Be")[0];
  Ca_ = param.get("Ca")[0];
  linear_operator_->set_params(param);
}

void PfieldLinearSolver::prepare(const double dt, const size_t cell, const bool restricted)
{
  dt_ = dt;
  linear_operator_->prepare(dt_, cell, restricted);
}

PfieldLinearSolver::VectorType PfieldLinearSolver::apply(const VectorType& rhs, const size_t cell) const
{
  if (is_schur_solver_)
    return apply_schur_solver(rhs, cell);
  else
    return apply_system_matrix_solver(rhs, cell);
}

// creates sparsity pattern of phasefield system matrix
XT::LA::SparsityPatternDefault
PfieldLinearSolver::system_matrix_pattern(const size_t n, const XT::LA::SparsityPatternDefault& submatrix_pattern)
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

// Calling this method will result in ret = M^{-1} rhs.
// Note: Writing rhs_mu.backend() = solver->solve(rhs_mu.backend()) gives wrong results for CG solver, ret may not be
// the same vector as rhs!
void PfieldLinearSolver::apply_inverse_mass_matrix(const EigenVectorType& rhs,
                                                   EigenVectorType& ret,
                                                   const EigenVectorType* initial_guess) const
{
  if (mass_matrix_solver_type_ == PfieldMassMatrixSolverType::sparse_lu)
    ret.backend() = mass_matrix_lu_solver_->solve(rhs.backend());
  else if (mass_matrix_solver_type_ == PfieldMassMatrixSolverType::cg) {
    if (initial_guess)
      ret.backend() = mass_matrix_cg_solver_->solveWithGuess(rhs.backend(), initial_guess->backend());
    else
      ret.backend() = mass_matrix_cg_solver_->solve(rhs.backend());
  } else if (mass_matrix_solver_type_ == PfieldMassMatrixSolverType::cg_incomplete_cholesky) {
    if (initial_guess)
      ret.backend() =
          mass_matrix_cg_incomplete_cholesky_solver_->solveWithGuess(rhs.backend(), initial_guess->backend());
    else
      mass_matrix_cg_incomplete_cholesky_solver_->solve(rhs.backend());
  } else {
    DUNE_THROW(Dune::NotImplemented, "Unknown mass matrix solver type!");
  }
}

const std::set<size_t>& PfieldLinearSolver::dirichlet_dofs() const
{
  return phi_dirichlet_dofs_;
}

void PfieldLinearSolver::fill_S() const
{
  S_00_ = M_;
  S_00_.axpy(dt_, D_);
  S_01_ = M_ell_;
  S_01_ *= gamma_ * dt_;
  S_10_ = G_;
  S_12_ = M_ell_;
  S_12_ *= 1. / Be_;
  S_12_.axpy(1. / (Be_ * std::pow(epsilon_, 2)), M_nonlin_);
  S_12_.axpy(1. / Ca_, M_);
  S_20_ = M_ell_;
  S_20_ *= epsilon_;
  S_20_.axpy(1. / epsilon_, M_nonlin_);
  S_20_.axpy(-epsilon_, A_boundary_);
  for (const auto& DoF : phi_dirichlet_dofs_) {
    S_00_.unit_row(DoF);
    S_01_.clear_row(DoF);
  }
  if (!is_schur_solver_) {
    *S_colmajor_ = S_->backend();
    S_colmajor_->makeCompressed();
  }
}

PfieldLinearSolver::VectorType PfieldLinearSolver::apply_system_matrix_solver(const VectorType& rhs,
                                                                              const size_t /*cell*/) const
{
  fill_S();
  auto rhs_eig = XT::Common::convert_to<EigenVectorType>(rhs);
  EigenVectorType update(rhs.size());
  Dune::InverseOperatorResult res;
  switch (solver_type_) {
    case PfieldLinearSolverType::direct:
      direct_solver_->factorize(*S_colmajor_);
      update.backend() = direct_solver_->solve(rhs_eig.backend());
      return XT::Common::convert_to<VectorType>(update);
    case PfieldLinearSolverType::gmres:
    case PfieldLinearSolverType::fgmres_gmres:
    case PfieldLinearSolverType::fgmres_bicgstab:
      outer_solver_->apply(last_pfield_update_, rhs_eig, res);
      update = last_pfield_update_;
    case PfieldLinearSolverType::schur_gmres:
    case PfieldLinearSolverType::schur_fgmres_bicgstab:
    case PfieldLinearSolverType::schur_fgmres_gmres:
      DUNE_THROW(Dune::NotImplemented, "Unknown solver type for system matrix!");
  }
  return XT::Common::convert_to<VectorType>(update);
}

PfieldLinearSolver::VectorType PfieldLinearSolver::apply_schur_solver(const VectorType& rhs, const size_t cell) const
{
  VectorType ret(3 * size_phi_, 0., 0);
  VectorViewType ret_phi(ret, 0, size_phi_);
  VectorViewType ret_phinat(ret, size_phi_, 2 * size_phi_);
  VectorViewType ret_mu(ret, 2 * size_phi_, 3 * size_phi_);
  ConstVectorViewType r0(rhs, 0, size_phi_);
  ConstVectorViewType r1(rhs, size_phi_, 2 * size_phi_);
  ConstVectorViewType r2(rhs, 2 * size_phi_, 3 * size_phi_);
  // apply inverse of first factor matrix to rhs
  // compute rhs_mu = C^{-1} r2
  auto& rhs_mu = phi_tmp_eigen4_;
  rhs_mu = r2;
  apply_inverse_mass_matrix(rhs_mu, phi_tmp_eigen3_);
  rhs_mu = phi_tmp_eigen3_;
  // compute rhs_phinat = H^{-1}(r1 - J C^{-1} r2) = H^{-1}(r1 - J rhs_mu)
  auto& rhs_phinat = phi_tmp_eigen3_;
  rhs_phinat = r1;
  M_.mv(rhs_mu, phi_tmp_eigen_);
  rhs_phinat.axpy(-1. / Ca_, phi_tmp_eigen_);
  M_ell_.mv(rhs_mu, phi_tmp_eigen_);
  rhs_phinat.axpy(-1. / Be_, phi_tmp_eigen_);
  M_nonlin_.mv(rhs_mu, phi_tmp_eigen_);
  rhs_phinat.axpy(-1. / (Be_ * std::pow(epsilon_, 2)), phi_tmp_eigen_);
  apply_inverse_mass_matrix(rhs_phinat, phi_tmp_eigen_);
  rhs_phinat = phi_tmp_eigen_;
  // compute rhs_phi = r0 - dt E H^{-1} (r1 - J C^{-1} r2) = r0 - dt E rhs_phinat
  auto& rhs_phi = phi_tmp_eigen2_;
  M_ell_.mv(rhs_phinat, rhs_phi);
  rhs_phi *= -gamma_ * dt_;
  // E has empty rows for the dirichlet dofs
  for (const auto& DoF : phi_dirichlet_dofs_)
    rhs_phi[DoF] = 0.;
  rhs_phi += r0;
  // now solve for phi
  Dune::InverseOperatorResult res;
  outer_solver_->apply(last_phi_update_[cell], rhs_phi, res);
  const auto& phi = last_phi_update_[cell];
  ret_phi = phi;
  // store C^{-1} A phi in phi_tmp_eigen2_
  M_ell_.mv(phi, phi_tmp_eigen2_);
  phi_tmp_eigen2_ *= epsilon_;
  M_nonlin_.mv(phi, phi_tmp_eigen_);
  phi_tmp_eigen2_.axpy(1. / epsilon_, phi_tmp_eigen_);
  A_boundary_.mv(phi, phi_tmp_eigen_);
  phi_tmp_eigen2_.axpy(-epsilon_, phi_tmp_eigen_);
  apply_inverse_mass_matrix(phi_tmp_eigen2_, phi_tmp_eigen_);
  phi_tmp_eigen2_ = phi_tmp_eigen_;
  // calculate mu
  ret_mu = rhs_mu;
  ret_mu -= phi_tmp_eigen2_;
  // now calculate phinat
  // calculate H^{-1} (G - J C^{-1} A) phi
  G_.mv(phi, phi_tmp_eigen4_);
  M_.mv(phi_tmp_eigen2_, phi_tmp_eigen_);
  phi_tmp_eigen4_.axpy(-1. / Ca_, phi_tmp_eigen_);
  M_ell_.mv(phi_tmp_eigen2_, phi_tmp_eigen_);
  phi_tmp_eigen4_.axpy(-1. / Be_, phi_tmp_eigen_);
  M_nonlin_.mv(phi_tmp_eigen2_, phi_tmp_eigen_);
  phi_tmp_eigen4_.axpy(-1. / (Be_ * std::pow(epsilon_, 2)), phi_tmp_eigen_);
  apply_inverse_mass_matrix(phi_tmp_eigen4_, phi_tmp_eigen_);
  phi_tmp_eigen4_ = phi_tmp_eigen_;
  ret_phinat = rhs_phinat;
  ret_phinat -= phi_tmp_eigen4_;
  // return
  return ret;
}

std::shared_ptr<PfieldLinearSolver::LinearOperatorType> PfieldLinearSolver::create_linear_operator()
{
  if (is_schur_solver_) {
    return std::make_shared<SchurMatrixLinearOperatorType>(
        M_, D_, M_ell_, G_, M_nonlin_, A_boundary_, *this, gamma_, epsilon_, Be_, Ca_);
  } else {
    return std::make_shared<SystemMatrixLinearOperatorType>(*S_);
  }
}

std::shared_ptr<Dune::ScalarProduct<PfieldLinearSolver::EigenVectorType>> PfieldLinearSolver::create_scalar_product()
{
  if (is_schur_solver_)
    return std::make_shared<PhiScalarProductType>(M_);
  else
    return std::make_shared<PfieldScalarProductType>(M_);
}

std::shared_ptr<PfieldLinearSolver::IterativeSolverType>
PfieldLinearSolver::create_preconditioner_solver(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                                 ScalarProductType& scalar_product,
                                                 const R inner_reduction,
                                                 const int inner_verbose,
                                                 const int inner_maxit)
{
  switch (solver_type_) {
    case PfieldLinearSolverType::fgmres_gmres:
    case PfieldLinearSolverType::schur_fgmres_gmres:
      return std::make_shared<RestartedGMResSolverType>(*linear_op,
                                                        scalar_product,
                                                        *identity_preconditioner_,
                                                        inner_reduction,
                                                        10000, // do not restart
                                                        inner_maxit,
                                                        inner_verbose);
    case PfieldLinearSolverType::fgmres_bicgstab:
    case PfieldLinearSolverType::schur_fgmres_bicgstab:
      return std::make_shared<BiCGSTABSolverType>(
          *linear_op, scalar_product, *identity_preconditioner_, inner_reduction, inner_maxit, inner_verbose);
    case PfieldLinearSolverType::direct:
    case PfieldLinearSolverType::gmres:
    case PfieldLinearSolverType::schur_gmres:
      return nullptr;
  }
}

std::shared_ptr<PfieldLinearSolver::IterativeSolverType>
PfieldLinearSolver::create_iterative_solver(std::shared_ptr<DuneLinearOperatorType> linear_op,
                                            ScalarProductType& scalar_product,
                                            const R outer_reduction,
                                            const size_t outer_restart,
                                            const int outer_verbose)
{
  switch (solver_type_) {
    case PfieldLinearSolverType::gmres:
    case PfieldLinearSolverType::schur_gmres:
      return std::make_shared<RestartedGMResSolverType>(
          *linear_op, scalar_product, *identity_preconditioner_, outer_reduction, outer_restart, 10000, outer_verbose);
    case PfieldLinearSolverType::fgmres_gmres:
    case PfieldLinearSolverType::fgmres_bicgstab:
    case PfieldLinearSolverType::schur_fgmres_gmres:
    case PfieldLinearSolverType::schur_fgmres_bicgstab:
      return std::make_shared<FGMResSolverType>(
          *linear_op, scalar_product, *preconditioner_, outer_reduction, outer_restart, 10000, outer_verbose);
    case PfieldLinearSolverType::direct:
      return nullptr;
  }
}


} // namespace Dune