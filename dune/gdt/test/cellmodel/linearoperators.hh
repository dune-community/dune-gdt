// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_CELLMODEL_LINEAROPERATORS_HH
#define DUNE_GDT_TEST_CELLMODEL_LINEAROPERATORS_HH

#include <dune/istl/operators.hh>
#include <dune/xt/common/parameter.hh>

namespace Dune {


template <class VectorType>
class LinearOperatorWrapper : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;

public:
  using Vector = VectorType;
  using Field = typename VectorType::ScalarType;

  LinearOperatorWrapper(const size_t vector_length)
    : tmp_vec_(vector_length, 0.)
  {}

  using BaseType::apply;

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  virtual void set_params(const XT::Common::Parameter& /*param*/) {}

  virtual void prepare(const double /*dt*/, const size_t /*cell*/, const bool /*restricted*/ = false) {}

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override
  {
    return SolverCategory::Category::sequential;
  }

private:
  mutable VectorType tmp_vec_;
};


template <class VectorType, class MatrixType, class CellModelSolverType>
class OfieldMatrixLinearPartOperator : public LinearOperatorWrapper<VectorType>
{
  using BaseType = LinearOperatorWrapper<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  OfieldMatrixLinearPartOperator(const Matrix& M,
                                 const Matrix& A,
                                 const Matrix& C_lin,
                                 const double kappa,
                                 const CellModelSolverType* cellmodel_solver)
    : BaseType(2 * M.rows())
    , M_(M)
    , A_(A)
    , C_lin_(C_lin)
    , kappa_(kappa)
    , cellmodel_solver_(cellmodel_solver)
    , size_P_(M_.rows())
    , x_P_(size_P_, 0., 0)
    , x_Pnat_(size_P_, 0., 0)
    , y_P_(size_P_, 0., 0)
    , y_Pnat_(size_P_, 0., 0)
    , tmp_vec_(size_P_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    const auto& input_dofs = cellmodel_solver_->ofield_deim_input_dofs_[cell_];
    const auto& Pnat_begin = cellmodel_solver_->Pnat_deim_input_dofs_begin_[cell_];
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    if (!restricted_) {
      for (size_t ii = 0; ii < size_P_; ++ii) {
        x_P_[ii] = x[ii];
        x_Pnat_[ii] = x[size_P_ + ii];
      }
    } else {
      for (size_t ii = 0; ii < Pnat_begin; ++ii)
        x_P_[input_dofs[ii]] = x[input_dofs[ii]];
      for (size_t ii = Pnat_begin; ii < input_dofs.size(); ++ii)
        x_Pnat_[input_dofs[ii] - size_P_] = x[input_dofs[ii]];
    }
    // apply matrices
    const auto mv = cellmodel_solver_->template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_->template axpy_func<Vector>(restricted_);
    const auto add = cellmodel_solver_->template add_func<Vector>(restricted_);
    // M+dtA
    const auto& P_dofs = cellmodel_solver_->P_deim_output_dofs_[cell_];
    mv(M_, x_P_, y_P_, P_dofs);
    mv(A_, x_P_, tmp_vec_, P_dofs);
    axpy(y_P_, dt_, tmp_vec_, P_dofs);
    // dt B and D
    const auto& Pnat_dofs = cellmodel_solver_->Pnat_deim_output_dofs_[cell_];
    mv(M_, x_Pnat_, y_Pnat_, P_dofs);
    if (restricted_)
      mv(M_, x_Pnat_, y_Pnat_, Pnat_dofs);
    axpy(y_P_, dt_ / kappa_, y_Pnat_, P_dofs);
    // linear part of C
    mv(C_lin_, x_P_, tmp_vec_, Pnat_dofs);
    add(y_Pnat_, tmp_vec_, Pnat_dofs);
    // copy to result vector
    if (!restricted_) {
      for (size_t ii = 0; ii < size_P_; ++ii) {
        y[ii] = y_P_[ii];
        y[size_P_ + ii] = y_Pnat_[ii];
      }
    } else {
      for (const auto& dof : P_dofs)
        y[dof] = y_P_[dof];
      for (const auto& dof : Pnat_dofs)
        y[size_P_ + dof] = y_Pnat_[dof];
    }
  }

  void prepare(const double dt, const size_t cell, const bool restricted = false) override final
  {
    dt_ = dt;
    cell_ = cell;
    restricted_ = restricted;
  }

private:
  const Matrix& M_;
  const Matrix& A_;
  const Matrix& C_lin_;
  const double kappa_;
  const CellModelSolverType* cellmodel_solver_;
  double dt_;
  size_t cell_;
  bool restricted_;
  const size_t size_P_;
  // vectors to store intermediate results
  mutable Vector x_P_;
  mutable Vector x_Pnat_;
  mutable Vector y_P_;
  mutable Vector y_Pnat_;
  mutable Vector tmp_vec_;
};

template <class VectorType, class MatrixType, class SolverType>
class PfieldPhiMatrixOperator : public LinearOperatorWrapper<VectorType>
{
  using BaseType = LinearOperatorWrapper<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  PfieldPhiMatrixOperator(const Matrix& M,
                          const Matrix& D,
                          const Matrix& M_ell,
                          const Matrix& G,
                          const Matrix& M_nonlin,
                          const Matrix& A_boundary,
                          const SolverType& solver,
                          const double gamma,
                          const double epsilon,
                          const double Be,
                          const double Ca)
    : BaseType(M.rows())
    , M_(M)
    , D_(D)
    , M_ell_(M_ell)
    , G_(G)
    , M_nonlin_(M_nonlin)
    , A_boundary_(A_boundary)
    , solver_(solver)
    , A_(M_.backend(), 0)
    , J_(M_.backend(), 0)
    , gamma_(gamma)
    , epsilon_(epsilon)
    , Be_(Be)
    , Ca_(Ca)
    , tmp_vec_(M_.rows(), 0., 0)
    , tmp_vec2_(M_.rows(), 0., 0)
    , last_A_result_(M_.rows(), 0., 0)
    , last_E_result_(M_.rows(), 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    // apply matrices
    // apply M + dt D
    M_.mv(x, y);
    D_.mv(x, tmp_vec_);
    y.axpy(dt_, tmp_vec_);
    // compute (G - J C^{-1} A) x
    A_.mv(x, tmp_vec_);
    solver_.apply_inverse_mass_matrix(tmp_vec_, tmp_vec2_, &last_A_result_);
    last_A_result_ = tmp_vec2_;
    J_.mv(tmp_vec2_, tmp_vec_);
    G_.mv(x, tmp_vec2_);
    tmp_vec2_ -= tmp_vec_;
    // apply -dt E H^{-1}
    // tmp_vec_.backend() = M_inv_.solve(tmp_vec2_.backend());
    solver_.apply_inverse_mass_matrix(tmp_vec2_, tmp_vec_, &last_E_result_);
    last_E_result_ = tmp_vec_.backend();
    M_ell_.mv(tmp_vec_, tmp_vec2_);
    y.axpy(-dt_ * gamma_, tmp_vec2_);
    for (const auto& DoF : solver_.dirichlet_dofs())
      y[DoF] = x[DoF];
  }

  void set_params(const XT::Common::Parameter& param) override final
  {
    gamma_ = param.get("gamma")[0];
    epsilon_ = param.get("epsilon")[0];
    Be_ = param.get("Be")[0];
    Ca_ = param.get("Ca")[0];
  }

  void prepare(const double dt, const size_t /*cell*/, const bool /*restricted*/) override final
  {
    dt_ = dt;
    // precompute A and J (M_nonlin_ may have changed)
    // saves three mvs in each apply, can be dropped if memory is an issue
    J_.backend() = M_nonlin_.backend();
    J_ *= 1. / (Be_ * std::pow(epsilon_, 2));
    J_.axpy(1. / Ca_, M_);
    J_.axpy(1. / Be_, M_ell_);
    A_.backend() = M_nonlin_.backend();
    A_ *= 1. / epsilon_;
    A_.axpy(epsilon_, M_ell_);
    A_.axpy(-epsilon_, A_boundary_);
  }

private:
  const Matrix& M_;
  const Matrix& D_;
  const Matrix& M_ell_;
  const Matrix& G_;
  const Matrix& M_nonlin_;
  const Matrix& A_boundary_;
  const SolverType& solver_;
  MatrixType A_;
  MatrixType J_;
  double gamma_;
  double epsilon_;
  double Be_;
  double Ca_;
  double dt_;
  // vectors to store intermediate results
  mutable Vector tmp_vec_;
  mutable Vector tmp_vec2_;
  mutable Vector last_A_result_;
  mutable Vector last_E_result_;
};

template <class VectorType, class MatrixType>
class MatrixToLinearOperator : public LinearOperatorWrapper<VectorType>
{
  using BaseType = LinearOperatorWrapper<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  MatrixToLinearOperator(const Matrix& M)
    : BaseType(M.rows())
    , M_(M)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    M_.mv(x, y);
  }

private:
  const Matrix& M_;
};

template <class VectorType, class MatrixType, class DirichletConstraintsType, class CellModelSolverType>
class PfieldMatrixLinearPartOperator : public LinearOperatorWrapper<VectorType>
{
  using BaseType = LinearOperatorWrapper<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Dirichlet = DirichletConstraintsType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  PfieldMatrixLinearPartOperator(const Matrix& M,
                                 const Matrix& D,
                                 const Matrix& M_ell,
                                 const Matrix& A_boundary,
                                 const Dirichlet& dirichlet,
                                 const double gamma,
                                 const double eps,
                                 const double Be,
                                 const CellModelSolverType* cellmodel_solver)
    : BaseType(3 * M.rows())
    , M_(M)
    , D_(D)
    , M_ell_(M_ell)
    , A_boundary_(A_boundary)
    , dirichlet_(dirichlet)
    , gamma_(gamma)
    , epsilon_(eps)
    , Be_(Be)
    , cellmodel_solver_(cellmodel_solver)
    , size_phi_(M_.rows())
    , x_phi_(size_phi_, 0., 0)
    , x_phinat_(size_phi_, 0., 0)
    , x_mu_(size_phi_, 0., 0)
    , y_phi_(size_phi_, 0., 0)
    , y_phinat_(size_phi_, 0., 0)
    , y_mu_(size_phi_, 0., 0)
    , tmp_vec_(size_phi_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    const auto& input_dofs = cellmodel_solver_->pfield_deim_input_dofs_[cell_];
    const auto& phinat_begin = cellmodel_solver_->phinat_deim_input_dofs_begin_[cell_];
    const auto& mu_begin = cellmodel_solver_->mu_deim_input_dofs_begin_[cell_];
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi_; ++ii) {
        x_phi_[ii] = x[ii];
        x_phinat_[ii] = x[size_phi_ + ii];
        x_mu_[ii] = x[2 * size_phi_ + ii];
      }
    } else {
      for (size_t ii = 0; ii < phinat_begin; ++ii)
        x_phi_[input_dofs[ii]] = x[input_dofs[ii]];
      for (size_t ii = phinat_begin; ii < mu_begin; ++ii)
        x_phinat_[input_dofs[ii] - size_phi_] = x[input_dofs[ii]];
      for (size_t ii = mu_begin; ii < input_dofs.size(); ++ii)
        x_mu_[input_dofs[ii] - 2 * size_phi_] = x[input_dofs[ii]];
    } // if (!restricted)
    // apply matrices
    const auto mv = cellmodel_solver_->template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_->template axpy_func<Vector>(restricted_);
    // first row
    const auto& phi_dofs = cellmodel_solver_->phi_deim_output_dofs_[cell_];
    mv(M_, x_phi_, y_phi_, phi_dofs);
    mv(D_, x_phi_, tmp_vec_, phi_dofs);
    axpy(y_phi_, dt_, tmp_vec_, phi_dofs);
    mv(M_ell_, x_phinat_, tmp_vec_, phi_dofs);
    axpy(y_phi_, dt_ * gamma_, tmp_vec_, phi_dofs);
    for (const auto& DoF : dirichlet_.dirichlet_DoFs())
      y_phi_[DoF] = x_phi_[DoF];
    // second row
    const auto& phinat_dofs = cellmodel_solver_->phinat_deim_output_dofs_[cell_];
    mv(M_, x_phinat_, y_phinat_, phinat_dofs);
    mv(M_ell_, x_mu_, tmp_vec_, phinat_dofs);
    axpy(y_phinat_, 1. / Be_, tmp_vec_, phinat_dofs);
    // third row
    const auto& mu_dofs = cellmodel_solver_->mu_deim_output_dofs_[cell_];
    mv(M_, x_mu_, y_mu_, mu_dofs);
    mv(M_ell_, x_phi_, tmp_vec_, mu_dofs);
    axpy(y_mu_, epsilon_, tmp_vec_, mu_dofs);
    mv(A_boundary_, x_phi_, tmp_vec_, mu_dofs);
    axpy(y_mu_, -epsilon_, tmp_vec_, mu_dofs);
    // dirichlet_.apply(y_mu_);
    // copy to result vector
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi_; ++ii) {
        y[ii] = y_phi_[ii];
        y[size_phi_ + ii] = y_phinat_[ii];
        y[2 * size_phi_ + ii] = y_mu_[ii];
      }
    } else {
      for (const auto& dof : phi_dofs)
        y[dof] = y_phi_[dof];
      for (const auto& dof : phinat_dofs)
        y[size_phi_ + dof] = y_phinat_[dof];
      for (const auto& dof : mu_dofs)
        y[2 * size_phi_ + dof] = y_mu_[dof];
    } // if (!restricted)
  }

  void set_params(const XT::Common::Parameter& param) override final
  {
    gamma_ = param.get("gamma")[0];
    epsilon_ = param.get("epsilon")[0];
    Be_ = param.get("Be")[0];
  }

  void prepare(const double dt, const size_t cell, const bool restricted = false) override final
  {
    dt_ = dt;
    cell_ = cell;
    restricted_ = restricted;
  }

private:
  const Matrix& M_;
  const Matrix& D_;
  const Matrix& M_ell_;
  const Matrix& A_boundary_;
  const Dirichlet& dirichlet_;
  double gamma_;
  double epsilon_;
  double Be_;
  const CellModelSolverType* cellmodel_solver_;
  size_t cell_;
  bool restricted_;
  double dt_;
  const size_t size_phi_;
  // vectors to store intermediate results
  mutable Vector x_phi_;
  mutable Vector x_phinat_;
  mutable Vector x_mu_;
  mutable Vector y_phi_;
  mutable Vector y_phinat_;
  mutable Vector y_mu_;
  mutable Vector tmp_vec_;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_LINEAROPERATORS_HH
