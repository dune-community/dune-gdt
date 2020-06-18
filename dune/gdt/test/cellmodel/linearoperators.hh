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


template <class MatrixType, class VectorType>
class LinearOperatorWrapper : public Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>
{
  using BaseType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;

public:
  using Matrix = MatrixType;
  using Vector = VectorType;
  using Field = typename VectorType::ScalarType;

  LinearOperatorWrapper(const MatrixType& M, const size_t vector_length)
    : M_(M)
    , tmp_vec_(vector_length, 0.)
  {}

  using BaseType::apply;

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  virtual void set_params(const XT::Common::Parameter& /*param*/, const bool /*restricted*/ = false) {}

  virtual void prepare(const size_t /*cell*/, const bool /*restricted*/ = false) {}

  virtual const Matrix& getmat() const override
  {
    return M_;
  }

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override
  {
    return SolverCategory::Category::sequential;
  }

private:
  const MatrixType& M_;
  mutable VectorType tmp_vec_;
};


template <class VectorType, class MatrixType, class CellModelSolverType>
class OfieldMatrixLinearPartOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  OfieldMatrixLinearPartOperator(const Matrix& M,
                                 const Matrix& A,
                                 const Matrix& C_lin,
                                 const double dt,
                                 const double kappa,
                                 const CellModelSolverType* cellmodel_solver)
    : BaseType(M, 2 * M.rows())
    , M_(M)
    , A_(A)
    , C_lin_(C_lin)
    , dt_(dt)
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
    const auto& source_dofs = cellmodel_solver_->ofield_deim_source_dofs_[cell_][1];
    const auto& Pnat_begin = cellmodel_solver_->Pnat_deim_source_dofs_begin_[cell_];
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    if (!restricted_) {
      for (size_t ii = 0; ii < size_P_; ++ii) {
        x_P_[ii] = x[ii];
        x_Pnat_[ii] = x[size_P_ + ii];
      }
    } else {
      for (size_t ii = 0; ii < Pnat_begin; ++ii)
        x_P_[source_dofs[ii]] = x[source_dofs[ii]];
      for (size_t ii = Pnat_begin; ii < source_dofs.size(); ++ii)
        x_Pnat_[source_dofs[ii] - size_P_] = x[source_dofs[ii]];
    }
    // apply matrices
    const auto mv = cellmodel_solver_->template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_->template vector_axpy_func<Vector>(restricted_);
    const auto add = cellmodel_solver_->template add_func<Vector>(restricted_);
    // M+dtA
    const auto& P_dofs = cellmodel_solver_->P_deim_range_dofs_[cell_];
    mv(M_, x_P_, y_P_, P_dofs);
    mv(A_, x_P_, tmp_vec_, P_dofs);
    axpy(y_P_, dt_, tmp_vec_, P_dofs);
    // dt B and D
    const auto& Pnat_dofs = cellmodel_solver_->Pnat_deim_range_dofs_[cell_];
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

  void prepare(const size_t cell, const bool restricted = false) override final
  {
    cell_ = cell;
    restricted_ = restricted;
  }

  virtual const Matrix& getmat() const override final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return M_;
  }

private:
  const Matrix& M_;
  const Matrix& A_;
  const Matrix& C_lin_;
  const double dt_;
  const double kappa_;
  const CellModelSolverType* cellmodel_solver_;
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
class PfieldPhiMatrixOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

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
                          const SolverType& solver,
                          const double dt,
                          const double gamma,
                          const double epsilon,
                          const double Be,
                          const double Ca)
    : BaseType(M, M.rows())
    , M_(M)
    , D_(D)
    , M_ell_(M_ell)
    , G_(G)
    , M_nonlin_(M_nonlin)
    , solver_(solver)
    , A_(M_.backend(), 0)
    , J_(M_.backend(), 0)
    , dt_(dt)
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
    solver_.apply_inverse_mass_matrix(tmp_vec2_, tmp_vec_, &last_E_result_);
    last_E_result_ = tmp_vec_.backend();
    M_ell_.mv(tmp_vec_, tmp_vec2_);
    y.axpy(-dt_ * gamma_, tmp_vec2_);
  }

  void set_params(const XT::Common::Parameter& param, const bool /*restricted*/ = false) override final
  {
    gamma_ = param.get("gamma")[0];
    epsilon_ = param.get("epsilon")[0];
    Be_ = param.get("Be")[0];
    Ca_ = param.get("Ca")[0];
  }

  void prepare(const size_t /*cell*/, const bool /*restricted*/) override final
  {
    DUNE_THROW(Dune::NotImplemented,
               "This prepare has to be called after M_nonlin_changed (i.e. in each Newton iteration), currently it is "
               "only called once before starting the Newton scheme.");
    // precompute A and J (M_nonlin_ may have changed)
    // saves three mvs in each apply, can be dropped if memory is an issue
    J_.backend() = M_nonlin_.backend();
    J_ *= 1. / (Be_ * std::pow(epsilon_, 2));
    J_.axpy(1. / Ca_, M_);
    J_.axpy(1. / Be_, M_ell_);
    A_.backend() = M_nonlin_.backend();
    A_ *= 1. / epsilon_;
    A_.axpy(epsilon_, M_ell_);
  }

  virtual const Matrix& getmat() const override final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return M_;
  }

private:
  const Matrix& M_;
  const Matrix& D_;
  const Matrix& M_ell_;
  const Matrix& G_;
  const Matrix& M_nonlin_;
  const SolverType& solver_;
  MatrixType A_;
  MatrixType J_;
  const double dt_;
  double gamma_;
  double epsilon_;
  double Be_;
  double Ca_;
  // vectors to store intermediate results
  mutable Vector tmp_vec_;
  mutable Vector tmp_vec2_;
  mutable Vector last_A_result_;
  mutable Vector last_E_result_;
};

template <class VectorType, class MatrixType>
class MatrixToLinearOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  MatrixToLinearOperator(const Matrix& M)
    : BaseType(M, M.rows())
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

template <class VectorType, class MatrixType, class CellModelSolverType>
class PfieldMatrixLinearPartOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  PfieldMatrixLinearPartOperator(const Matrix& M,
                                 const Matrix& D,
                                 const Matrix& M_ell,
                                 const double dt,
                                 const double gamma,
                                 const double eps,
                                 const double Be,
                                 const CellModelSolverType* cellmodel_solver)
    : BaseType(M, 3 * M.rows())
    , M_(M)
    , D_(D)
    , M_ell_(M_ell)
    , dt_(dt)
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
    const auto& source_dofs = cellmodel_solver_->pfield_deim_source_dofs_[cell_][0];
    const auto& phinat_begin = cellmodel_solver_->phinat_deim_source_dofs_begin_[cell_];
    const auto& mu_begin = cellmodel_solver_->mu_deim_source_dofs_begin_[cell_];
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi_; ++ii) {
        x_phi_[ii] = x[ii];
        x_phinat_[ii] = x[size_phi_ + ii];
        x_mu_[ii] = x[2 * size_phi_ + ii];
      }
    } else {
      for (size_t ii = 0; ii < phinat_begin; ++ii)
        x_phi_[source_dofs[ii]] = x[source_dofs[ii]];
      for (size_t ii = phinat_begin; ii < mu_begin; ++ii)
        x_phinat_[source_dofs[ii] - size_phi_] = x[source_dofs[ii]];
      for (size_t ii = mu_begin; ii < source_dofs.size(); ++ii)
        x_mu_[source_dofs[ii] - 2 * size_phi_] = x[source_dofs[ii]];
    } // if (!restricted)
    // apply matrices
    const auto mv = cellmodel_solver_->template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_->template vector_axpy_func<Vector>(restricted_);
    const auto scal = cellmodel_solver_->template scal_func<Vector>(restricted_);
    // first row
    const auto& phi_dofs = cellmodel_solver_->phi_deim_range_dofs_[cell_];
    mv(M_, x_phi_, y_phi_, phi_dofs);
    mv(D_, x_phi_, tmp_vec_, phi_dofs);
    axpy(y_phi_, dt_, tmp_vec_, phi_dofs);
    mv(M_ell_, x_phinat_, tmp_vec_, phi_dofs);
    axpy(y_phi_, dt_ * gamma_, tmp_vec_, phi_dofs);
    // second row
    const auto& phinat_dofs = cellmodel_solver_->phinat_deim_range_dofs_[cell_];
    mv(M_, x_phinat_, y_phinat_, phinat_dofs);
    scal(y_phinat_, 1., phinat_dofs);
    mv(M_ell_, x_mu_, tmp_vec_, phinat_dofs);
    axpy(y_phinat_, 1. / Be_, tmp_vec_, phinat_dofs);
    // third row
    const auto& mu_dofs = cellmodel_solver_->mu_deim_range_dofs_[cell_];
    mv(M_, x_mu_, y_mu_, mu_dofs);
    mv(M_ell_, x_phi_, tmp_vec_, mu_dofs);
    axpy(y_mu_, epsilon_, tmp_vec_, mu_dofs);
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

  void set_params(const XT::Common::Parameter& param, const bool /*restricted*/ = false) override final
  {
    gamma_ = param.get("gamma")[0];
    epsilon_ = param.get("epsilon")[0];
    Be_ = param.get("Be")[0];
  }

  void prepare(const size_t cell, const bool restricted = false) override final
  {
    cell_ = cell;
    restricted_ = restricted;
  }

  virtual const Matrix& getmat() const override final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return M_;
  }

private:
  const Matrix& M_;
  const Matrix& D_;
  const Matrix& M_ell_;
  double dt_;
  double gamma_;
  double epsilon_;
  double Be_;
  const CellModelSolverType* cellmodel_solver_;
  size_t cell_;
  bool restricted_;
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


template <class VectorType, class MatrixType, bool use_incomplete_lut = false>
class Matrix2InverseOperator : public Dune::InverseOperator<VectorType, VectorType>
{
  using BaseType = Dune::InverseOperator<VectorType, VectorType>;

public:
  using Vector = VectorType;
  using Field = typename VectorType::ScalarType;
  using ColMajorMatrixType = ::Eigen::SparseMatrix<Field, ::Eigen::ColMajor>;
  using LUSolverType = ::Eigen::SparseLU<ColMajorMatrixType>;
  using IncompleteLUTSolverType = ::Eigen::IncompleteLUT<Field>;
  using SolverType = std::conditional_t<use_incomplete_lut, IncompleteLUTSolverType, LUSolverType>;

  Matrix2InverseOperator(const MatrixType& matrix)
    : matrix_(matrix)
    , solver_(std::make_shared<SolverType>())
  {}

  virtual void apply(VectorType& x, VectorType& b, InverseOperatorResult& /*res*/) override final
  {
    x.backend() = solver_->solve(b.backend());
  }

  virtual void apply(VectorType& x, VectorType& b, double /*reduction*/, InverseOperatorResult& res) override final
  {
    apply(x, b, res);
  }

  void prepare()
  {
    matrix_colmajor_ = matrix_.backend();
    matrix_colmajor_.makeCompressed();
    solver_->analyzePattern(matrix_colmajor_);
    solver_->factorize(matrix_colmajor_);
  }

  //! Category of the solver (see SolverCategory::Category)
  virtual SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

private:
  const MatrixType& matrix_;
  ColMajorMatrixType matrix_colmajor_;
  std::shared_ptr<SolverType> solver_;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_LINEAROPERATORS_HH
