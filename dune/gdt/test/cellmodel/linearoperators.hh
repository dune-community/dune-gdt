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
#include <dune/istl/solver.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/parameter.hh>

#include <Eigen/UmfPackSupport>
#include <Eigen/IterativeLinearSolvers>

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

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const final
  {
    auto& Sx = tmp_vec_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  virtual void prepare(const size_t /*cell*/, const bool /*restricted*/ = false) {}

  const Matrix& getmat() const override
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
  OfieldMatrixLinearPartOperator(const CellModelSolverType& cellmodel_solver)
    : BaseType(cellmodel_solver.M_ofield_, 2 * cellmodel_solver.M_ofield_.rows())
    , cellmodel_solver_(cellmodel_solver)
    , x_P_(cellmodel_solver.size_P_, 0., 0)
    , x_Pnat_(cellmodel_solver.size_P_, 0., 0)
    , y_P_(cellmodel_solver.size_P_, 0., 0)
    , y_Pnat_(cellmodel_solver.size_P_, 0., 0)
    , tmp_vec_(cellmodel_solver.size_P_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const final
  {
    // get some variables from cellmodel_solver_
    const auto dt = cellmodel_solver_.dt_;
    const auto kappa = cellmodel_solver_.kappa_;
    const auto Pa = cellmodel_solver_.Pa_;
    const auto size_P = cellmodel_solver_.size_P_;
    const auto& P_dofs = cellmodel_solver_.P_deim_range_dofs_[cell_];
    const auto& Pnat_dofs = cellmodel_solver_.Pnat_deim_range_dofs_[cell_];
    const auto mv = cellmodel_solver_.template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_.template vector_axpy_func<Vector>(restricted_);
    const auto add = cellmodel_solver_.template add_func<Vector>(restricted_);
    const MatrixType& M = cellmodel_solver_.M_ofield_;
    const MatrixType& E = cellmodel_solver_.E_ofield_;
    const MatrixType& B = cellmodel_solver_.B_ofield_;
    const MatrixType& C_incl_coeffs_and_sign = cellmodel_solver_.C_ofield_incl_coeffs_and_sign_;
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    if (!restricted_) {
      for (size_t ii = 0; ii < size_P; ++ii) {
        x_P_[ii] = x[ii];
        x_Pnat_[ii] = x[size_P + ii];
      }
    } else {
      const auto& source_dofs = cellmodel_solver_.ofield_deim_source_dofs_[cell_][1];
      const auto& Pnat_begin = cellmodel_solver_.Pnat_deim_source_dofs_begin_[cell_];
      for (size_t ii = 0; ii < Pnat_begin; ++ii)
        x_P_[source_dofs[ii]] = x[source_dofs[ii]];
      for (size_t ii = Pnat_begin; ii < source_dofs.size(); ++ii)
        x_Pnat_[source_dofs[ii] - size_P] = x[source_dofs[ii]];
    }
    // apply matrices
    // S_00
    mv(M, x_P_, y_P_, P_dofs);
    mv(B, x_P_, tmp_vec_, P_dofs);
    axpy(y_P_, dt, tmp_vec_, P_dofs);
    // S_01, S_11
    mv(M, x_Pnat_, y_Pnat_, P_dofs);
    if (restricted_)
      mv(M, x_Pnat_, y_Pnat_, Pnat_dofs);
    axpy(y_P_, dt / kappa, y_Pnat_, P_dofs);
    // linear part of S_10
    mv(C_incl_coeffs_and_sign, x_P_, tmp_vec_, Pnat_dofs);
    add(y_Pnat_, tmp_vec_, Pnat_dofs);
    mv(E, x_P_, tmp_vec_, Pnat_dofs);
    axpy(y_Pnat_, -1. / Pa, tmp_vec_, Pnat_dofs);
    // copy to result vector
    if (!restricted_) {
      for (size_t ii = 0; ii < size_P; ++ii) {
        y[ii] = y_P_[ii];
        y[size_P + ii] = y_Pnat_[ii];
      }
    } else {
      for (const auto& dof : P_dofs)
        y[dof] = y_P_[dof];
      for (const auto& dof : Pnat_dofs)
        y[size_P + dof] = y_Pnat_[dof];
    }
  }

  void prepare(const size_t cell, const bool restricted = false) final
  {
    cell_ = cell;
    restricted_ = restricted;
  }

  virtual const Matrix& getmat() const final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return cellmodel_solver_.M_ofield_;
  }

  const CellModelSolverType& cellmodel_solver_;
  size_t cell_;
  bool restricted_;
  // vectors to store intermediate results
  mutable Vector x_P_;
  mutable Vector x_Pnat_;
  mutable Vector y_P_;
  mutable Vector y_Pnat_;
  mutable Vector tmp_vec_;
};

template <class VectorType, class MatrixType, class CellModelSolverType>
class StokesMatrixLinearOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  StokesMatrixLinearOperator(const CellModelSolverType& cellmodel_solver)
    : BaseType(cellmodel_solver.A_stokes_, cellmodel_solver.size_u_ + cellmodel_solver.size_p_)
    , cellmodel_solver_(cellmodel_solver)
    , x_u_(cellmodel_solver.size_u_, 0., 0)
    , x_p_(cellmodel_solver.size_p_, 0., 0)
    , y_u_(cellmodel_solver.size_u_, 0., 0)
    , y_p_(cellmodel_solver.size_p_, 0., 0)
    , tmp_vec_u_(cellmodel_solver.size_u_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const final
  {
    // get some variables from cellmodel_solver_
    const auto size_u = cellmodel_solver_.size_u_;
    const auto size_p = cellmodel_solver_.size_p_;
    const auto& u_dofs = cellmodel_solver_.u_deim_range_dofs_;
    const auto& p_dofs = cellmodel_solver_.p_deim_range_dofs_;
    const auto mv = cellmodel_solver_.template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_.template vector_axpy_func<Vector>(restricted_);
    const auto add = cellmodel_solver_.template add_func<Vector>(restricted_);
    const MatrixType& A = cellmodel_solver_.A_stokes_;
    const MatrixType& BT = cellmodel_solver_.BT_stokes_;
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    Vector x_u_restricted_, x_p_restricted_, y_p_restricted_;
    if (!restricted_) {
      for (size_t ii = 0; ii < size_u; ++ii)
        x_u_[ii] = x[ii];
      for (size_t ii = 0; ii < size_p; ++ii)
        x_p_[ii] = x[size_u + ii];
    } else {
      const auto& source_dofs = cellmodel_solver_.stokes_deim_source_dofs_[2];
      const auto p_begin = cellmodel_solver_.p_deim_source_dofs_begin_;
      assert(x.size() >= p_begin);
      assert(x.size() == source_dofs.size());
      x_u_restricted_ = Vector(p_begin);
      x_p_restricted_ = Vector(source_dofs.size() - p_begin);
      y_p_restricted_ = Vector(p_dofs.size());
      for (size_t ii = 0; ii < p_begin; ++ii) {
        x_u_[source_dofs[ii]] = x[ii];
        x_u_restricted_[ii] = x[ii];
      }
      for (size_t ii = p_begin; ii < source_dofs.size(); ++ii) {
        x_p_[source_dofs[ii] - size_u] = x[ii];
        x_p_restricted_[ii - p_begin] = x[ii];
      }
    }
    // apply matrices
    // first row
    mv(A, x_u_, y_u_, u_dofs);
    mv(BT, x_p_, tmp_vec_u_, u_dofs);
    add(y_u_, tmp_vec_u_, u_dofs);
    if (restricted_)
      cellmodel_solver_.B_stokes_restricted_.mv(x_u_restricted_, y_p_restricted_);
    else
      BT.mtv(x_u_, y_p_);
    // copy to result vector
    if (!restricted_) {
      y_p_[0] = x_p_[0];
      for (size_t ii = 0; ii < size_u; ++ii)
        y[ii] = y_u_[ii];
      for (size_t ii = 0; ii < size_p; ++ii)
        y[size_u + ii] = y_p_[ii];
    } else {
      for (const auto& dof : u_dofs)
        y[dof] = y_u_[dof];
      for (size_t ii = 0; ii < p_dofs.size(); ++ii)
        y[size_u + p_dofs[ii]] = (p_dofs[ii] == 0) ? x_p_restricted_[0] : y_p_restricted_[ii];
    }
  }

  void prepare(const size_t /*cell*/, const bool restricted = false) final
  {
    restricted_ = restricted;
  }

  const Matrix& getmat() const final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return cellmodel_solver_.A_stokes_;
  }

private:
  const CellModelSolverType& cellmodel_solver_;
  bool restricted_;
  // vectors to store intermediate results
  mutable Vector x_u_;
  mutable Vector x_p_;
  mutable Vector y_u_;
  mutable Vector y_p_;
  mutable Vector tmp_vec_u_;
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
  void apply(const Vector& x, Vector& y) const final
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
  PfieldMatrixLinearPartOperator(const CellModelSolverType& cellmodel_solver, const bool exclude_mass_matrices = false)
    : BaseType(cellmodel_solver.M_pfield_, 3 * cellmodel_solver.M_pfield_.rows())
    , cellmodel_solver_(cellmodel_solver)
    , x_phi_(cellmodel_solver.size_phi_, 0., 0)
    , x_phinat_(cellmodel_solver.size_phi_, 0., 0)
    , x_mu_(cellmodel_solver.size_phi_, 0., 0)
    , y_phi_(cellmodel_solver.size_phi_, 0., 0)
    , y_phinat_(cellmodel_solver.size_phi_, 0., 0)
    , y_mu_(cellmodel_solver.size_phi_, 0., 0)
    , tmp_vec_(cellmodel_solver.size_phi_, 0., 0)
    , exclude_mass_matrices_(exclude_mass_matrices)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const final
  {
    // get some variables from cellmodel_solver_
    const MatrixType& M = cellmodel_solver_.M_pfield_;
    const MatrixType& B = cellmodel_solver_.B_pfield_;
    const MatrixType& E = cellmodel_solver_.E_pfield_;
    const auto size_phi = cellmodel_solver_.size_phi_;
    const auto dt = cellmodel_solver_.dt_;
    const auto gamma = cellmodel_solver_.gamma_;
    const auto epsilon = cellmodel_solver_.epsilon_;
    const auto Be = cellmodel_solver_.Be_;
    const auto Ca = cellmodel_solver_.Ca_;
    const auto mv = cellmodel_solver_.template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_.template vector_axpy_func<Vector>(restricted_);
    const auto scal = cellmodel_solver_.template scal_func<Vector>(restricted_);
    const auto& phi_dofs = cellmodel_solver_.phi_deim_range_dofs_[cell_];
    const auto& phinat_dofs = cellmodel_solver_.phinat_deim_range_dofs_[cell_];
    const auto& mu_dofs = cellmodel_solver_.mu_deim_range_dofs_[cell_];
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi; ++ii) {
        x_phi_[ii] = x[ii];
        x_phinat_[ii] = x[size_phi + ii];
        x_mu_[ii] = x[2 * size_phi + ii];
      }
    } else {
      const auto& source_dofs = cellmodel_solver_.pfield_deim_source_dofs_[cell_][0];
      const auto& phinat_begin = cellmodel_solver_.phinat_deim_source_dofs_begin_[cell_];
      const auto& mu_begin = cellmodel_solver_.mu_deim_source_dofs_begin_[cell_];
      for (size_t ii = 0; ii < phinat_begin; ++ii)
        x_phi_[source_dofs[ii]] = x[source_dofs[ii]];
      for (size_t ii = phinat_begin; ii < mu_begin; ++ii)
        x_phinat_[source_dofs[ii] - size_phi] = x[source_dofs[ii]];
      for (size_t ii = mu_begin; ii < source_dofs.size(); ++ii)
        x_mu_[source_dofs[ii] - 2 * size_phi] = x[source_dofs[ii]];
    } // if (!restricted)
    // apply matrices
    // first row
    if (!exclude_mass_matrices_)
      mv(M, x_phi_, y_phi_, phi_dofs);
    else
      scal(y_phi_, 0, phi_dofs);
    mv(B, x_phi_, tmp_vec_, phi_dofs);
    axpy(y_phi_, -dt, tmp_vec_, phi_dofs);
    mv(E, x_phinat_, tmp_vec_, phi_dofs);
    axpy(y_phi_, dt * gamma, tmp_vec_, phi_dofs);
    // second row
    if (!exclude_mass_matrices_)
      mv(M, x_phinat_, y_phinat_, phinat_dofs);
    else
      scal(y_phinat_, 0, phinat_dofs);
    mv(E, x_mu_, tmp_vec_, phinat_dofs);
    axpy(y_phinat_, 1. / Be, tmp_vec_, phinat_dofs);
    // third row (and M * mu from second row)
    // even if we exclude the mass matrices, we apply the mass matrix first...
    mv(M, x_mu_, y_mu_, mu_dofs);
    if (restricted_)
      mv(M, x_mu_, y_mu_, phinat_dofs);
    // ... to use the result here for the second row, ...
    axpy(y_phinat_, 1. / Ca, y_mu_, phinat_dofs);
    if (exclude_mass_matrices_) {
      // ...and set y_mu to 0 again here
      scal(y_mu_, 0, mu_dofs);
    }
    mv(E, x_phi_, tmp_vec_, mu_dofs);
    axpy(y_mu_, epsilon, tmp_vec_, mu_dofs);
    // copy to result vector
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi; ++ii) {
        y[ii] = y_phi_[ii];
        y[size_phi + ii] = y_phinat_[ii];
        y[2 * size_phi + ii] = y_mu_[ii];
      }
    } else {
      for (const auto& dof : phi_dofs)
        y[dof] = y_phi_[dof];
      for (const auto& dof : phinat_dofs)
        y[size_phi + dof] = y_phinat_[dof];
      for (const auto& dof : mu_dofs)
        y[2 * size_phi + dof] = y_mu_[dof];
    } // if (!restricted)
  }

  void prepare(const size_t cell, const bool restricted = false) final
  {
    cell_ = cell;
    restricted_ = restricted;
  }

  virtual const Matrix& getmat() const final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return cellmodel_solver_.M_pfield_;
  }

private:
  const CellModelSolverType& cellmodel_solver_;
  size_t cell_;
  bool restricted_;
  // vectors to store intermediate results
  mutable Vector x_phi_;
  mutable Vector x_phinat_;
  mutable Vector x_mu_;
  mutable Vector y_phi_;
  mutable Vector y_phinat_;
  mutable Vector y_mu_;
  mutable Vector tmp_vec_;
  const bool exclude_mass_matrices_;
};

template <class VectorType, class MatrixType, class CellModelSolverType>
class PfieldDiagonalMassMatricesOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  PfieldDiagonalMassMatricesOperator(const CellModelSolverType& cellmodel_solver)
    : BaseType(cellmodel_solver.M_pfield_, 3 * cellmodel_solver.M_pfield_.rows())
    , cellmodel_solver_(cellmodel_solver)
    , x_phi_(cellmodel_solver.size_phi_, 0., 0)
    , x_phinat_(cellmodel_solver.size_phi_, 0., 0)
    , x_mu_(cellmodel_solver.size_phi_, 0., 0)
    , y_phi_(cellmodel_solver.size_phi_, 0., 0)
    , y_phinat_(cellmodel_solver.size_phi_, 0., 0)
    , y_mu_(cellmodel_solver.size_phi_, 0., 0)
    , tmp_vec_(cellmodel_solver.size_phi_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const final
  {
    // get some variables from cellmodel_solver_
    const MatrixType& M = cellmodel_solver_.M_pfield_;
    const auto size_phi = cellmodel_solver_.size_phi_;
    const auto mv = cellmodel_solver_.template mv_func<Vector>(restricted_);
    const auto axpy = cellmodel_solver_.template vector_axpy_func<Vector>(restricted_);
    const auto scal = cellmodel_solver_.template scal_func<Vector>(restricted_);
    const auto& phi_dofs = cellmodel_solver_.phi_deim_range_dofs_[cell_];
    const auto& phinat_dofs = cellmodel_solver_.phinat_deim_range_dofs_[cell_];
    const auto& mu_dofs = cellmodel_solver_.mu_deim_range_dofs_[cell_];
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi; ++ii) {
        x_phi_[ii] = x[ii];
        x_phinat_[ii] = x[size_phi + ii];
        x_mu_[ii] = x[2 * size_phi + ii];
      }
    } else {
      const auto& source_dofs = cellmodel_solver_.pfield_deim_source_dofs_[cell_][0];
      const auto& phinat_begin = cellmodel_solver_.phinat_deim_source_dofs_begin_[cell_];
      const auto& mu_begin = cellmodel_solver_.mu_deim_source_dofs_begin_[cell_];
      for (size_t ii = 0; ii < phinat_begin; ++ii)
        x_phi_[source_dofs[ii]] = x[source_dofs[ii]];
      for (size_t ii = phinat_begin; ii < mu_begin; ++ii)
        x_phinat_[source_dofs[ii] - size_phi] = x[source_dofs[ii]];
      for (size_t ii = mu_begin; ii < source_dofs.size(); ++ii)
        x_mu_[source_dofs[ii] - 2 * size_phi] = x[source_dofs[ii]];
    } // if (!restricted)
    // apply matrices
    // first row
    mv(M, x_phi_, y_phi_, phi_dofs);
    mv(M, x_phinat_, y_phinat_, phinat_dofs);
    mv(M, x_mu_, y_mu_, mu_dofs);
    // copy to result vector
    if (!restricted_) {
      for (size_t ii = 0; ii < size_phi; ++ii) {
        y[ii] = y_phi_[ii];
        y[size_phi + ii] = y_phinat_[ii];
        y[2 * size_phi + ii] = y_mu_[ii];
      }
    } else {
      for (const auto& dof : phi_dofs)
        y[dof] = y_phi_[dof];
      for (const auto& dof : phinat_dofs)
        y[size_phi + dof] = y_phinat_[dof];
      for (const auto& dof : mu_dofs)
        y[2 * size_phi + dof] = y_mu_[dof];
    } // if (!restricted)
  }

  void prepare(const size_t cell, const bool restricted = false) final
  {
    cell_ = cell;
    restricted_ = restricted;
  }

  virtual const Matrix& getmat() const final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return cellmodel_solver_.M_pfield_;
  }

private:
  const CellModelSolverType& cellmodel_solver_;
  size_t cell_;
  bool restricted_;
  // vectors to store intermediate results
  mutable Vector x_phi_;
  mutable Vector x_phinat_;
  mutable Vector x_mu_;
  mutable Vector y_phi_;
  mutable Vector y_phinat_;
  mutable Vector y_mu_;
  mutable Vector tmp_vec_;
};

// For a saddle point matrix (A B^T; B C) this models the Schur complement (B A^{-1} B^T - C)
// We assume that A is spd and C has a single entry -1 in the upper left corner (to fix the value of p to 0 at that DoF)
template <class VectorType, class MatrixType, class CellModelSolverType>
class StokesSchurComplementOperator : public LinearOperatorWrapper<MatrixType, VectorType>
{
  using BaseType = LinearOperatorWrapper<MatrixType, VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  StokesSchurComplementOperator(const CellModelSolverType& cellmodel_solver)
    : BaseType(cellmodel_solver.A_stokes_, cellmodel_solver.size_u_)
    , cellmodel_solver_(cellmodel_solver)
    , tmp_vec_u1_(cellmodel_solver.size_u_)
    , tmp_vec_u2_(cellmodel_solver.size_u_)
  {}

  void apply(const Vector& x, Vector& y) const final
  {
    // we want to calculate y = (B A^{-1} B^T - C) x
    // calculate B^T x
    auto& BT_x = tmp_vec_u1_;
    cellmodel_solver_.BT_stokes_.mv(x, BT_x);
    // calculate A^{-1} B1 x
    auto& Ainv_BT_x = tmp_vec_u2_;
    Ainv_BT_x.backend() = cellmodel_solver_.stokes_A_solver_->solve(BT_x.backend());
    // apply B
    cellmodel_solver_.BT_stokes_.mtv(Ainv_BT_x, y);
    // add -Cx
    y[0] = x[0];
  }

  const Matrix& getmat() const final
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return cellmodel_solver_.A_stokes_;
  }

private:
  const CellModelSolverType& cellmodel_solver_;
  // vectors to store intermediate results
  mutable Vector tmp_vec_u1_;
  mutable Vector tmp_vec_u2_;
};


template <class VectorType, class MatrixType, bool use_incomplete_lut = false>
class Matrix2InverseOperator : public Dune::InverseOperator<VectorType, VectorType>
{
  using BaseType = Dune::InverseOperator<VectorType, VectorType>;

public:
  using Vector = VectorType;
  using Field = typename VectorType::ScalarType;
  using ColMajorMatrixType = ::Eigen::SparseMatrix<Field, ::Eigen::ColMajor>;
  using LUSolverType = ::Eigen::UmfPackLU<ColMajorMatrixType>;
  using IncompleteLUTSolverType = ::Eigen::IncompleteLUT<Field>;
  using SolverType = std::conditional_t<use_incomplete_lut, IncompleteLUTSolverType, LUSolverType>;

  Matrix2InverseOperator(const std::shared_ptr<MatrixType>& matrix)
    : matrix_(matrix)
    , solver_(std::make_shared<SolverType>())
  {}

  virtual void apply(VectorType& x, VectorType& b, InverseOperatorResult& /*res*/) final
  {
    x.backend() = solver_->solve(b.backend());
  }

  virtual void apply(VectorType& x, VectorType& b, double /*reduction*/, InverseOperatorResult& res) final
  {
    apply(x, b, res);
  }

  void prepare()
  {
    if constexpr (use_incomplete_lut) {
      // solver_->setDroptol(1e-6);
      solver_->setFillfactor(DXTC_CONFIG_GET("pfield_ilu_fillfactor", 80));
    }
    solver_->compute(matrix_->backend());
  }

  //! Category of the solver (see SolverCategory::Category)
  virtual SolverCategory::Category category() const final
  {
    return SolverCategory::Category::sequential;
  }

private:
  const std::shared_ptr<MatrixType>& matrix_;
  std::shared_ptr<SolverType> solver_;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_LINEAROPERATORS_HH
