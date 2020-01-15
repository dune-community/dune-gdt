// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_EXAMPLES_CELLMODEL_HH
#define DUNE_GDT_EXAMPLES_CELLMODEL_HH

#include <chrono>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-market.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/solver/istl/preconditioners.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/filters/element.hh>
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
#include <dune/gdt/local/integrands/symmetrized-laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/gradient-value.hh>
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
  static const size_t d = DiscreteFunc::d;
  const auto& grid_view = func.space().grid_view();
  const auto local_func = func.local_function();
  // write one file per MPI rank
  std::ofstream rankfile(get_rankfilename(prefix, func, step, grid_view.comm().rank()));
  for (const auto& entity : elements(grid_view, Partitions::interiorBorder)) {
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

template <class VectorType, class MatrixType, class CellModelSolverType>
class OfieldMatrixLinearPartOperator : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

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
    : M_(M)
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
    , tmp_vec2_(2 * size_P_, 0., 0)
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
    // M+dtA
    M_.mv(x_P_, y_P_);
    A_.mv(x_P_, tmp_vec_);
    y_P_.axpy(dt_, tmp_vec_);
    // dt B and D
    M_.mv(x_Pnat_, y_Pnat_);
    y_P_.axpy(dt_ / kappa_, y_Pnat_);
    // linear part of C
    C_lin_.mv(x_P_, tmp_vec_);
    y_Pnat_ += tmp_vec_;
    // copy to result vector
    const auto& P_dofs = cellmodel_solver_->P_deim_output_dofs_[cell_];
    const auto& Pnat_dofs = cellmodel_solver_->Pnat_deim_output_dofs_[cell_];
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

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec2_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

  void prepare(const double dt, const size_t cell, const bool restricted = false)
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
  mutable Vector tmp_vec2_;
};

template <class VectorType, class MatrixType>
class OfieldSchurComplementOperator : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  OfieldSchurComplementOperator(const Matrix& S)
    : S_(S)
    , tmp_vec_(S_.rows(), 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    S_.mv(x, y);
  }

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

private:
  const Matrix& S_;
  mutable Vector tmp_vec_;
};

template <class VectorType, class MatrixType, class SolverType, class DirichletConstraintsType>
class PfieldPhiMuMatrixOperator : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

public:
  using Vector = VectorType;
  using Matrix = MatrixType;
  using Dirichlet = DirichletConstraintsType;
  using Field = typename VectorType::ScalarType;

  // Matrix dimensions are
  // A: m x m, B1, B2: m x n, C: n x n
  PfieldPhiMuMatrixOperator(const Matrix& M,
                            const SolverType& M_inv,
                            const Matrix& D,
                            const Matrix& M_ell,
                            const Matrix& G,
                            const Matrix& M_nonlin,
                            const Dirichlet& dirichlet,
                            const double gamma,
                            const double eps,
                            const double Be,
                            const double Ca)
    : M_(M)
    , M_inv_(M_inv)
    , D_(D)
    , M_ell_(M_ell)
    , G_(G)
    , M_nonlin_(M_nonlin)
    , A_(M_.backend(), 0)
    , J_(M_.backend(), 0)
    , dirichlet_(dirichlet)
    , gamma_(gamma)
    , eps_(eps)
    , Be_(Be)
    , Ca_(Ca)
    , size_phi_(M_.rows())
    , x_phi_(size_phi_, 0., 0)
    , x_mu_(size_phi_, 0., 0)
    , y_phi_(size_phi_, 0., 0)
    , y_mu_(size_phi_, 0., 0)
    , tmp_vec_(size_phi_, 0., 0)
    , tmp_vec2_(size_phi_, 0., 0)
    , tmp_vec3_(2 * size_phi_, 0., 0)
  {}

  /*! \brief apply operator to x:  \f$ y = S(x) \f$
        The input vector is consistent and the output must also be
     consistent on the interior+border partition.
   */
  void apply(const Vector& x, Vector& y) const override final
  {
    // copy to temporary vectors (we do not use vector views to improve performance of mv)
    for (size_t ii = 0; ii < size_phi_; ++ii) {
      x_phi_[ii] = x[ii];
      x_mu_[ii] = x[ii + size_phi_];
    }
    // apply matrices
    // apply M + dt D
    M_.mv(x_phi_, y_phi_);
    D_.mv(x_phi_, tmp_vec_);
    y_phi_.axpy(dt_, tmp_vec_);
    // compute H^{-1} (G phi + J mu)
    G_.mv(x_phi_, tmp_vec_);
    J_.mv(x_mu_, tmp_vec2_);
    tmp_vec_ += tmp_vec2_;
    tmp_vec_.backend() = M_inv_.solve(tmp_vec_.backend());
    // apply dt E
    M_ell_.mv(tmp_vec_, tmp_vec2_);
    y_phi_.axpy(-dt_ * gamma_, tmp_vec2_);
    // apply A
    A_.mv(x_phi_, y_mu_);
    // apply C
    M_.mv(x_mu_, tmp_vec_);
    y_mu_ += tmp_vec_;
    // apply constraints
    for (const auto& DoF : dirichlet_.dirichlet_DoFs())
      y_mu_[DoF] = x_phi_[DoF];
    // copy to result vector
    for (size_t ii = 0; ii < size_phi_; ++ii) {
      y[ii] = y_phi_[ii];
      y[ii + size_phi_] = y_mu_[ii];
    }
  }

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec3_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

  void set_params(const double gamma, const double eps, const double Be, const double Ca)
  {
    gamma_ = gamma;
    eps_ = eps;
    Be_ = Be;
    Ca_ = Ca;
  }

  void prepare(const double dt)
  {
    dt_ = dt;
    // precompute A and J (M_nonlin_ may have changed)
    // saves three mvs in each apply, can be dropped if memory is an issue
    J_.backend() = M_nonlin_.backend();
    J_ *= 1. / (Be_ * std::pow(eps_, 2));
    J_.axpy(1. / Ca_, M_);
    J_.axpy(1. / Be_, M_ell_);
    A_.backend() = M_nonlin_.backend();
    A_ *= 1. / eps_;
    A_.axpy(eps_, M_ell_);
  }

private:
  const Matrix& M_;
  const SolverType& M_inv_;
  const Matrix& D_;
  const Matrix& M_ell_;
  const Matrix& G_;
  const Matrix& M_nonlin_;
  MatrixType A_;
  MatrixType J_;
  const Dirichlet& dirichlet_;
  double gamma_;
  double eps_;
  double Be_;
  double Ca_;
  double dt_;
  const size_t size_phi_;
  // vectors to store intermediate results
  mutable Vector x_phi_;
  mutable Vector x_mu_;
  mutable Vector y_phi_;
  mutable Vector y_mu_;
  mutable Vector tmp_vec_;
  mutable Vector tmp_vec2_;
  mutable Vector tmp_vec3_;
};


template <class VectorType, class MatrixType, class DirichletConstraintsType, class CellModelSolverType>
class PfieldMatrixLinearPartOperator : public Dune::LinearOperator<VectorType, VectorType>
{
  using BaseType = Dune::LinearOperator<VectorType, VectorType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

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
                                 const Dirichlet& dirichlet,
                                 const double gamma,
                                 const double eps,
                                 const double Be,
                                 const CellModelSolverType* cellmodel_solver)
    : M_(M)
    , D_(D)
    , M_ell_(M_ell)
    , dirichlet_(dirichlet)
    , gamma_(gamma)
    , eps_(eps)
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
    , tmp_vec2_(3 * size_phi_, 0., 0)
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
    // second row
    const auto& phinat_dofs = cellmodel_solver_->phinat_deim_output_dofs_[cell_];
    mv(M_, x_phinat_, y_phinat_, phinat_dofs);
    mv(M_ell_, x_mu_, tmp_vec_, phinat_dofs);
    axpy(y_phinat_, 1. / Be_, tmp_vec_, phinat_dofs);
    // third row
    const auto& mu_dofs = cellmodel_solver_->mu_deim_output_dofs_[cell_];
    mv(M_, x_mu_, y_mu_, mu_dofs);
    mv(M_ell_, x_phi_, tmp_vec_, mu_dofs);
    axpy(y_mu_, eps_, tmp_vec_, mu_dofs);
    dirichlet_.apply(y_mu_);
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

  void applyscaleadd(Field alpha, const Vector& x, Vector& y) const override final
  {
    auto& Sx = tmp_vec2_;
    apply(x, Sx);
    y.axpy(alpha, Sx);
  }

  void set_params(const double gamma, const double eps, const double Be)
  {
    gamma_ = gamma;
    eps_ = eps;
    Be_ = Be;
  }

  //! Category of the linear operator (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return SolverCategory::Category::sequential;
  }

  void prepare(const double dt, const size_t cell, const bool restricted = false)
  {
    dt_ = dt;
    cell_ = cell;
    restricted_ = restricted;
  }

private:
  const Matrix& M_;
  const Matrix& D_;
  const Matrix& M_ell_;
  const Dirichlet& dirichlet_;
  double gamma_;
  double eps_;
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
  mutable Vector tmp_vec2_;
};


struct CellModelSolver
{
  // using G = ALU_2D_SIMPLEX_CONFORMING;
  using G = YASP_2D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  using GV = typename G::LeafGridView;
  using PGV = XT::Grid::PeriodicGridView<GV>;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using PI = XT::Grid::extract_intersection_t<PGV>;
  using MatrixType = XT::LA::EigenRowMajorSparseMatrix<double>;
  using DenseMatrixType = XT::LA::EigenDenseMatrix<double>;
  using VectorType = XT::LA::CommonDenseVector<double>;
  using EigenVectorType = XT::LA::EigenDenseVector<double>;
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using EigenVectorViewType = XT::LA::VectorView<EigenVectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;
  using ConstEigenVectorViewType = XT::LA::ConstVectorView<EigenVectorType>;
  using R = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, PGV, 1, 1, R>;
  using LocalDiscreteFunctionType = typename DiscreteFunctionType::LocalFunctionType;
  using VectorDiscreteFunctionType = DiscreteFunction<VectorType, PGV, d, 1, R>;
  using VectorLocalDiscreteFunctionType = typename VectorDiscreteFunctionType::LocalFunctionType;
  using ViewDiscreteFunctionType = DiscreteFunction<VectorViewType, PGV, 1, 1, R>;
  using ViewLocalDiscreteFunctionType = typename ViewDiscreteFunctionType::LocalFunctionType;
  using ViewVectorDiscreteFunctionType = DiscreteFunction<VectorViewType, PGV, d, 1, R>;
  using ViewVectorLocalDiscreteFunctionType = typename ViewVectorDiscreteFunctionType::LocalFunctionType;
  using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;
  using DomainRetType = XT::Common::FieldVector<R, d>;
  using JacobianRetType = XT::Common::FieldMatrix<R, d, d>;
  using ColMajorBackendType = ::Eigen::SparseMatrix<R, ::Eigen::ColMajor>;
  using RowMajorBackendType = typename MatrixType::BackendType;
  using DirectSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using SolverType = ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::IncompleteLUT<R>>;
  using DiagSolverType = ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::DiagonalPreconditioner<R>>;
  using PhiDirichletConstraintsType = DirichletConstraints<PI, SpaceInterface<PGV, 1, 1, R>>;
  using PfieldPhiMuMatrixOperatorType =
      PfieldPhiMuMatrixOperator<EigenVectorType, MatrixType, DirectSolverType, PhiDirichletConstraintsType>;
  using PfieldPhiMuSolverType = Dune::RestartedGMResSolver<EigenVectorType>;
  // using PfieldPhiMuSolverType = Dune::BiCGSTABSolver<EigenVectorType>;
  using OfieldSchurSolverType = Dune::RestartedGMResSolver<EigenVectorType>;
  using PerThreadVectorLocalFunc = XT::Common::PerThreadValue<std::unique_ptr<VectorLocalDiscreteFunctionType>>;
  using PerThreadScalarLocalFuncs = XT::Common::PerThreadValue<std::vector<std::unique_ptr<LocalDiscreteFunctionType>>>;
  using PerThreadVectorLocalFuncs =
      XT::Common::PerThreadValue<std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>>>;

  CellModelSolver(const std::string testcase = "single_cell",
                  const double t_end = 1.,
                  const unsigned int num_elements_x = 50,
                  const unsigned int num_elements_y = 50,
                  const bool use_tbb = true,
                  const double Be = 0.3, // bending capillary number, ratio of viscous forces to bending forces
                  const double Ca = 0.1, // capillary number, ratio of viscous forces to surface tension forces
                  const double Pa = 1, // polarization elasticity number
                  const double Re = 5e-13, // Reynolds number
                  const double Fa = 1., // active force number
                  const double xi = 1.1, // alignment of P with the flow, > 0 for rod-like cells and < 0 for oblate ones
                  const double kappa = 1.65, // eta_rot/eta, scaling factor between rotational and dynamic viscosity
                  const double c_1 = 5., // double well shape parameter
                  const double beta = 0., // alignment of P with the boundary of cell
                  const double gamma = 0.025, // phase field mobility coefficient
                  const double epsilon = 0.21, // phase field parameter
                  const double In = 1., // interaction parameter
                  const bool linearize = false,
                  const double pol_order = 2)
    : lower_left_(get_lower_left(testcase))
    , upper_right_(get_upper_right(testcase))
    , t_end_(t_end)
    , t_(0.)
    , use_tbb_(use_tbb)
    , Re_(Re)
    , Fa_inv_(1. / Fa)
    , xi_(xi)
    , kappa_(kappa)
    , c_1_(c_1)
    , Pa_(Pa)
    , last_pfield_Pa_(Pa_)
    , last_ofield_Pa_(Pa_)
    , beta_(beta)
    , gamma_(gamma)
    , Be_(Be)
    , Ca_(Ca)
    , epsilon_(epsilon)
    , In_(In)
    , vol_domain_((upper_right_[0] - lower_left_[0]) * (upper_right_[1] - lower_left_[1]))
    , num_cells_(get_num_cells(testcase))
    , linearize_(linearize)
    // do a global refine once, this makes simplicial grids look more symmetric
    , grid_(XT::Grid::make_cube_grid<G>(lower_left_, upper_right_, {num_elements_x, num_elements_y}, 1))
    , nonperiodic_grid_view_(grid_.leaf_view())
    , grid_view_(nonperiodic_grid_view_, std::bitset<d>(get_periodic_directions(testcase)))
    , u_space_(make_continuous_lagrange_space<d>(grid_view_, pol_order))
    , p_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order - 1))
    , phi_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order))
    , size_u_(u_space_.mapper().size())
    , size_p_(p_space_.mapper().size())
    , size_phi_(phi_space_.mapper().size())
    , num_mutexes_u_(use_tbb ? 0 : size_u_ / 100)
    , num_mutexes_ofield_(use_tbb ? 0 : size_u_ / 100)
    , num_mutexes_pfield_(use_tbb ? 0 : size_phi_ / 100)
    , stokes_vector_(size_u_ + size_p_, 0., 0)
    , ofield_vectors_(num_cells_, VectorType(2 * size_u_, 0., 0))
    , pfield_vectors_(num_cells_, VectorType(3 * size_phi_, 0., 0))
    , u_view_(stokes_vector_, 0, size_u_)
    , p_view_(stokes_vector_, size_u_, size_u_ + size_p_)
    , u_(u_space_, u_view_, "u")
    , p_(p_space_, p_view_, "p")
    , S_stokes_(size_u_ + size_p_, size_u_ + size_p_, create_stokes_pattern(u_space_, p_space_), 100)
    , A_stokes_(S_stokes_, 0, size_u_, 0, size_u_)
    , B_stokes_(S_stokes_, 0, size_u_, size_u_, size_u_ + size_p_)
    , BT_stokes_(S_stokes_, size_u_, size_u_ + size_p_, 0, size_u_)
    , M_p_stokes_(size_p_, size_p_, make_element_sparsity_pattern(p_space_, p_space_, grid_view_), 100)
    , A_stokes_op_(std::make_shared<MatrixOperator<MatrixViewType, PGV, d>>(grid_view_, u_space_, u_space_, A_stokes_))
    , stokes_rhs_vector_(size_u_ + size_p_, 0., num_mutexes_u_)
    , stokes_f_vector_(stokes_rhs_vector_, 0, size_u_)
    , stokes_g_vector_(stokes_rhs_vector_, size_u_, size_u_ + size_p_)
    , p_basis_integrated_vector_(size_p_)
    , u_dirichlet_constraints_(make_dirichlet_constraints(u_space_, boundary_info_))
    , phi_dirichlet_constraints_(make_dirichlet_constraints(phi_space_, boundary_info_))
    , ofield_submatrix_pattern_(make_element_sparsity_pattern(u_space_, u_space_, grid_view_))
    // , S_ofield_(2 * size_u_, 2 * size_u_, create_ofield_pattern(size_u_, ofield_submatrix_pattern_), 100)
    // , S_ofield_00_(S_ofield_, 0, size_u_, 0, size_u_)
    // , S_ofield_01_(S_ofield_, 0, size_u_, size_u_, 2 * size_u_)
    // , S_ofield_10_(S_ofield_, size_u_, 2 * size_u_, 0, size_u_)
    // , S_ofield_11_(S_ofield_, size_u_, 2 * size_u_, size_u_, 2 * size_u_)
    , M_ofield_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
    , A_ofield_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
    , C_ofield_elliptic_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
    , C_ofield_linear_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
    , C_ofield_nonlinear_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
    , S_schur_ofield_linear_part_(size_u_, size_u_, ofield_submatrix_pattern_, 0)
    , S_schur_ofield_(size_u_, size_u_, ofield_submatrix_pattern_, 0)
    , M_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, M_ofield_))
    , A_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, A_ofield_))
    , C_ofield_linear_part_op_(
          std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, C_ofield_linear_part_))
    , C_ofield_nonlinear_part_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(
          grid_view_, u_space_, u_space_, C_ofield_nonlinear_part_))
    , ofield_jac_linear_op_(M_ofield_, A_ofield_, C_ofield_linear_part_, kappa_, this)
    , ofield_schur_op_(S_schur_ofield_)
    , ofield_rhs_vector_(2 * size_u_, 0., num_mutexes_ofield_)
    , ofield_f_vector_(ofield_rhs_vector_, 0, size_u_)
    , ofield_g_vector_(ofield_rhs_vector_, size_u_, 2 * size_u_)
    , P_eigen_(num_cells_, EigenVectorType(size_u_, 0., 0))
    , P_tmp_eigen_(size_u_, 0., 0)
    , P_tmp_eigen2_(size_u_, 0., 0)
    // , ofield_old_result_(num_cells_, EigenVectorType(2 * size_u_, 0.))
    , stokes_solver_(std::make_shared<DirectSolverType>())
    , ofield_solver_(std::make_shared<SolverType>())
    , ofield_mass_matrix_solver_(std::make_shared<DirectSolverType>())
    , identity_prec_(SolverCategory::Category::sequential)
    , ofield_schur_solver_(
          std::make_shared<OfieldSchurSolverType>(ofield_schur_op_, identity_prec_, 1e-10, 20, 10000, false))
    , ofield_tmp_vec_(2 * size_u_, 0., 0)
    // , ofield_schur_solver_(std::make_shared<DiagSolverType>())
    , ofield_deim_input_dofs_(num_cells_)
    , Pnat_deim_input_dofs_begin_(num_cells_)
    , ofield_deim_output_dofs_(num_cells_)
    , ofield_deim_unique_output_dofs_(num_cells_)
    , P_deim_output_dofs_(num_cells_)
    , Pnat_deim_output_dofs_(num_cells_)
    , ofield_deim_entities_(num_cells_)
    , pfield_submatrix_pattern_(make_element_sparsity_pattern(phi_space_, phi_space_, grid_view_))
    // , S_pfield_(3 * size_phi_, 3 * size_phi_, create_pfield_pattern(size_phi_, pfield_submatrix_pattern_), 100)
    // , S_pfield_00_(S_pfield_, 0, size_phi_, 0, size_phi_)
    // , S_pfield_01_(S_pfield_, 0, size_phi_, size_phi_, 2 * size_phi_)
    // , S_pfield_10_(S_pfield_, size_phi_, 2 * size_phi_, 0, size_phi_)
    // , S_pfield_11_(S_pfield_, size_phi_, 2 * size_phi_, size_phi_, 2 * size_phi_)
    // , S_pfield_12_(S_pfield_, size_phi_, 2 * size_phi_, 2 * size_phi_, 3 * size_phi_)
    // , S_pfield_20_(S_pfield_, 2 * size_phi_, 3 * size_phi_, 0, size_phi_)
    // , S_pfield_22_(S_pfield_, 2 * size_phi_, 3 * size_phi_, 2 * size_phi_, 3 * size_phi_)
    , M_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
    , D_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
    , M_ell_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
    , M_nonlin_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
    , G_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
    , pfield_solver_(std::make_shared<SolverType>())
    , pfield_mass_matrix_solver_(std::make_shared<DirectSolverType>())
    , pfield_phimu_matrixop_(M_pfield_,
                             *pfield_mass_matrix_solver_,
                             D_pfield_,
                             M_ell_pfield_,
                             G_pfield_,
                             M_nonlin_pfield_,
                             phi_dirichlet_constraints_,
                             gamma_,
                             epsilon_,
                             Be_,
                             Ca_)
    , pfield_jac_linear_op_(
          M_pfield_, D_pfield_, M_ell_pfield_, phi_dirichlet_constraints_, gamma_, epsilon_, Be_, this)
    , pfield_phimu_solver_(
          std::make_shared<PfieldPhiMuSolverType>(pfield_phimu_matrixop_, identity_prec_, 1e-10, 20, 10000, false))
    , D_pfield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, D_pfield_))
    , G_pfield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, G_pfield_))
    , M_nonlin_pfield_op_(
          std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, M_nonlin_pfield_))
    , pfield_rhs_vector_(3 * size_phi_, 0., num_mutexes_pfield_)
    , pfield_g_vector_(pfield_rhs_vector_, 0, size_phi_)
    , pfield_h_vector_(pfield_rhs_vector_, size_phi_, 2 * size_phi_)
    , pfield_f_vector_(pfield_rhs_vector_, 2 * size_phi_, 3 * size_phi_)
    , phi_tmp_eigen_(size_phi_, 0., 0)
    , phi_tmp_eigen2_(size_phi_, 0., 0)
    , pfield_tmp_phimu_(num_cells_, EigenVectorType(2 * size_phi_, 0., 0))
    , pfield_tmp_rhs_phimu_(2 * size_phi_, 0., 0)
    , pfield_tmp_rhs_phimu_phi_(pfield_tmp_rhs_phimu_, 0, size_phi_)
    , pfield_tmp_rhs_phimu_mu_(pfield_tmp_rhs_phimu_, size_phi_, 2 * size_phi_)
    // , pfield_old_result_(num_cells_, EigenVectorType(3 * size_phi_, 0.))
    , pfield_deim_input_dofs_(num_cells_)
    , phinat_deim_input_dofs_begin_(num_cells_)
    , mu_deim_input_dofs_begin_(num_cells_)
    , pfield_deim_output_dofs_(num_cells_)
    , pfield_deim_unique_output_dofs_(num_cells_)
    , phi_deim_output_dofs_(num_cells_)
    , phinat_deim_output_dofs_(num_cells_)
    , mu_deim_output_dofs_(num_cells_)
    , both_mu_and_phi_deim_output_dofs_(num_cells_)
    , pfield_deim_entities_(num_cells_)
    , pfield_tmp_vec_(3 * size_phi_, 0., 0)
    , pfield_tmp_vec2_(3 * size_phi_, 0., 0)
    , phi_tmp_vec_(size_phi_, 0., 0)
    , u_tmp_vec_(size_u_, 0., 0)
    , u_tmp_(u_space_)
    , u_tmp_local_(std::make_shared<PerThreadVectorLocalFunc>())
    , P_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
    , Pnat_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
    , phi_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
    , phinat_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
    , mu_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
  {

    /************************** create and project initial values*****************************************
     ************************** we only need initial values for P and phi ********************************
     ************************** mu_initial is only needed if linearization is used ***********************/

    std::shared_ptr<const XT::Functions::FunctionInterface<d, d>> u_initial_func;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> phi_initial_funcs;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> mu_initial_funcs;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d, d>>> P_initial_funcs;

    if (testcase == "single_cell") {
      // Initially, cell is circular with Radius R=5 and placed in the center of the domain
      // \Omega = [0, 160] \times [0, 40].
      // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
      // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
      FieldVector<double, d> center{upper_right_[0] / 2., upper_right_[1] / 2.};
      auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
      phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
          50,
          /*evaluate=*/
          [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2.) * epsilon)); },
          /*name=*/"phi_initial"));
      mu_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d>>(
          50,
          /*evaluate=*/
          [& phi_in = phi_initial_funcs[0], epsilon](const auto& x, const auto& param) {
            // TODO: add approximation of laplacian term
            const auto phi = phi_in->evaluate(x, param);
            return 1. / epsilon * (std::pow(phi, 3) - phi);
          },
          /*name=*/"mu_initial"));

      // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
      // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
      std::srand(1); // set seed for std::rand to 1
      P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
          50,
          /*evaluate=*/
          [& phi_in = phi_initial_funcs[0]](const auto& x, const auto& param) {
            // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
            // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
            // auto ret = FieldVector<double, d>({1. + rand1, 0. +
            // rand2});
            auto ret = FieldVector<double, d>({1., 0.});
            ret *= (phi_in->evaluate(x, param) + 1.) / 2.;
            return ret;
          },
          /*name=*/"P_initial"));
      u_initial_func = std::make_shared<const XT::Functions::ConstantFunction<d, d>>(0.);

      // interpolate initial and boundary values
    } else if (testcase == "two_cells") {
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
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    }

    /*************************************************************************************************
     ******************************* create variables, set initial values ****************************
     *************************************************************************************************/

    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    const XT::Functions::ConstantFunction<d> minus_one(-1.);
    XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
    default_interpolation(*u_initial_func, u_);
    // create system and temporary vectors, DiscreteFunctions, etc.
    for (size_t kk = 0; kk < num_cells_; kk++) {
      P_view_.emplace_back(ofield_vectors_[kk], 0, size_u_);
      Pnat_view_.emplace_back(ofield_vectors_[kk], size_u_, 2 * size_u_);
      phi_view_.emplace_back(pfield_vectors_[kk], 0, size_phi_);
      phinat_view_.emplace_back(pfield_vectors_[kk], size_phi_, 2 * size_phi_);
      mu_view_.emplace_back(pfield_vectors_[kk], 2 * size_phi_, 3 * size_phi_);
      const auto kk_str = XT::Common::to_string(kk);
      P_.emplace_back(make_discrete_function(u_space_, P_view_[kk], "P_" + kk_str));
      Pnat_.emplace_back(make_discrete_function(u_space_, Pnat_view_[kk], "Pnat_" + kk_str));
      phi_.emplace_back(make_discrete_function(phi_space_, phi_view_[kk], "phi_" + kk_str));
      phinat_.emplace_back(make_discrete_function(phi_space_, phinat_view_[kk], "phinat_" + kk_str));
      mu_.emplace_back(make_discrete_function(phi_space_, mu_view_[kk], "mu_" + kk_str));
      P_tmp_.emplace_back(u_space_);
      Pnat_tmp_.emplace_back(u_space_);
      phi_tmp_.emplace_back(phi_space_);
      phinat_tmp_.emplace_back(phi_space_);
      mu_tmp_.emplace_back(phi_space_);
      default_interpolation(*phi_initial_funcs[kk], phi_[kk]);
      boundary_interpolation(minus_one, phi_[kk], all_dirichlet_boundary_info, XT::Grid::DirichletBoundary{});
      default_interpolation(*mu_initial_funcs[kk], mu_[kk]);
      default_interpolation(*P_initial_funcs[kk], P_[kk]);
      pfield_tmp_phimu_phi_.emplace_back(pfield_tmp_phimu_[kk], 0, size_phi_);
      pfield_tmp_phimu_mu_.emplace_back(pfield_tmp_phimu_[kk], size_phi_, 2 * size_phi_);
    }

    /*************************************************************************************************
     *************************************** Stokes **************************************************
     *************************************************************************************************/
    if (Re_ > 1e-2)
      DUNE_THROW(NotImplemented, "No Navier-Stokes solver implemented yet!");

    MatrixOperator<MatrixViewType, PGV, 1, 1, d> B_stokes_op(grid_view_, p_space_, u_space_, B_stokes_);
    MatrixOperator<MatrixType, PGV, 1, 1, 1> M_p_stokes_op(grid_view_, p_space_, p_space_, M_p_stokes_);
    // calculate A_{ij} as \int \nabla v_i \nabla v_j
    // A_stokes_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalSymmetricEllipticIntegrand<E>(1.)));
    A_stokes_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>()));
    // calculate B_{ij} as \int \nabla p_i div(v_j)
    B_stokes_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
        LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
    M_p_stokes_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>()));

    auto p_basis_integrated_functional = make_vector_functional(p_space_, p_basis_integrated_vector_);
    const XT::Functions::ConstantGridFunction<E> one_function(1);
    p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), one_function)));
    B_stokes_op.append(p_basis_integrated_functional);

    // Dirichlet constrainst for u
    A_stokes_op_->append(u_dirichlet_constraints_);
    // assemble everything
    A_stokes_op_->append(B_stokes_op);
    A_stokes_op_->append(M_p_stokes_op);
    A_stokes_op_->assemble(use_tbb_);

    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the p_size_-th row of
    // [A B; B^T 0] to the unit vector.
    const size_t dof_index = 0;
    S_stokes_.set_entry(size_u_ + dof_index, size_u_ + dof_index, 1.);
    B_stokes_.clear_col(dof_index);
    stokes_g_vector_.set_entry(dof_index, 0.);

    u_dirichlet_constraints_.apply(A_stokes_, false, true);
    for (const auto& DoF : u_dirichlet_constraints_.dirichlet_DoFs())
      B_stokes_.clear_row(DoF);

    // Set B^T
    const auto B_pattern = B_stokes_.pattern();
    for (size_t ii = 0; ii < size_u_; ii++)
      for (const auto& jj : B_pattern.inner(ii))
        BT_stokes_.set_entry(jj, ii, B_stokes_.get_entry(ii, jj));

    S_colmajor_ = S_stokes_.backend();
    stokes_solver_->analyzePattern(S_colmajor_);
    stokes_solver_->factorize(S_colmajor_);

    /*************************************************************************************************
     ************************************ Orientationfield *******************************************
     *************************************************************************************************/
    // calculate M_{ij} as \int \psi_i phi_j
    M_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(1.)));
    C_ofield_elliptic_part_ *= 0.;
    MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, u_space_, u_space_, C_ofield_elliptic_part_);
    ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
    M_ofield_op_->append(ofield_elliptic_op);
    M_ofield_op_->assemble(use_tbb_);
    // S_schur_ofield_colmajor_ = S_schur_ofield_.backend();
    // ofield_solver_->analyzePattern(S_ofield_.backend());
    // ofield_schur_solver_->analyzePattern(S_schur_ofield_.backend());
    M_ofield_colmajor_ = M_ofield_.backend();
    ofield_mass_matrix_solver_->analyzePattern(M_ofield_colmajor_);
    ofield_mass_matrix_solver_->factorize(M_ofield_colmajor_);


    /*************************************************************************************************
     **************************************** Phasefield *********************************************
     *************************************************************************************************/

    MatrixOperator<MatrixType, PGV, 1> M_pfield_op(grid_view_, phi_space_, phi_space_, M_pfield_);
    M_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(1.)));
    MatrixOperator<MatrixType, PGV, 1> M_ell_pfield_op(grid_view_, phi_space_, phi_space_, M_ell_pfield_);
    M_ell_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
    M_pfield_op.append(phi_dirichlet_constraints_);
    M_pfield_op.append(M_ell_pfield_op);
    M_pfield_op.assemble(use_tbb_);

    M_pfield_colmajor_ = M_pfield_.backend();
    pfield_mass_matrix_solver_->analyzePattern(M_pfield_colmajor_);
    pfield_mass_matrix_solver_->factorize(M_pfield_colmajor_);
  } // constructor

  void pfield_jacobian()
  {
    // Set matrix S_{22} = C = M
    // S_pfield_22_ = M_pfield_;
    // // apply Dirichlet constraints to linear part
    // for (const auto& DoF : phi_dirichlet_constraints_.dirichlet_DoFs())
    //   S_pfield_22_.clear_row(DoF);
    // // Set matrix S_{11} = H = M
    // S_pfield_11_ = M_pfield_;
    // // Set matrix S_{01} = E
    // S_pfield_01_ = M_ell_pfield_;
    // S_pfield_01_ *= gamma_;
    // pfield_solver_->analyzePattern(S_pfield_.backend());
    // XT::LA::write_matrix_market(M_pfield_, "phasefield_mass_matrix.mtx");
  }

  size_t num_cells() const
  {
    return num_cells_;
  }

  bool linear() const
  {
    return linearize_;
  }

  bool finished() const
  {
    return XT::Common::FloatCmp::eq(t_end_, t_);
  }

  //******************************************************************************************************************
  //********************************* Solve methods for whole system of equations ************************************
  //******************************************************************************************************************

  // Solves whole system of equations using the values stored in stokes_vector_, ofield_vectors_ and pfield_vectors_  as
  // initial values Returns the whole trajectory, i.e., ret[i] contains the results in the i-th timestep. The first
  // num_cells_ entries of ret[i] correspond to the phasefield for each cell, the next num_cells entries are the
  // orientation field vectors, and the last one is the stokes vector. dt: Time step length. write: Whether to write
  // .vtu and .txt files. write_step: Time interval at which results should be written. If negative, all steps are
  // written. Ignored if write = false. filename: Prefix for .vtu and .txt files. Ignored if write = false. subsampling:
  // Whether to use subsampling for visualization. Ignored if write = false.
  std::vector<std::vector<VectorType>> solve(const double dt,
                                             const bool write,
                                             const double write_step,
                                             const std::string filename = "cellmodel",
                                             const bool subsampling = true)
  {
    std::vector<std::vector<VectorType>> ret(1 + 2 * num_cells_);
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      ret[kk].push_back(pfield_vectors_[kk]);
      ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
    }
    ret[2 * num_cells_].push_back(stokes_vector_);
    // implicit Euler timestepping
    if (write)
      visualize(filename, 0, t_, subsampling);
    assert(Dune::XT::Common::FloatCmp::ge(t_end_, t_));
    double next_save_time = t_ + write_step > t_end_ ? t_end_ : t_ + write_step;
    size_t save_step_counter = 1;

    while (Dune::XT::Common::FloatCmp::lt(t_, t_end_)) {
      double max_dt = dt;
      // match saving times and t_end_ exactly
      if (Dune::XT::Common::FloatCmp::gt(t_ + dt, t_end_))
        max_dt = t_end_ - t_;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t_ << std::endl;
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        prepare_pfield_op(dt, kk);
        ret[kk].push_back(apply_inverse_pfield_op(ret[kk].back(), kk));
        set_pfield_vec(kk, ret[kk].back());
        std::cout << "Pfield " << kk << " done" << std::endl;
        prepare_ofield_op(dt, kk);
        ret[num_cells_ + kk].push_back(apply_inverse_ofield_op(ret[num_cells_ + kk].back(), kk));
        set_ofield_vec(kk, ret[num_cells_ + kk].back());
        std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      prepare_stokes_op();
      ret[2 * num_cells_].push_back(apply_inverse_stokes_op());
      set_stokes_vec(ret[2 * num_cells_].back());
      std::cout << "Stokes done" << std::endl;

      t_ += actual_dt;

      // check if data should be written in this timestep (and write)
      if (write) {
        if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t_, next_save_time)) {
          visualize(filename, save_step_counter, t_, subsampling);
          next_save_time += write_step;
          ++save_step_counter;
        }
      }
    } // while (t_ < t_end_)
    return ret;
  }

  // Like solve, but only computes and returns the next n timesteps
  std::vector<std::vector<VectorType>> next_n_timesteps(const size_t n, const double dt)
  {
    std::vector<std::vector<VectorType>> ret(2 * num_cells_ + 1);
    size_t count = 0;
    if (XT::Common::is_zero(t_)) {
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        ret[kk].push_back(pfield_vectors_[kk]);
        ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
      }
      ret[2 * num_cells_].push_back(stokes_vector_);
      // Hack to avoid adding initial_values twice
      ++count;
      t_ = 1e-100;
    }
    // Undo hack to avoid adding initial_values twice
    if (XT::Common::is_zero(t_ - 1e-100) && count == 0)
      t_ = 0.;

    assert(Dune::XT::Common::FloatCmp::ge(t_end_, t_));

    // implicit Euler timestepping
    while (Dune::XT::Common::FloatCmp::lt(t_, t_end_) && count < n) {
      double max_dt = dt;
      // match saving times and t_end_ exactly
      if (Dune::XT::Common::FloatCmp::gt(t_ + dt, t_end_))
        max_dt = t_end_ - t_;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        prepare_pfield_op(dt, kk);
        pfield_vectors_[kk] = apply_inverse_pfield_op(pfield_vectors_[kk], kk);
        ret[kk].push_back(pfield_vectors_[kk]);
        // std::cout << "Pfield " << kk << " done" << std::endl;
        prepare_ofield_op(dt, kk);
        ofield_vectors_[kk] = apply_inverse_ofield_op(ofield_vectors_[kk], kk);
        ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
        // std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      prepare_stokes_op();
      stokes_vector_ = apply_inverse_stokes_op();
      ret[2 * num_cells_].push_back(stokes_vector_);
      // std::cout << "Stokes done" << std::endl;

      ++count;
      t_ += actual_dt;
    } // while (t_ < t_end_)
    return ret;
  }

  //******************************************************************************************************************
  //********************************* Product operators (mass matrix application) ************************************
  //******************************************************************************************************************

  // applies the pfield mass matrix to phi, phinat, mu
  // To calculate the sum of the squared L2 products of phi, phinat and mu, calculate the inner product of the result
  // with vec.
  VectorType apply_pfield_product_op(const VectorType& vec) const
  {
    VectorType ret(3 * size_phi_);
    ConstVectorViewType phi_view(vec, 0, size_phi_);
    ConstVectorViewType phinat_view(vec, size_phi_, 2 * size_phi_);
    ConstVectorViewType mu_view(vec, 2 * size_phi_, 3 * size_phi_);
    VectorViewType phi_ret_view(ret, 0, size_phi_);
    VectorViewType phinat_ret_view(ret, size_phi_, 2 * size_phi_);
    VectorViewType mu_ret_view(ret, 2 * size_phi_, 3 * size_phi_);
    M_pfield_.mv(phi_view, phi_ret_view);
    M_pfield_.mv(phinat_view, phinat_ret_view);
    M_pfield_.mv(mu_view, mu_ret_view);
    return ret;
  }

  // applies the ofield mass matrix to P, Pnat
  // To calculate the sum of the squared L2 products of P and Pnat, calculate the inner product of the result with vec.
  VectorType apply_ofield_product_op(const VectorType& vec) const
  {
    VectorType ret(2 * size_u_);
    ConstVectorViewType P_view(vec, 0, size_u_);
    ConstVectorViewType Pnat_view(vec, size_u_, 2 * size_u_);
    VectorViewType P_ret_view(ret, 0, size_u_);
    VectorViewType Pnat_ret_view(ret, size_u_, 2 * size_u_);
    M_ofield_.mv(P_view, P_ret_view);
    M_ofield_.mv(Pnat_view, Pnat_ret_view);
    return ret;
  }

  // applies the ofield mass matrix to P, Pnat
  // To calculate the sum of the squared L2 products of P and Pnat, calculate the inner product of the result with vec.
  VectorType apply_stokes_product_op(const VectorType& vec) const
  {
    VectorType ret(size_u_ + size_p_);
    ConstVectorViewType u_view(vec, 0, size_u_);
    ConstVectorViewType p_view(vec, size_u_, size_u_ + size_p_);
    VectorViewType u_ret_view(ret, 0, size_u_);
    VectorViewType p_ret_view(ret, size_u_, size_u_ + size_p_);
    // The Orientation field variables and u have the same basis so use M_ofield_
    M_ofield_.mv(u_view, u_ret_view);
    M_p_stokes_.mv(p_view, p_ret_view);
    return ret;
  }

  //******************************************************************************************************************
  //*****************************************  Visualization   *******************************************************
  //******************************************************************************************************************

  // Visualizes given vector as phasefield finite element vector
  void visualize_pfield(const std::string& filename, const VectorType& vec, const bool subsampling = true) const
  {
    auto vtk_writer = phi_[0].create_vtkwriter(phi_space_.grid_view(), subsampling);
    const ConstVectorViewType phi_vec(vec, 0, size_phi_);
    const ConstVectorViewType phinat_vec(vec, size_phi_, 2 * size_phi_);
    const ConstVectorViewType mu_vec(vec, 2 * size_phi_, 3 * size_phi_);
    const auto phi_func = make_discrete_function(phi_space_, phi_vec, "phi");
    const auto phinat_func = make_discrete_function(phi_space_, phinat_vec, "phinat");
    const auto mu_func = make_discrete_function(phi_space_, mu_vec, "mu");
    phi_func.add_to_vtkwriter(*vtk_writer);
    phinat_func.add_to_vtkwriter(*vtk_writer);
    mu_func.add_to_vtkwriter(*vtk_writer);
    phi_[0].write_visualization(*vtk_writer, filename);
  } // void visualize_pfield(...)

  // Visualizes given vector as orientation field finite element vector
  void visualize_ofield(const std::string& filename, const VectorType& vec, const bool subsampling = true) const
  {
    auto vtk_writer = P_[0].create_vtkwriter(u_space_.grid_view(), subsampling);
    const ConstVectorViewType P_vec(vec, 0, size_u_);
    const ConstVectorViewType Pnat_vec(vec, size_u_, 2 * size_u_);
    const auto P_func = make_discrete_function(u_space_, P_vec, "P");
    const auto Pnat_func = make_discrete_function(u_space_, Pnat_vec, "Pnat");
    P_func.add_to_vtkwriter(*vtk_writer);
    Pnat_func.add_to_vtkwriter(*vtk_writer);
    P_[0].write_visualization(*vtk_writer, filename);
  } // void visualize_ofield(...)

  // Visualizes given vector as stokes finite element vector
  void visualize_stokes(const std::string& filename, const VectorType& vec, const bool subsampling = true) const
  {
    auto vtk_writer = u_.create_vtkwriter(u_.space().grid_view(), subsampling);
    const ConstVectorViewType u_vec(vec, 0, size_u_);
    const ConstVectorViewType p_vec(vec, size_u_, size_u_ + size_p_);
    const auto u_func = make_discrete_function(u_space_, u_vec, "u");
    const auto p_func = make_discrete_function(p_space_, p_vec, "p");
    u_.add_to_vtkwriter(*vtk_writer);
    p_func.add_to_vtkwriter(*vtk_writer);
    u_.write_visualization(*vtk_writer, filename);
  } // void visualize_stokes(...)

  // Visualizes variables currently stored in this class.
  // If txt = true, also writes textfiles containing the values.
  void visualize(const std::string& prefix,
                 const size_t step,
                 const double t,
                 const bool subsampling = true,
                 const bool vtu = true,
                 const bool txt = false) const
  {
    auto vtk_writer = u_.create_vtkwriter(u_.space().grid_view(), subsampling);
    std::string postfix = "_" + XT::Common::to_string(step);
    if (vtu) {
      u_.add_to_vtkwriter(*vtk_writer);
      p_.add_to_vtkwriter(*vtk_writer);
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        // std::cout << "phi l2_norm: " << l2_norm(phi_[kk].space().grid_view(), phi_[kk]) << std::endl;
        // std::cout << "phinat l2_norm: " << l2_norm(phinat_[kk].space().grid_view(), phinat_[kk]) << std::endl;
        // std::cout << "mu l2_norm: " << l2_norm(mu_[kk].space().grid_view(), mu_[kk]) << std::endl;
        P_[kk].add_to_vtkwriter(*vtk_writer);
        Pnat_[kk].add_to_vtkwriter(*vtk_writer);
        phi_[kk].add_to_vtkwriter(*vtk_writer);
        phinat_[kk].add_to_vtkwriter(*vtk_writer);
        mu_[kk].add_to_vtkwriter(*vtk_writer);
        phi_[kk].add_gradient_to_vtkwriter(*vtk_writer);
        phinat_[kk].add_gradient_to_vtkwriter(*vtk_writer);
        mu_[kk].add_gradient_to_vtkwriter(*vtk_writer);
      }
      u_.write_visualization(*vtk_writer, prefix + postfix);
    } // if (vtu)
    if (txt) {
      write_to_textfile(u_, prefix, step, t);
      write_to_textfile(p_, prefix, step, t);
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        write_to_textfile(P_[kk], prefix, step, t);
        write_to_textfile(Pnat_[kk], prefix, step, t);
        write_to_textfile(phi_[kk], prefix, step, t);
        write_to_textfile(phinat_[kk], prefix, step, t);
        write_to_textfile(mu_[kk], prefix, step, t);
      } // kk
    } // if (txt)
  } // void visualize(...)

  //******************************************************************************************************************
  //*******************************  Methods to get and set variable values   ****************************************
  //******************************************************************************************************************

  // Sets stokes vector to stokes_vec
  void set_stokes_vec(const VectorType& stokes_vec)
  {
    DUNE_THROW_IF(
        stokes_vec.size() != size_u_ + size_p_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
    stokes_vector_ = stokes_vec;
  }

  // Sets orientation field vector belonging to cell to pfield_vec
  void set_ofield_vec(const size_t cell, const VectorType& ofield_vec)
  {
    DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
    DUNE_THROW_IF(ofield_vec.size() != 2 * size_u_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
    ofield_vectors_[cell] = ofield_vec;
  }

  // Sets phasefield vector belonging to cell to pfield_vec
  void set_pfield_vec(const size_t cell, const VectorType& pfield_vec)
  {
    DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
    DUNE_THROW_IF(
        pfield_vec.size() != 3 * size_phi_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
    pfield_vectors_[cell] = pfield_vec;
  }

  // Get stokes finite element vector
  const VectorType& stokes_vec()
  {
    return stokes_vector_;
  }

  // Get orientation field finite element vector belonging to cell
  const VectorType& ofield_vec(const size_t cell)
  {
    return ofield_vectors_[cell];
  }

  // Get phase field finite element vector belonging to cell
  const VectorType& pfield_vec(const size_t cell)
  {
    return pfield_vectors_[cell];
  }

  //******************************************************************************************************************
  //****** Prepare methods (calculate everything that is linear for the respective operator, but depends on **********
  //****** the values of other variables, so cannot be computed once and for all in the constructor )       **********
  //******************************************************************************************************************

  void prepare_stokes_op()
  {
    u_tmp_.dofs().vector() = u_.dofs().vector();
    for (size_t kk = 0; kk < num_cells_; kk++) {
      phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
      phinat_tmp_[kk].dofs().vector() = phinat_[kk].dofs().vector();
      P_tmp_[kk].dofs().vector() = P_[kk].dofs().vector();
      Pnat_tmp_[kk].dofs().vector() = Pnat_[kk].dofs().vector();
    }
    assemble_stokes_rhs();
  }

  void prepare_ofield_op(const double dt, const size_t cell, const bool restricted = false)
  {
    u_tmp_.dofs().vector() = u_.dofs().vector();
    P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
    phi_tmp_[cell].dofs().vector() = phi_[cell].dofs().vector();
    dt_ = dt;
    ofield_jac_linear_op_.prepare(dt, cell, restricted);
    assemble_ofield_rhs(dt, cell);
    assemble_ofield_linear_jacobian(dt, cell);
  }

  void prepare_pfield_op(const double dt, const size_t cell, const bool restricted = false)
  {
    u_tmp_.dofs().vector() = u_.dofs().vector();
    P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
    for (size_t kk = 0; kk < num_cells_; kk++) {
      phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
      mu_tmp_[kk].dofs().vector() = mu_[kk].dofs().vector();
    }
    assemble_pfield_rhs(dt, cell, restricted);
    assemble_pfield_linear_jacobian(dt, cell, restricted);
    pfield_jac_linear_op_.prepare(dt, cell, restricted);
    dt_ = dt;
  }

  void compute_restricted_ofield_dofs(const std::vector<size_t>& output_dofs, const size_t cell)
  {
    if (!ofield_deim_output_dofs_[cell] || *ofield_deim_output_dofs_[cell] != output_dofs) {
      const auto& pattern = create_ofield_pattern(size_u_, ofield_submatrix_pattern_);
      // We need to keep the original output_dofs which is unordered and may contain duplicates, as the restricted
      // operator will return exactly these dofs. For computations, however, we often need unique dofs.
      ofield_deim_output_dofs_[cell] = std::make_shared<std::vector<size_t>>(output_dofs);
      auto& unique_output_dofs = ofield_deim_unique_output_dofs_[cell];
      unique_output_dofs = output_dofs;
      std::sort(unique_output_dofs.begin(), unique_output_dofs.end());
      unique_output_dofs.erase(std::unique(unique_output_dofs.begin(), unique_output_dofs.end()),
                               unique_output_dofs.end());
      // sort output into dofs belonging to phi, phinat and mu
      auto& P_output_dofs = P_deim_output_dofs_[cell];
      auto& Pnat_output_dofs = Pnat_deim_output_dofs_[cell];
      P_output_dofs.clear();
      Pnat_output_dofs.clear();
      for (const auto& dof : unique_output_dofs) {
        if (dof < size_u_)
          P_output_dofs.push_back(dof);
        else
          Pnat_output_dofs.push_back(dof);
      }
      for (auto& dof : Pnat_output_dofs)
        dof -= size_u_;
      // get input dofs corresponding to output dofs
      auto& input_dofs = ofield_deim_input_dofs_[cell];
      input_dofs.clear();
      for (const auto& dof : unique_output_dofs) {
        const auto& new_input_dofs = pattern.inner(dof);
        input_dofs.insert(input_dofs.end(), new_input_dofs.begin(), new_input_dofs.end());
      }
      // sort and remove duplicate entries
      std::sort(input_dofs.begin(), input_dofs.end());
      input_dofs.erase(std::unique(input_dofs.begin(), input_dofs.end()), input_dofs.end());
      Pnat_deim_input_dofs_begin_[cell] =
          std::lower_bound(input_dofs.begin(), input_dofs.end(), size_u_) - input_dofs.begin();
      // store all entities that contain an output dof
      const auto& mapper = u_space_.mapper();
      DynamicVector<size_t> global_indices;
      ofield_deim_entities_[cell].clear();
      for (const auto& entity : Dune::elements(grid_view_)) {
        mapper.global_indices(entity, global_indices);
        maybe_add_entity(entity, global_indices, *ofield_deim_output_dofs_[cell], ofield_deim_entities_[cell]);
      } // entities
    } // if (not already computed)
  } // void compute_restricted_pfield_dofs(...)

  void compute_restricted_pfield_dofs(const std::vector<size_t>& output_dofs, const size_t cell)
  {
    if (!pfield_deim_output_dofs_[cell] || *pfield_deim_output_dofs_[cell] != output_dofs) {
      const auto& pattern = create_pfield_pattern(size_phi_, pfield_submatrix_pattern_);
      // We need to keep the original output_dofs which is unordered and may contain duplicates, as the restricted
      // operator will return exactly these dofs. For computations, however, we often need unique dofs.
      pfield_deim_output_dofs_[cell] = std::make_shared<std::vector<size_t>>(output_dofs);
      auto& unique_output_dofs = pfield_deim_unique_output_dofs_[cell];
      unique_output_dofs = output_dofs;
      std::sort(unique_output_dofs.begin(), unique_output_dofs.end());
      unique_output_dofs.erase(std::unique(unique_output_dofs.begin(), unique_output_dofs.end()),
                               unique_output_dofs.end());
      // sort output into dofs belonging to phi, phinat and mu
      auto& phi_output_dofs = phi_deim_output_dofs_[cell];
      auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
      auto& mu_output_dofs = mu_deim_output_dofs_[cell];
      auto& phinat_mu_output_dofs = both_mu_and_phi_deim_output_dofs_[cell];
      phi_output_dofs.clear();
      phinat_output_dofs.clear();
      mu_output_dofs.clear();
      phinat_mu_output_dofs.clear();
      for (const auto& dof : unique_output_dofs) {
        if (dof < size_phi_)
          phi_output_dofs.push_back(dof);
        else if (dof < 2 * size_phi_)
          phinat_output_dofs.push_back(dof);
        else
          mu_output_dofs.push_back(dof);
      }
      for (auto& dof : phinat_output_dofs)
        dof -= size_phi_;
      for (auto& dof : mu_output_dofs)
        dof -= 2 * size_phi_;
      for (const auto& dof : phinat_output_dofs)
        if (std::find(mu_output_dofs.begin(), mu_output_dofs.end(), dof) != mu_output_dofs.end())
          phinat_mu_output_dofs.push_back(dof);
      // get input dofs corresponding to output dofs
      auto& input_dofs = pfield_deim_input_dofs_[cell];
      input_dofs.clear();
      for (const auto& dof : unique_output_dofs) {
        const auto& new_input_dofs = pattern.inner(dof);
        input_dofs.insert(input_dofs.end(), new_input_dofs.begin(), new_input_dofs.end());
      }
      // sort and remove duplicate entries
      std::sort(input_dofs.begin(), input_dofs.end());
      input_dofs.erase(std::unique(input_dofs.begin(), input_dofs.end()), input_dofs.end());
      phinat_deim_input_dofs_begin_[cell] =
          std::lower_bound(input_dofs.begin(), input_dofs.end(), size_phi_) - input_dofs.begin();
      mu_deim_input_dofs_begin_[cell] =
          std::lower_bound(input_dofs.begin(), input_dofs.end(), 2 * size_phi_) - input_dofs.begin();
      // store all entities that contain an output dof
      const auto& mapper = phi_space_.mapper();
      DynamicVector<size_t> global_indices;
      pfield_deim_entities_[cell].clear();
      for (const auto& entity : Dune::elements(grid_view_)) {
        mapper.global_indices(entity, global_indices);
        maybe_add_entity(entity, global_indices, *pfield_deim_output_dofs_[cell], pfield_deim_entities_[cell]);
      } // entities
    } // if (not already computed)
  } // void compute_restricted_pfield_dofs(...)

  //******************************************************************************************************************
  //*********************************************** Apply operators **************************************************
  //******************************************************************************************************************

  // Applies stokes operator (applies the F if Stokes equation is F(y) = 0)
  VectorType apply_stokes_op(VectorType y) const
  {
    VectorType ret(size_u_ + size_p_, 0.);
    S_stokes_.mv(y, ret);
    ret -= stokes_rhs_vector_;
    return ret;
  }

  VectorType apply_stokes_op_with_param(VectorType y, const XT::Common::Parameter& /*param*/) const
  {
    // TODO: implement parameter handling
    return apply_stokes_op(y);
  }

  // Applies cell-th orientation field operator (applies F if the orientation field equation is F(y) = 0)
  VectorType apply_ofield_op(const VectorType& y, const size_t cell, const bool restricted = false) const
  {
    // linear part
    VectorType ret(y.size());
    ofield_jac_linear_op_.apply(y, ret);
    // subtract rhs
    ret -= ofield_rhs_vector_;
    if (!linearize_) {
      fill_tmp_ofield(cell, y, restricted);
      // nonlinear part
      VectorViewType res1_vec(ret, size_u_, 2 * size_u_);
      auto nonlinear_res_functional = make_vector_functional(u_space_, res1_vec);
      XT::Functions::GenericGridFunction<E, d, 1> nonlinear_res_pf(
          /*order = */ 3 * u_space_.max_polorder(),
          /*post_bind_func*/
          [cell, this](const E& element) { this->bind_P(cell, element); },
          /*evaluate_func*/
          [cell, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
            // evaluate P, divP
            auto ret = this->eval_P(cell, x_local, param);
            ret *= factor * (ret * ret);
            return ret;
          });
      nonlinear_res_functional.append(LocalElementIntegralFunctional<E, d>(
          local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, d>(), nonlinear_res_pf)));
      if (!restricted)
        nonlinear_res_functional.assemble(use_tbb_);
      else
        nonlinear_res_functional.assemble_range(ofield_deim_entities_[cell]);
    }
    return ret;
  }

  VectorType apply_ofield_op_with_param(const VectorType& y, const size_t cell, const double Pa)
  {
    update_ofield_parameters(cell, Pa);
    return apply_ofield_op(y, cell);
  }

  void update_ofield_parameters(const size_t cell, const double Pa)
  {
    // Pa may have been set to a new value already (via update_pfield_parameters)
    if (XT::Common::FloatCmp::ne(Pa, last_ofield_Pa_)) {
      std::cout << "Ofield params updated, old Pa = " << last_ofield_Pa_ << ", new Pa = " << Pa << std::endl;
      Pa_ = Pa;
      last_ofield_Pa_ = Pa_;
      C_ofield_elliptic_part_ *= 0.;
      MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, u_space_, u_space_, C_ofield_elliptic_part_);
      ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
      ofield_elliptic_op.assemble(use_tbb_);
      prepare_ofield_op(dt_, cell);
    }
  }

  // Applies cell-th phase field operator (applies F if phase field equation is F(y) = 0)
  VectorType apply_pfield_op(const VectorType& y, const size_t cell, const bool restricted)
  {
    const auto& output_dofs = *pfield_deim_output_dofs_[cell];
    const auto& unique_output_dofs = pfield_deim_unique_output_dofs_[cell];
    const auto& input_dofs = pfield_deim_input_dofs_[cell];
    auto& source = pfield_tmp_vec_;
    auto& residual = pfield_tmp_vec2_;
    // copy values to high-dimensional vector
    if (restricted)
      copy_ld_to_hd_vec(input_dofs, y, source);
    else
      source = y;
    // linear part
    pfield_jac_linear_op_.apply(source, residual);
    // subtract rhs
    const auto sub = sub_func<VectorType>(restricted);
    sub(residual, pfield_rhs_vector_, unique_output_dofs);
    // nonlinear part
    fill_tmp_pfield(cell, source, restricted);
    assemble_nonlinear_part_of_pfield_residual(residual, cell, restricted);
    if (restricted) {
      VectorType ret(output_dofs.size());
      for (size_t ii = 0; ii < output_dofs.size(); ++ii)
        ret[ii] = residual[output_dofs[ii]];
      return ret;
    } else {
      return residual;
    }
  }

  void
  update_pfield_parameters(const size_t cell, const bool restricted, const double Be, const double Ca, const double Pa)
  {
    if (XT::Common::FloatCmp::ne(Be, Be_) || XT::Common::FloatCmp::ne(Ca, Ca_)
        || XT::Common::FloatCmp::ne(Pa, last_pfield_Pa_)) {
      std::cout << "Pfield params updated, old (Be, Ca, Pa) = " << Be_ << ", " << Ca_ << ", " << last_pfield_Pa_
                << ", new = " << Be << ", " << Ca << ", " << Pa << std::endl;
      Be_ = Be;
      Ca_ = Ca;
      Pa_ = Pa;
      last_pfield_Pa_ = Pa_;
      // TODO: we do not need to reassemble rhs if only Ca or gamma is changed and we do not need to prepare the phimu
      // operator again if only c_1 or Pa have changed
      assemble_pfield_rhs(dt_, cell, restricted);
      pfield_jac_linear_op_.set_params(gamma_, epsilon_, Be);
      pfield_phimu_matrixop_.set_params(gamma_, epsilon_, Be, Ca);
      pfield_phimu_matrixop_.prepare(dt_);
    }
  }

  VectorType apply_pfield_op_with_param(const VectorType& y,
                                        const size_t cell,
                                        const double Be,
                                        const double Ca,
                                        const double Pa,
                                        const bool restricted = false)
  {
    // std::cout << "Pfield: Be " << Be << ", Ca: " << Ca << ", Pa: " << Pa << std::endl;
    update_pfield_parameters(cell, restricted, Be, Ca, Pa);
    return apply_pfield_op(y, cell, restricted);
  }

  //******************************************************************************************************************
  //******************************************* Apply inverse operators **********************************************
  //******************************************************************************************************************

  // Applies inverse stokes operator (solves F(y) = 0)
  VectorType apply_inverse_stokes_op() const
  {
    // now solve the system
    // auto begin = std::chrono::steady_clock::now();
    EigenVectorType ret(size_u_ + size_p_);
    ret.backend() = stokes_solver_->solve(stokes_rhs_vector_.backend());
    // std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    // std::cout << "Solving Stokes took: " << time.count() << " s!" << std::endl;

    // ensure int_\Omega p = 0 (TODO: remove (?), not necessary as p is not used anywhere)
    // auto p_integral = p_basis_integrated_vector_ * p_.dofs().vector();
    // auto p_correction = make_discrete_function<VectorType>(p_space_, "p_corr");
    // XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain_);
    // default_interpolation(const_p_integral_func, p_correction);
    // p_ -= p_correction;

    return XT::Common::convert_to<VectorType>(ret);
  }

  // Applies inverse orientation field operator (solves F(y) = 0)
  // y_guess is the initial guess for the Newton iteration
  VectorType apply_inverse_ofield_op(const VectorType& y_guess, const size_t cell)
  {
    if (linearize_) {
      return solve_ofield_linear_system(ofield_rhs_vector_, cell);
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 200;
      const auto max_dampening_iter = 1000;

      auto l2_norm_P = l2_norm(grid_view_, P_[cell]);
      auto l2_norm_Pnat = l2_norm(grid_view_, Pnat_[cell]);

      // ********* compute residual *********
      auto begin = std::chrono::steady_clock::now();
      auto residual = apply_ofield_op(y_guess, cell);
      auto res_norm = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
      std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
      // std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

      size_t iter = 0;
      VectorType y_n = y_guess;
      VectorType y_n_plus_1 = y_guess;
      VectorType update;
      while (true) {
        if (res_norm < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_ofield_nonlinear_jacobian(y_n, cell);
        time = std::chrono::steady_clock::now() - begin;
        // std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        residual *= -1.;
        update = solve_ofield_linear_system(residual, cell);

        DUNE_THROW_IF(iter >= max_iter,
                      Exceptions::operator_error,
                      "max iterations in ofield reached!\n|residual|_l2 = " << res_norm << ", param: "
                                                                            << "(" << Re_ << ", " << 1. / Fa_inv_
                                                                            << ", " << xi_ << ")" << std::endl);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        // revert_ofield_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res_norm) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res_norm << "\nl = " << iter << "\n");
          y_n_plus_1 = y_n + update * lambda;
          residual = apply_ofield_op(y_n_plus_1, cell);
          candidate_res = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
          // std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        y_n = y_n_plus_1;
        res_norm = candidate_res;
        // std::cout << "Current res: " << candidate_res << std::endl;
        iter += 1;
      } // while (true)
      return y_n;
    }
  }

  VectorType apply_inverse_ofield_op_with_param(const VectorType& y_guess, const size_t cell, const double Pa)
  {
    // std::cout << "Ofield inverse params: Pa: " << Pa << std::endl;
    update_ofield_parameters(cell, Pa);
    return apply_inverse_ofield_op(y_guess, cell);
  }

  // Applies inverse phase field operator (solves F(y) = 0)
  // y_guess is the initial guess for the Newton iteration
  VectorType apply_inverse_pfield_op(const VectorType& y_guess, const size_t cell)
  {
    if (linearize_) {
      return solve_pfield_linear_system(pfield_rhs_vector_, cell);
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 200;
      const auto max_dampening_iter = 1000;

      const auto l2_norm_phi = l2_norm(grid_view_, phi_[cell]);
      const auto l2_norm_phinat = l2_norm(grid_view_, phinat_[cell]);
      const auto l2_norm_mu = l2_norm(grid_view_, mu_[cell]);

      // ********* compute residual *********
      auto begin = std::chrono::steady_clock::now();
      auto residual = apply_pfield_op(y_guess, cell, false);
      auto res_norm = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
      std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
      // std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

      size_t iter = 0;
      VectorType x_n = y_guess;
      VectorType x_n_plus_1 = y_guess;
      VectorType update;
      while (true) {
        if (res_norm < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_pfield_nonlinear_jacobian(x_n, cell, false);
        time = std::chrono::steady_clock::now() - begin;
        // std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        residual *= -1.;
        update = solve_pfield_linear_system(residual, cell);

        DUNE_THROW_IF(iter >= max_iter,
                      Exceptions::operator_error,
                      "max iterations in pfield reached!\n|residual|_l2 = " << res_norm << ", param: "
                                                                            << "(" << Re_ << ", " << 1. / Fa_inv_
                                                                            << ", " << xi_ << ")" << std::endl);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        // revert_pfield_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res_norm) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res_norm << "\nl = " << iter << "\n");
          x_n_plus_1 = x_n + update * lambda;
          residual = apply_pfield_op(x_n_plus_1, cell, false);
          candidate_res = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
          // std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        x_n = x_n_plus_1;
        res_norm = candidate_res;
        // std::cout << "Current res: " << candidate_res << std::endl;
        iter += 1;
      } // while (true)
      return x_n;
    }
  }

  VectorType apply_inverse_pfield_op_with_param(
      const VectorType& y_guess, const size_t cell, const double Be, const double Ca, const double Pa)
  {
    // std::cout << "Pfield inverse params: Be " << Be << ", Ca: " << Ca << ", Pa: " << Pa << std::endl;
    update_pfield_parameters(cell, false, Be, Ca, Pa);
    return apply_inverse_pfield_op(y_guess, cell);
  }

  //******************************************************************************************************************
  //********************************************** Apply jacobians ***************************************************
  //******************************************************************************************************************

  // Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
  // As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
  void
  apply_pfield_jacobian(const VectorType& source, VectorType& range, const size_t cell, const bool restricted = false)
  {
    const auto& output_dofs = *pfield_deim_output_dofs_[cell];
    const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
    const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
    VectorType& full_range = restricted ? pfield_tmp_vec_ : range;
    VectorViewType range_phi(full_range, 0, size_phi_);
    VectorViewType range_phinat(full_range, size_phi_, 2 * size_phi_);
    VectorViewType range_mu(full_range, 2 * size_phi_, 3 * size_phi_);
    const ConstVectorViewType source_phi(source, 0, size_phi_);
    const ConstVectorViewType source_phinat(source, size_phi_, 2 * size_phi_);
    const ConstVectorViewType source_mu(source, 2 * size_phi_, 3 * size_phi_);
    // linear part
    pfield_jac_linear_op_.apply(source, full_range);
    // nonlinear_part
    fill_tmp_pfield(cell, source, restricted);
    assemble_M_nonlin_pfield(cell, restricted);
    auto& tmp_vec = phi_tmp_vec_;
    const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
    const auto axpy = axpy_func<VectorViewType, VectorType>(restricted);
    const auto scal = scal_func<VectorType>(restricted);
    const auto add = add_func<VectorViewType, VectorType>(restricted);
    // apply missing parts of J (including the linear 1./Ca_ part)
    mv(M_nonlin_pfield_, source_mu, tmp_vec, phinat_output_dofs);
    axpy(range_phinat, 1. / (Be_ * std::pow(epsilon_, 2)), tmp_vec, phinat_output_dofs);
    mv(M_pfield_, source_mu, tmp_vec, phinat_output_dofs);
    axpy(range_phinat, 1. / Ca_, tmp_vec, phinat_output_dofs);
    // apply missing parts of A
    mv(M_nonlin_pfield_, source_phi, tmp_vec, mu_output_dofs);
    scal(tmp_vec, 1. / epsilon_, mu_output_dofs);
    // TODO: Only apply if DoF is both in dirichlet_dofs and in mu_output_dofs in restricted case
    for (const auto& DoF : phi_dirichlet_constraints_.dirichlet_DoFs())
      tmp_vec[DoF] = source_phi[DoF];
    add(range_mu, tmp_vec, mu_output_dofs);
    // apply G
    assemble_G_pfield(cell, restricted);
    mv(G_pfield_, source_phi, tmp_vec, phinat_output_dofs);
    add(range_phinat, tmp_vec, phinat_output_dofs);

    if (restricted)
      for (size_t ii = 0; ii < output_dofs.size(); ++ii)
        range.set_entry(ii, full_range.get_entry(output_dofs[ii]));
  }

  void
  apply_inverse_pfield_jacobian(const VectorType& source, const VectorType& rhs, VectorType& range, const size_t cell)
  {
    assemble_pfield_nonlinear_jacobian(source, cell, false);
    range = solve_pfield_linear_system(rhs, cell);
  }

  // Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
  // As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
  void
  apply_ofield_jacobian(const VectorType& source, VectorType& range, const size_t cell, const bool restricted = false)
  {
    const auto& output_dofs = *ofield_deim_output_dofs_[cell];
    const auto& Pnat_output_dofs = Pnat_deim_output_dofs_[cell];
    VectorType& full_range = restricted ? ofield_tmp_vec_ : range;
    VectorViewType range_P(full_range, 0, size_u_);
    VectorViewType range_Pnat(full_range, size_u_, 2 * size_u_);
    const ConstVectorViewType source_P(source, 0, size_u_);
    const ConstVectorViewType source_Pnat(source, size_u_, 2 * size_u_);
    // linear part
    ofield_jac_linear_op_.apply(source, full_range);
    // nonlinear_part
    fill_tmp_ofield(cell, source, restricted);
    assemble_C_ofield_nonlinear_part(cell, restricted);
    auto& tmp_vec = u_tmp_vec_;
    const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
    const auto add = add_func<VectorViewType, VectorType>(restricted);
    mv(C_ofield_nonlinear_part_, source_P, tmp_vec, Pnat_output_dofs);
    add(range_Pnat, tmp_vec, Pnat_output_dofs);
    if (restricted)
      for (size_t ii = 0; ii < output_dofs.size(); ++ii)
        range.set_entry(ii, full_range.get_entry(output_dofs[ii]));
  }

  void
  apply_inverse_ofield_jacobian(const VectorType& source, const VectorType& rhs, VectorType& range, const size_t cell)
  {
    assemble_ofield_nonlinear_jacobian(source, cell, false);
    range = solve_ofield_linear_system(rhs, cell);
  }

  void apply_stokes_jacobian(const VectorType& source, VectorType& range, const bool /*restricted*/ = false)
  {
    S_stokes_.mv(source, range);
  }

  void apply_inverse_stokes_jacobian(const VectorType& rhs, VectorType& range, const size_t cell)
  {
    EigenVectorType rhs_eigen = XT::Common::convert_to<EigenVectorType>(rhs);
    EigenVectorType ret(stokes_solver_->solve(rhs_eigen.backend()));
    range = XT::Common::convert_to<VectorType>(ret);
  }

  //******************************************************************************************************************
  //**************************** Methods to assemble rhs, residuals and jacobians ************************************
  //******************************************************************************************************************

  // Computes stokes rhs using currently stored values of variables and stores in stokes_rhs_vector_
  void assemble_stokes_rhs()
  {
    auto f_functional = make_vector_functional(u_space_, stokes_f_vector_);
    // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
    f_functional.append(LocalElementIntegralFunctional<E, d>(
        /*order*/ [& u_space =
                       u_space_](const auto& test_basis,
                                 const auto& param) { return 3 * u_space.max_polorder() + test_basis.order(param); },
        /*evaluate_func*/
        [this](const auto& test_basis,
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
          for (size_t kk = 0; kk < this->num_cells_; ++kk) {
            const auto P = this->eval_P(kk, x_local, param);
            const auto Pnat = this->eval_Pnat(kk, x_local, param);
            const auto phi = this->eval_phi(kk, x_local, param);
            const auto phinat = this->eval_phinat(kk, x_local, param);
            const auto grad_P = this->grad_P(kk, x_local, param);
            const auto grad_phi = this->grad_phi(kk, x_local, param);
            const auto phi_tilde = (phi + 1.) / 2.;

            // evaluate rhs terms
            const auto phinat_grad_phi = grad_phi * phinat;
            auto grad_P_T_times_Pnat = P;
            grad_P.mtv(Pnat, grad_P_T_times_Pnat);
            for (size_t ii = 0; ii < sz; ++ii) {
              for (size_t mm = 0; mm < d; ++mm) {
                result[ii] += (phinat_grad_phi[mm] + grad_P_T_times_Pnat[mm]) * test_basis_values_[ii][mm];
                for (size_t nn = 0; nn < d; ++nn)
                  result[ii] += (-this->Fa_inv_ * phi_tilde * P[mm] * P[nn] - 0.5 * (this->xi_ + 1) * Pnat[mm] * P[nn]
                                 - 0.5 * (this->xi_ - 1) * P[mm] * Pnat[nn])
                                * test_basis_grads_[ii][mm][nn];
              } // mm
            } // ii
          } // kk
        },
        /*post_bind_func*/
        [this](const E& element) {
          for (size_t kk = 0; kk < this->num_cells_; ++kk) {
            this->bind_phi(kk, element);
            this->bind_phinat(kk, element);
            this->bind_P(kk, element);
            this->bind_Pnat(kk, element);
          } // kk
        }));
    stokes_f_vector_ *= 0.;
    f_functional.assemble(use_tbb_);
    // apply dirichlet constraints for u.
    u_dirichlet_constraints_.apply(stokes_f_vector_);
  }

  // Computes orientation field rhs using currently stored values of variables and stores in ofield_rhs_vector_
  void assemble_ofield_rhs(const double /*dt*/, const size_t cell)
  {
    auto g_functional = make_vector_functional(u_space_, ofield_g_vector_);
    ofield_g_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, d> g(
        /*order = */ 3 * u_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) {
          this->bind_phi(cell, element);
          if (this->linearize_) {
            this->bind_P(cell, element);
          }
        },
        /*evaluate_func*/
        [cell, factor1 = beta_ / Pa_, factor2 = -2. * c_1_ / Pa_, this](const DomainType& x_local,
                                                                        const XT::Common::Parameter& param) {
          // evaluate rhs terms
          const auto grad_phi = this->grad_phi(cell, x_local, param);
          auto ret = grad_phi;
          ret *= factor1;
          if (linearize_) {
            const auto P_n = this->eval_P(cell, x_local, param);
            auto ret2 = P_n;
            ret2 *= factor2 * (P_n * P_n);
            ret += ret2;
          }
          return ret;
        });
    g_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, d>(), g)));
    g_functional.assemble(use_tbb_);
    M_ofield_.mv(P_[cell].dofs().vector(), ofield_f_vector_);
  }

  // Computes phase field rhs using currently stored values of variables and stores in pfield_rhs_vector_
  void assemble_pfield_rhs(const double /*dt*/, const size_t cell, const bool restricted)
  {
    const auto& phi_output_dofs = phi_deim_output_dofs_[cell];
    const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
    auto f_functional = make_vector_functional(phi_space_, pfield_f_vector_);
    auto h_functional = make_vector_functional(phi_space_, pfield_h_vector_);
    // calculate f
    if (linearize_)
      pfield_f_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) { this->bind_phi(cell, element); },
        /*evaluate_func*/
        [cell, two_epsilon_inv = 2. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate phi_
          const R phi_n = this->eval_phi(cell, x_local, param);
          return two_epsilon_inv * std::pow(phi_n, 3);
        });
    f_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), f_pf)));
    if (linearize_)
      h_functional.append(f_functional);

    // calculate g
    const auto mv = mv_func<VectorViewType>(restricted);
    mv(M_pfield_, phi_[cell].dofs().vector(), pfield_g_vector_, phi_output_dofs);

    // calculate h
    const auto scal = scal_func<VectorViewType>(restricted);
    scal(pfield_h_vector_, 0., phinat_output_dofs);
    XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) {
          this->bind_P(cell, element);
          if (this->linearize_) {
            this->bind_phi(cell, element);
            this->bind_mu(cell, element);
          }
        },
        /*evaluate_func*/
        [cell,
         factor0 = 6. / (Be_ * std::pow(epsilon_, 2)),
         factor1 = -c_1_ / (2. * Pa_),
         factor2 = -beta_ / Pa_,
         this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto Pn = this->eval_P(cell, x_local, param);
          const auto grad_P = this->grad_P(cell, x_local, param);
          R div_P(0.);
          for (size_t ii = 0; ii < d; ++ii)
            div_P += grad_P[ii][ii];
          auto ret = factor1 * (Pn * Pn) + factor2 * div_P;
          if (this->linearize_) {
            const auto phi_n = this->eval_phi(cell, x_local, param);
            const auto mu_n = this->eval_mu(cell, x_local, param);
            ret += factor0 * std::pow(phi_n, 2) * mu_n;
          }
          return ret;
        });
    h_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(1.), h_pf)));
    // assemble rhs
    if (!restricted)
      h_functional.assemble(use_tbb_);
    else
      h_functional.assemble_range(pfield_deim_entities_[cell]);
  }

  // assembles linear part of orientation field jacobian and stores in S_ofield_
  void assemble_ofield_linear_jacobian(const double dt, const size_t cell)
  {
    // calculate A
    // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
    XT::Functions::GenericGridFunction<E, d, d> Omega_minus_xi_D_transposed(
        /*order = */ std::max(u_space_.max_polorder() - 1, 0),
        /*post_bind_func*/
        [this](const E& element) { this->bind_u(element); },
        /*evaluate_func*/
        [this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate \nabla u
          auto grad_u = this->grad_u(x_local, param);
          auto grad_u_T = grad_u;
          grad_u_T.transpose();
          auto& ret = grad_u;
          ret *= (1. - this->xi_) / 2.;
          grad_u_T *= (1. + this->xi_) / 2.;
          ret -= grad_u_T;
          return ret;
        });
    A_ofield_.set_to_zero();
    A_ofield_op_->clear();
    A_ofield_op_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(Omega_minus_xi_D_transposed)));
    A_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(u_)));
    A_ofield_op_->assemble(use_tbb_);

    // calculate linear part S_10 = C
    C_ofield_linear_part_op_->clear();
    C_ofield_linear_part_ = C_ofield_elliptic_part_;
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_inv_phi(
        /*order = */ u_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) { this->bind_phi(cell, element); },
        /*evaluate_func*/
        [cell, factor = c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi = this->eval_phi(cell, x_local, param);
          return factor * phi;
        });
    C_ofield_linear_part_op_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(c1_Pa_inv_phi)));
    C_ofield_linear_part_op_->assemble(use_tbb_);
    S_schur_ofield_linear_part_.backend() = M_ofield_.backend();
    S_schur_ofield_linear_part_.axpy(dt, A_ofield_);
    S_schur_ofield_linear_part_.axpy(-dt / kappa_, C_ofield_linear_part_);

    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_ofield_nonlinear_jacobian(ofield_vec(cell), cell);
  }

  // assembles nonlinear part of orientation field jacobian and adds to S_ofield_
  // if assemble_ofield_linear_jacobian has been called first, S_ofield now contains the whole orientation field
  // jacobian
  void assemble_ofield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted = false) const
  {
    fill_tmp_ofield(cell, y, restricted);
    assemble_C_ofield_nonlinear_part(cell, restricted);
    S_schur_ofield_.backend() = S_schur_ofield_linear_part_.backend();
    S_schur_ofield_.axpy(-dt_ / kappa_, C_ofield_nonlinear_part_);
  }

  void assemble_C_ofield_nonlinear_part(const size_t cell, const bool restricted = false) const
  {
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_P2(
        /*order = */ 2. * u_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) { this->bind_P(cell, element); },
        /*evaluate_func*/
        [cell, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = this->eval_P(cell, x_local, param);
          return factor * P_n.two_norm2();
        });
    C_ofield_nonlinear_part_.set_to_zero();
    C_ofield_nonlinear_part_op_->clear();
    C_ofield_nonlinear_part_op_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(c1_Pa_P2)));
    C_ofield_nonlinear_part_op_->append(LocalElementIntegralBilinearForm<E, d>(
        LocalElementOtimesMatrixIntegrand<E, d>(P_tmp_[cell], -2. * c_1_ / Pa_)));
    C_ofield_nonlinear_part_op_->assemble(use_tbb_);
  }

  // reverts S_ofield_ to the state directly after assemble_ofield_linear_jacobian has been called the last time
  // Assumes only blocks that are touched by assemble_ofield_nonlinear_jacobian are modified
  // void revert_ofield_jacobian_to_linear() const
  // {
  //   S_ofield_10_ = C_ofield_linear_part_;
  //   S_schur_ofield_ = S_schur_ofield_linear_part_;
  // }

  // assembles linear part of phase field jacobian
  void assemble_pfield_linear_jacobian(const double /*dt*/, const size_t cell, const bool restricted)
  {
    // assemble matrix S_{00} = M + dt D
    assemble_D_pfield(cell, restricted);
    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_pfield_nonlinear_jacobian(pfield_vec(cell), cell, restricted);
  }

  // void assemble_S_pfield()
  // {
  //   assemble_D_pfield(cell, false);
  //   S_pfield_00_ = M_pfield_;
  //   S_pfield_00_.axpy(dt, D_pfield_);
  //   S_pfield_01_ = M_ell_pfield_;
  //   S_pfield_01_ *= gamma_ * dt;
  //   // linear part of matrix S_{12} = J
  //   S_pfield_12_ = M_ell_pfield_;
  //   S_pfield_12_ *= 1./Be_;
  //   // linear part of matrix S_{20} = A
  //   S_pfield_20_ = M_ell_pfield_;
  //   S_pfield_20_ *= epsilon_;
  // }

  void set_mat_to_zero(MatrixType& mat,
                       const bool restricted,
                       const XT::LA::SparsityPatternDefault& pattern,
                       const std::vector<size_t>& rows)
  {
    if (!restricted) {
      mat.set_to_zero();
    } else {
      for (const auto& row : rows)
        for (const auto& col : pattern.inner(row))
          mat.set_entry(row, col, 0.);
    }
  }

  void assemble_D_pfield(const size_t cell, const bool restricted)
  {
    const auto& phi_output_dofs = phi_deim_output_dofs_[cell];
    set_mat_to_zero(D_pfield_, restricted, pfield_submatrix_pattern_, phi_output_dofs);
    XT::Functions::GenericGridFunction<E, d, 1> minus_u(
        /*order = */ u_space_.max_polorder(),
        /*post_bind_func*/
        [this](const E& element) { this->bind_u(element); },
        /*evaluate_func*/
        [this](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = this->eval_u(x_local, param);
          ret *= -1;
          return ret;
        });
    D_pfield_op_->clear();
    D_pfield_op_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u)));
    if (!restricted)
      D_pfield_op_->assemble(use_tbb_);
    else
      D_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
  }

  // assembles nonlinear part of phase field jacobian
  void assemble_pfield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
  {
    fill_tmp_pfield(cell, y, restricted);
    assemble_M_nonlin_pfield(cell, restricted);
    assemble_G_pfield(cell, restricted);
    pfield_phimu_matrixop_.prepare(dt_);
  }

  // reverts S_pfield_ to the state directly after assemble_pfield_linear_jacobian has been called the last time
  // Assumes only blocks that are touched by assemble_pfield_nonlinear_jacobian are modified
  // void revert_pfield_jacobian_to_linear()
  // {
  //   // clear S_{10} = G
  //   S_pfield_10_ *= 0.;
  //   // linear part of matrix S_{12} = J
  //   S_pfield_12_ = M_ell_pfield_;
  //   S_pfield_12_ *= 1./Be_;
  //   // linear part of matrix S_{20} = A
  //   S_pfield_20_ = M_ell_pfield_;
  //   S_pfield_20_ *= epsilon_;
  // }

  // void assemble_pfield_restricted_linear_jacobian(const double dt, const size_t cell)
  // {
  //   // assemble matrix S_{00} = M/dt + D
  //   const auto& phi_output_dofs = phi_deim_output_dofs_[cell];
  //   const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  //   const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
  //   partially_copy_mat(phi_output_dofs, S_pfield_00_, M_pfield_);
  //   partially_scal_mat(phi_output_dofs, S_pfield_00_, 1. / dt);
  //   assemble_D_pfield(dt, cell, true);
  //   // linear part of matrix S_{12} = J
  //   partially_copy_mat(phinat_output_dofs, S_pfield_12_, M_ell_pfield_);
  //   partially_scal_mat(phinat_output_dofs, S_pfield_12_, 1./Be_);
  //   // linear part of matrix S_{20} = A
  //   partially_copy_mat(mu_output_dofs, S_pfield_20_, M_ell_pfield_);
  //   partially_scal_mat(mu_output_dofs, S_pfield_20_, epsilon_);
  // }

  // void revert_restricted_pfield_jacobian_to_linear(const size_t cell)
  // {
  //   const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  //   const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
  //   // clear S_{10} = G
  //   partially_scal_mat(phinat_output_dofs, S_pfield_10_, 0.);
  //   // linear part of matrix S_{12} = J
  //   partially_copy_mat(phinat_output_dofs, S_pfield_12_, M_ell_pfield_);
  //   partially_scal_mat(phinat_output_dofs, S_pfield_12_, 1./Be_);
  //   // linear part of matrix S_{20} = A
  //   partially_copy_mat(mu_output_dofs, S_pfield_20_, M_ell_pfield_);
  //   partially_scal_mat(mu_output_dofs, S_pfield_20_, epsilon_);
  // }

  // stores matrix with entries \int (3 phi^2 - 1) varphi_i varphi_j in M_nonlin_pfield_
  void assemble_M_nonlin_pfield(const size_t cell, const bool restricted)
  {
    const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
    const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
    set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, mu_output_dofs);
    if (restricted)
      set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, phinat_output_dofs);
    XT::Functions::GenericGridFunction<E, 1, 1> M_nonlin_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) { this->bind_phi(cell, element); },
        /*evaluate_func*/
        [cell, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const R phi_n = this->eval_phi(cell, x_local, param);
          return (3. * phi_n * phi_n - 1.);
        });
    M_nonlin_pfield_op_->clear();
    M_nonlin_pfield_op_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(M_nonlin_prefactor)));
    if (!restricted)
      M_nonlin_pfield_op_->assemble(use_tbb_);
    else
      M_nonlin_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
  }

  // stores nonlinear part of block G of the phase field jacobian matrix in G_pfield_nonlinear_part_
  void assemble_G_pfield(const size_t cell, const bool restricted)
  {
    const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
    set_mat_to_zero(G_pfield_, restricted, pfield_submatrix_pattern_, phinat_output_dofs);
    XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) {
          this->bind_phi(cell, element);
          this->bind_mu(cell, element);
          if (this->num_cells_ > 1) {
            for (size_t kk = 0; kk < this->num_cells_; ++kk)
              this->bind_phi(kk, element);
          }
        },
        /*evaluate_func*/
        [cell, In_inv = 1. / In_, eps_inv = 1. / epsilon_, six_inv_Be_eps2 = 6. / (Be_ * std::pow(epsilon_, 2)), this](
            const DomainType& x_local, const XT::Common::Parameter& param) {
          const R phi_n = this->eval_phi(cell, x_local, param);
          const R mu_n = this->eval_mu(cell, x_local, param);
          auto ret = six_inv_Be_eps2 * phi_n * mu_n;
          if (this->num_cells_ > 1) {
            R wsum = 0.;
            R Bsum = 0.;
            for (size_t kk = 0; kk < this->num_cells_; ++kk) {
              if (kk != cell) {
                wsum += this->w_func(kk, x_local, param);
                Bsum += this->B_func(kk, x_local, param);
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
    G_pfield_op_->clear();
    G_pfield_op_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(G_prefactor)));
    if (!restricted)
      G_pfield_op_->assemble(use_tbb_);
    else
      G_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
  }

  // assembles nonlinear part of phasefield residual and adds to residual
  void assemble_nonlinear_part_of_pfield_residual(VectorType& residual, const size_t cell, const bool restricted)
  {
    VectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
    VectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);
    const auto res1 = make_discrete_function(phi_space_, res1_vec);
    const auto res2 = make_discrete_function(phi_space_, res2_vec);
    auto nonlinear_res1_functional = make_vector_functional(phi_space_, res1_vec);
    auto nonlinear_res2_functional = make_vector_functional(phi_space_, res2_vec);
    XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf1(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) {
          this->bind_phi(cell, element);
          this->bind_mu(cell, element);
          if (this->num_cells_ > 1) {
            for (size_t kk = 0; kk < this->num_cells_; ++kk)
              this->bind_phi(kk, element);
          }
        },
        /*evaluate_func*/
        [cell,
         Ca_inv = 1. / Ca_,
         In_inv = 1. / In_,
         eps_inv = 1. / epsilon_,
         num_cells = num_cells_,
         inv_Be_eps2 = 1. / (Be_ * std::pow(epsilon_, 2)),
         this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto phi_n = this->eval_phi(cell, x_local, param);
          const auto mu_n = this->eval_mu(cell, x_local, param);
          auto ret = (Ca_inv + inv_Be_eps2 * (3. * phi_n * phi_n - 1)) * mu_n;
          if (num_cells > 1) {
            R wsum = 0.;
            R Bsum = 0.;
            for (size_t kk = 0; kk < num_cells; ++kk) {
              if (kk != cell) {
                wsum += this->w_func(kk, x_local, param);
                Bsum += this->B_func(kk, x_local, param);
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
        [cell, this](const E& element) { this->bind_phi(cell, element); },
        /*evaluate_func*/
        [cell, inv_eps = 1. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto phi_n = this->eval_phi(cell, x_local, param);
          return inv_eps * (phi_n * phi_n - 1) * phi_n;
        });
    nonlinear_res1_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), nonlinear_res_pf1)));
    nonlinear_res2_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), nonlinear_res_pf2)));
    nonlinear_res1_functional.append(nonlinear_res2_functional);
    if (!restricted)
      nonlinear_res1_functional.assemble(use_tbb_);
    else
      nonlinear_res1_functional.walk_range(pfield_deim_entities_[cell]);
    // high-dimensional operation, TODO: replace
    phi_dirichlet_constraints_.apply(res2_vec);
  }

  //******************************************************************************************************************
  //******************************************* DEIM related methods *************************************************
  //******************************************************************************************************************

  // Dofs needed for evaluation of output_dofs provided in
  std::vector<size_t> pfield_deim_input_dofs(const size_t cell) const
  {
    return pfield_deim_input_dofs_[cell];
  }

  size_t pfield_deim_input_dofs_size(const size_t cell) const
  {
    return pfield_deim_input_dofs_[cell].size();
  }

  // Dofs needed for evaluation of output_dofs provided in
  std::vector<size_t> ofield_deim_input_dofs(const size_t cell) const
  {
    return ofield_deim_input_dofs_[cell];
  }

  // private:
  //******************************************************************************************************************
  //************ The following methods all bind or evaluate the respective temporary discrete function ***************
  //******************************************************************************************************************

  void bind_u(const E& element) const
  {
    auto& u_local = **u_tmp_local_;
    if (!u_local)
      u_local = u_tmp_.local_function();
    u_local->bind(element);
  }

  DomainRetType eval_u(const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    return (**u_tmp_local_)->evaluate(x_local, param);
  }

  JacobianRetType grad_u(const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    return (**u_tmp_local_)->jacobian(x_local, param);
  }

  void bind_P(const size_t cell, const E& element) const
  {
    auto& P_local_cell = (**P_tmp_local_)[cell];
    if (!P_local_cell)
      P_local_cell = P_tmp_[cell].local_function();
    P_local_cell->bind(element);
  }

  DomainRetType eval_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& P_local_cell = (**P_tmp_local_)[cell];
    return P_local_cell->evaluate(x_local, param);
  }

  JacobianRetType grad_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& P_local_cell = (**P_tmp_local_)[cell];
    return P_local_cell->jacobian(x_local, param);
  }

  void bind_Pnat(const size_t cell, const E& element) const
  {
    auto& Pnat_local_cell = (**Pnat_tmp_local_)[cell];
    if (!Pnat_local_cell)
      Pnat_local_cell = Pnat_tmp_[cell].local_function();
    Pnat_local_cell->bind(element);
  }

  DomainRetType eval_Pnat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& Pnat_local_cell = (**Pnat_tmp_local_)[cell];
    return Pnat_local_cell->evaluate(x_local, param);
  }

  void bind_phi(const size_t cell, const E& element) const
  {
    auto& phi_local_cell = (**phi_tmp_local_)[cell];
    if (!phi_local_cell)
      phi_local_cell = phi_tmp_[cell].local_function();
    phi_local_cell->bind(element);
  }

  R eval_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phi_local_cell = (**phi_tmp_local_)[cell];
    return phi_local_cell->evaluate(x_local, param)[0];
  }

  DomainRetType grad_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phi_local_cell = (**phi_tmp_local_)[cell];
    return phi_local_cell->jacobian(x_local, param)[0];
  }

  void bind_phinat(const size_t cell, const E& element) const
  {
    auto& phinat_local_cell = (**phinat_tmp_local_)[cell];
    if (!phinat_local_cell)
      phinat_local_cell = phinat_tmp_[cell].local_function();
    phinat_local_cell->bind(element);
  }

  R eval_phinat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phinat_local_cell = (**phinat_tmp_local_)[cell];
    return phinat_local_cell->evaluate(x_local, param)[0];
  }

  void bind_mu(const size_t cell, const E& element) const
  {
    auto& mu_local_cell = (**mu_tmp_local_)[cell];
    if (!mu_local_cell)
      mu_local_cell = mu_tmp_[cell].local_function();
    mu_local_cell->bind(element);
  }

  R eval_mu(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& mu_local_cell = (**mu_tmp_local_)[cell];
    return mu_local_cell->evaluate(x_local, param)[0];
  }

  //******************************************************************************************************************
  //****************  Linear algebra operations acting only on parts of the given matrices and vectors ***************
  //******************************************************************************************************************

  // Copys given row of source_mat to target_mat.
  // TODO: This functions could be much more efficient for row or col major matrices with the same pattern (which we
  // usually have) by directly copying the values array. This would probably need something like a RowMajorMatrixView.
  // void
  // partially_copy_mat(const std::vector<size_t>& rows, MatrixViewType& target_mat, const MatrixType& source_mat) const
  // {
  //   const auto& pattern = target_mat.get_pattern();
  //   for (const auto& row : rows)
  //     for (const auto& col : pattern.inner(row))
  //       target_mat.set_entry(row, col, source_mat.get_entry(row, col));
  // }

  // Adds given row of rhs to the respective rows of lhs.
  // TODO: This functions could be much more efficient for row or col major matrices with the same pattern by directly
  // modifying the values array. This would probably need something like a RowMajorMatrixView.
  // void partially_add_mats(const std::vector<size_t>& rows, MatrixViewType& lhs, const MatrixType& rhs) const
  // {
  //   const auto& pattern = lhs.get_pattern();
  //   for (const auto& row : rows)
  //     for (const auto& col : pattern.inner(row))
  //       lhs.add_to_entry(row, col, rhs.get_entry(row, col));
  // }

  // Multiplies given rows of mat by alpha.
  // TODO: This functions could be much more efficient for row or col major matrices by directly modifying the values
  // array. This would probably need something like a RowMajorMatrixView.
  // void partially_scal_mat(const std::vector<size_t>& rows, MatrixViewType& mat, const R& alpha) const
  // {
  //   const auto& pattern = mat.get_pattern();
  //   for (const auto& row : rows)
  //     for (const auto& col : pattern.inner(row))
  //       mat.set_entry(row, col, mat.get_entry(row, col) * alpha);
  // }

  // void partially_scal_mat(const std::vector<size_t>& rows,
  //                         MatrixType& mat,
  //                         const XT::LA::SparsityPatternDefault& pattern,
  //                         const R& alpha) const
  // {
  //   for (const auto& row : rows)
  //     for (const auto& col : pattern.inner(row))
  //       mat.set_entry(row, col, mat.get_entry(row, col) * alpha);
  // }

  // Matrix vector multiplication.
  // In the restricted case, only uses the matrix rows provided and does not touch the other entries of range.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(const MatrixType&, const VecType1&, VecType2&, const std::vector<size_t>&)>
  mv_func(const bool restricted) const
  {
    if (!restricted) {
      return [](const MatrixType& mat, const VecType1& source, VecType2& range, const std::vector<size_t>&) {
        mat.mv(source, range);
      };
    } else {
      return [](const MatrixType& mat, const VecType1& source, VecType2& range, const std::vector<size_t>& rows) {
        for (const auto& row : rows) {
          range[row] = 0.;
          for (typename MatrixType::BackendType::InnerIterator it(mat.backend(), row); it; ++it)
            range[row] += it.value() * source[it.col()];
        }
      };
    } // if (restricted)
  }

  // Multiplies given entries of vec by alpha.
  template <class VecType>
  std::function<void(VecType&, const R, const std::vector<size_t>&)> scal_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType& vec, const R alpha, const std::vector<size_t>&) { vec *= alpha; };
    } else {
      return [](VecType& vec, const R alpha, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          vec[dof] *= alpha;
      };
    } // if (restricted)
  }

  // Adds given entries of rhs to respective entries of lhs.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const VecType2&, const std::vector<size_t>&)> add_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>&) { lhs += rhs; };
    } else {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] += rhs[dof];
      };
    } // if (restricted)
  }

  // Subtracts given entries of rhs from respective entries of lhs.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const VecType2&, const std::vector<size_t>&)> sub_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>&) { lhs -= rhs; };
    } else {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] -= rhs[dof];
      };
    } // if (restricted)
  }

  // Computes lhs += alpha * rhs;
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const R, const VecType2&, const std::vector<size_t>&)>
  axpy_func(const bool restricted) const
  {
    if (!restricted) {
      return
          [](VecType1& lhs, const R alpha, const VecType2& rhs, const std::vector<size_t>&) { lhs.axpy(alpha, rhs); };
    } else {
      return [](VecType1& lhs, const R alpha, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] += rhs[dof] * alpha;
      };
    } // if (restricted)
  }

  // Copies low-dimensional vec to given entries of high-dimensional vec.
  void copy_ld_to_hd_vec(const std::vector<size_t> dofs, const VectorType& ld_vec, VectorType& hd_vec)
  {
    assert(ld_vec.size() == dofs.size());
    for (size_t ii = 0; ii < dofs.size(); ++ii)
      hd_vec[dofs[ii]] = ld_vec[ii];
  }

  //******************************************************************************************************************
  //*******************************************  Linear solvers ******************************************************
  //******************************************************************************************************************

#if 0
  VectorType solve_ofield_linear_system(const VectorType& rhs, const size_t cell) const
  {
    //    std::ofstream S_file("S_" + XT::Common::to_string(dt) + ".txt");
    //    S_file << S_ << std::endl;
    //    S_file.close();
    //    DUNE_THROW(NotImplemented, "");
    //    const auto ret = XT::LA::solve(S_, rhs_vector_, XT::LA::SolverOptions<MatrixType>::options("lu.umfpack"));
    //      ofield_update_ = XT::LA::solve(S_, ofield_residual_);
    EigenVectorType update(rhs.size());
    ofield_solver_->compute(S_ofield_.backend());
    const auto rhs_eig = XT::Common::convert_to<EigenVectorType>(rhs);
    update.backend() = ofield_solver_->solveWithGuess(rhs_eig.backend(), ofield_old_result_[cell].backend());
    ofield_old_result_[cell] = update;
    return XT::Common::convert_to<VectorType>(update);
  }
#endif

  VectorType solve_ofield_linear_system(const VectorType& rhs, const size_t cell) const
  {
    VectorType ret(2 * size_u_, 0., 0);
    ConstVectorViewType rhs_P(rhs, 0, size_u_);
    ConstVectorViewType rhs_Pnat(rhs, size_u_, 2 * size_u_);

    // compute P^{n+1} first
    P_tmp_eigen_ = rhs_P;
    P_tmp_eigen_.axpy(-dt_ / kappa_, rhs_Pnat);
    // ofield_schur_solver_->compute(S_schur_ofield_.backend());
    // P_eigen_[cell].backend() = ofield_schur_solver_->solveWithGuess(P_tmp_eigen_.backend(),
    // P_eigen_[cell].backend());
    Dune::InverseOperatorResult res;
    ofield_schur_solver_->apply(P_eigen_[cell], P_tmp_eigen_, res);

    // compute P^{natural,n+1}
    P_tmp_eigen_ = rhs_Pnat;
    C_ofield_linear_part_.mv(P_eigen_[cell], P_tmp_eigen2_);
    P_tmp_eigen_ -= P_tmp_eigen2_;
    C_ofield_nonlinear_part_.mv(P_eigen_[cell], P_tmp_eigen2_);
    P_tmp_eigen_ -= P_tmp_eigen2_;
    P_tmp_eigen_.backend() = ofield_mass_matrix_solver_->solve(P_tmp_eigen_.backend());

    // copy to return vector and return
    VectorViewType ret_P_view_(ret, 0., size_u_);
    VectorViewType ret_Pnat_view_(ret, size_u_, 2 * size_u_);
    ret_P_view_ = P_eigen_[cell];
    ret_Pnat_view_ = P_tmp_eigen_;
    return ret;
  }

  // VectorType solve_pfield_linear_system(const VectorType& rhs, const size_t cell)
  // {
  //   EigenVectorType update(rhs.size());
  //   pfield_solver_->compute(S_pfield_.backend());
  //   const auto rhs_eig = XT::Common::convert_to<EigenVectorType>(rhs);
  //   update.backend() = pfield_solver_->solveWithGuess(rhs_eig.backend(), pfield_old_result_[cell].backend());
  //   pfield_old_result_[cell] = update;
  //   return XT::Common::convert_to<VectorType>(update);
  // }

  VectorType solve_pfield_linear_system(const VectorType& rhs, const size_t cell)
  {
    VectorType ret(3 * size_phi_, 0.);
    VectorViewType ret_phi(ret, 0, size_phi_);
    VectorViewType ret_phinat(ret, size_phi_, 2 * size_phi_);
    VectorViewType ret_mu(ret, 2 * size_phi_, 3 * size_phi_);
    ConstVectorViewType rhs_phi(rhs, 0, size_phi_);
    ConstVectorViewType rhs_phinat(rhs, size_phi_, 2 * size_phi_);
    ConstVectorViewType rhs_mu(rhs, 2 * size_phi_, 3 * size_phi_);
    // apply inverse of first factor to rhs
    // phi part
    phi_tmp_eigen_ = rhs_phinat;
    phi_tmp_eigen_.backend() = pfield_mass_matrix_solver_->solve(phi_tmp_eigen_.backend());
    M_ell_pfield_.mv(phi_tmp_eigen_, phi_tmp_eigen2_);
    phi_tmp_eigen2_ *= -gamma_ * dt_;
    phi_tmp_eigen2_ += rhs_phi;
    pfield_tmp_rhs_phimu_phi_ = phi_tmp_eigen2_;
    // mu part
    pfield_tmp_rhs_phimu_mu_ = rhs_mu;
    Dune::InverseOperatorResult res;
    pfield_phimu_solver_->apply(pfield_tmp_phimu_[cell], pfield_tmp_rhs_phimu_, res);
    ret_phi = pfield_tmp_phimu_phi_[cell];
    ret_mu = pfield_tmp_phimu_mu_[cell];
    // now solve for phinat
    // compute mu part of rhs
    phi_tmp_eigen_ = rhs_phinat;
    G_pfield_.mv(ret_phi, phi_tmp_eigen2_);
    phi_tmp_eigen_ -= phi_tmp_eigen2_;
    M_pfield_.mv(ret_mu, phi_tmp_eigen2_);
    phi_tmp_eigen_.axpy(-1. / Ca_, phi_tmp_eigen2_);
    M_ell_pfield_.mv(ret_mu, phi_tmp_eigen2_);
    phi_tmp_eigen_.axpy(-1. / Be_, phi_tmp_eigen2_);
    M_nonlin_pfield_.mv(ret_mu, phi_tmp_eigen2_);
    phi_tmp_eigen_.axpy(-1. / (Be_ * std::pow(epsilon_, 2)), phi_tmp_eigen2_);
    phi_tmp_eigen_.backend() = pfield_mass_matrix_solver_->solve(phi_tmp_eigen_.backend());
    ret_phinat = phi_tmp_eigen_;
    return ret;
  }


  //***************************************************************************************************************
  //*********************************************  Helper methods  ************************************************
  //***************************************************************************************************************

  // get lower left of computational domain from testcase name
  static XT::Common::FieldVector<R, d> get_lower_left(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return {{0., 0.}};
    else if (testcase == "two_cells")
      return {{0., 0.}};
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return FieldVector<R, d>();
  }

  // get upper right of computational domain from testcase name
  static XT::Common::FieldVector<R, d> get_upper_right(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return {{160., 40.}};
    else if (testcase == "two_cells")
      return {{50., 50.}};
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return FieldVector<R, d>();
  }

  // get directions in which domain is periodic from testcase name
  static std::string get_periodic_directions(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return "01";
    else if (testcase == "two_cells")
      return "00";
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return "";
  }

  // get number of cells from testcase name
  static size_t get_num_cells(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return 1;
    else if (testcase == "two_cells")
      return 2;
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return 0;
  }

  // creates sparsity pattern of orientation field system matrix
  static XT::LA::SparsityPatternDefault create_ofield_pattern(const size_t n,
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

  // creates sparsity pattern of phasefield system matrix
  static XT::LA::SparsityPatternDefault create_pfield_pattern(const size_t n,
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

  // creates sparsity pattern of stokes system matrix
  static XT::LA::SparsityPatternDefault create_stokes_pattern(const SpaceInterface<PGV, d, 1, R>& u_space,
                                                              const SpaceInterface<PGV, 1, 1, R>& p_space)
  {
    const auto pattern_A = make_element_sparsity_pattern(u_space, u_space, u_space.grid_view());
    const auto pattern_B = make_element_sparsity_pattern(u_space, p_space, u_space.grid_view());
    const auto m = u_space.mapper().size();
    const auto n = p_space.mapper().size();
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

  // appends entity to input_entities if one of its global_indices is in output_dofs
  void maybe_add_entity(const E& entity,
                        const DynamicVector<size_t>& global_indices,
                        const std::vector<size_t>& output_dofs,
                        std::vector<E>& input_entities) const
  {
    for (const auto& output_dof : output_dofs) {
      const size_t dof = output_dof % size_phi_;
      for (size_t jj = 0; jj < global_indices.size(); ++jj) {
        if (global_indices[jj] == dof) {
          input_entities.push_back(entity);
          return;
        }
      } // jj
    } // dof
  }

  // sets temporary orientation field discrete functions to source values
  void fill_tmp_ofield(const size_t cell, const VectorType& source, const bool restricted = false) const
  {
    if (!restricted) {
      ConstVectorViewType P_vec(source, 0, size_u_);
      P_tmp_[cell].dofs().vector() = P_vec;
    } else {
      const auto& input_dofs = ofield_deim_input_dofs_[cell];
      for (size_t ii = 0; ii < Pnat_deim_input_dofs_begin_[cell]; ++ii)
        P_tmp_[cell].dofs().vector().set_entry(input_dofs[ii], source[input_dofs[ii]]);
    }
  }

  // sets temporary phase field discrete functions to source values
  void fill_tmp_pfield(const size_t cell, const VectorType& source, const bool restricted) const
  {
    if (!restricted) {
      ConstVectorViewType phi_vec(source, 0, size_phi_);
      ConstVectorViewType mu_vec(source, 2 * size_phi_, 3 * size_phi_);
      phi_tmp_[cell].dofs().vector() = phi_vec;
      mu_tmp_[cell].dofs().vector() = mu_vec;
    } else {
      const auto& input_dofs = pfield_deim_input_dofs_[cell];
      for (size_t ii = 0; ii < phinat_deim_input_dofs_begin_[cell]; ++ii)
        phi_tmp_[cell].dofs().vector().set_entry(input_dofs[ii], source[input_dofs[ii]]);
      for (size_t ii = mu_deim_input_dofs_begin_[cell]; ii < input_dofs.size(); ++ii)
        mu_tmp_[cell].dofs().vector().set_entry(input_dofs[ii] - 2 * size_phi_, source[input_dofs[ii]]);
    }
  }

  // error norm used in orientation field Newton iteration
  // TODO: use appropriate norm
  double ofield_residual_norm(const VectorType& residual, double l2_ref_P, double l2_ref_Pnat) const
  {
    l2_ref_P = l2_ref_P < 1. ? 1. : l2_ref_P;
    l2_ref_Pnat = l2_ref_Pnat < 1. ? 1. : l2_ref_Pnat;
    ConstVectorViewType res0_vec(residual, 0, size_u_);
    ConstVectorViewType res1_vec(residual, size_u_, 2 * size_u_);
    const auto res0 = make_discrete_function(u_space_, res0_vec);
    const auto res1 = make_discrete_function(u_space_, res1_vec);
    return l2_norm(grid_view_, res0) / l2_ref_P + l2_norm(grid_view_, res1) / l2_ref_Pnat;
  }

  // error norm used in phase field Newton iteration
  // TODO: use appropriate norm
  double
  pfield_residual_norm(const VectorType& residual, double l2_ref_phi, double l2_ref_phinat, double l2_ref_mu) const
  {
    l2_ref_phi = l2_ref_phi < 1. ? 1. : l2_ref_phi;
    l2_ref_phinat = l2_ref_phinat < 1. ? 1. : l2_ref_phinat;
    l2_ref_mu = l2_ref_mu < 1. ? 1. : l2_ref_mu;
    ConstVectorViewType res0_vec(residual, 0, size_phi_);
    ConstVectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
    ConstVectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);
    const auto res0 = make_discrete_function(phi_space_, res0_vec);
    const auto res1 = make_discrete_function(phi_space_, res1_vec);
    const auto res2 = make_discrete_function(phi_space_, res2_vec);
    return l2_norm(grid_view_, res0) / l2_ref_phi + l2_norm(grid_view_, res1) / l2_ref_phinat
           + l2_norm(grid_view_, res2) / l2_ref_mu;
  }

  R B_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    const R phi_n = eval_phi(kk, x_local, param);
    return 1 / epsilon_ * std::pow(std::pow(phi_n, 2) - 1, 2);
  }

  R w_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    const R phi_n = eval_phi(kk, x_local, param);
    if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
      return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
    else
      return 0.;
  }

  //******************************************************************************************************************
  //*******************************************  Member variables ****************************************************
  //******************************************************************************************************************

  // Model parameters
  XT::Common::FieldVector<R, d> lower_left_;
  XT::Common::FieldVector<R, d> upper_right_;
  const double t_end_;
  double t_;
  const bool use_tbb_;
  const double Re_;
  double Fa_inv_;
  double xi_;
  double kappa_;
  double c_1_;
  double Pa_;
  double last_pfield_Pa_;
  double last_ofield_Pa_;
  double beta_;
  double gamma_;
  double Be_;
  double Ca_;
  double epsilon_;
  double In_;
  const double vol_domain_;
  const size_t num_cells_;
  const bool linearize_;
  // Grid and grid views
  XT::Grid::GridProvider<G> grid_;
  const GV nonperiodic_grid_view_;
  const PGV grid_view_;
  // Finite element function spaces
  const ContinuousLagrangeSpace<PGV, d, R> u_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> p_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> phi_space_;
  // Size of finite element vectors
  const size_t size_u_;
  const size_t size_p_;
  const size_t size_phi_;
  const size_t num_mutexes_u_;
  const size_t num_mutexes_ofield_;
  const size_t num_mutexes_pfield_;
  // Finite element vectors for phase field, orientation field and stokes system
  // There is one phase field and orientation field per cell
  VectorType stokes_vector_;
  std::vector<VectorType> ofield_vectors_;
  std::vector<VectorType> pfield_vectors_;
  // Views on parts of the system vectors corresponding to the respective variables
  VectorViewType u_view_;
  VectorViewType p_view_;
  std::vector<VectorViewType> P_view_;
  std::vector<VectorViewType> Pnat_view_;
  std::vector<VectorViewType> phi_view_;
  std::vector<VectorViewType> phinat_view_;
  std::vector<VectorViewType> mu_view_;
  // DiscreteFunctions corresponding to view vectors
  ViewVectorDiscreteFunctionType u_;
  ViewDiscreteFunctionType p_;
  std::vector<ViewVectorDiscreteFunctionType> P_;
  std::vector<ViewVectorDiscreteFunctionType> Pnat_;
  std::vector<ViewDiscreteFunctionType> phi_;
  std::vector<ViewDiscreteFunctionType> phinat_;
  std::vector<ViewDiscreteFunctionType> mu_;
  // Stokes system matrix S = (A B; B^T 0) and views on matrix blocks
  MatrixType S_stokes_;
  MatrixViewType A_stokes_;
  MatrixViewType B_stokes_;
  MatrixViewType BT_stokes_;
  // pressure mass matrix
  MatrixType M_p_stokes_;
  // Matrix operator for A_stokes_
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, d>> A_stokes_op_;
  // finite element vector rhs = (f; g) for stokes system and views on velocity and pressure parts f and g
  EigenVectorType stokes_rhs_vector_;
  EigenVectorViewType stokes_f_vector_;
  EigenVectorViewType stokes_g_vector_;
  // vector containing integrals of pressure basis functions (for normalizing such that \int p = 0)
  VectorType p_basis_integrated_vector_;
  // Dirichlet constraints
  const XT::Grid::AllDirichletBoundaryInfo<PI> boundary_info_;
  DirichletConstraints<PI, SpaceInterface<PGV, d, 1, R>> u_dirichlet_constraints_;
  PhiDirichletConstraintsType phi_dirichlet_constraints_;
  // Sparsity pattern of one block of orientation field system matrix
  XT::LA::SparsityPatternDefault ofield_submatrix_pattern_;
  // Orientation field system matrix S = (M/dt+A B; C D)
  // MatrixType S_ofield_;
  // mutable MatrixViewType S_ofield_00_;
  // mutable MatrixViewType S_ofield_01_;
  // mutable MatrixViewType S_ofield_10_;
  // mutable MatrixViewType S_ofield_11_;
  // Orientation field mass matrix
  MatrixType M_ofield_;
  MatrixType A_ofield_;
  // Part of C that is independent of phi and dt
  MatrixType C_ofield_elliptic_part_;
  // Whole linear part of C (elliptic part + phi-dependent part)
  MatrixType C_ofield_linear_part_;
  // Whole linear part of C (elliptic part + phi-dependent part)
  mutable MatrixType C_ofield_nonlinear_part_;
  // Linear part of ofield schur matrix M/dt + A - 1/kappa C
  MatrixType S_schur_ofield_linear_part_;
  mutable MatrixType S_schur_ofield_;
  // Matrix operators for orientation field matrices
  std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> M_ofield_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> A_ofield_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> C_ofield_linear_part_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> C_ofield_nonlinear_part_op_;
  OfieldMatrixLinearPartOperator<VectorType, MatrixType, CellModelSolver> ofield_jac_linear_op_;
  OfieldSchurComplementOperator<EigenVectorType, MatrixType> ofield_schur_op_;
  // finite element vector rhs = (f; g) for ofield system and views on P and Pnat parts f and g
  VectorType ofield_rhs_vector_;
  VectorViewType ofield_f_vector_;
  VectorViewType ofield_g_vector_;
  mutable std::vector<EigenVectorType> P_eigen_;
  mutable EigenVectorType P_tmp_eigen_;
  mutable EigenVectorType P_tmp_eigen2_;
  // Vectors for orientation field Newton scheme
  // TODO: do we need all these vectors? Do we need the same vectors for pfield or can we use the same ones for both?
  // mutable std::vector<EigenVectorType> ofield_old_result_;
  // VectorType ofield_residual_;
  // Eigen column major matrix backend, needed for direct solvers
  mutable ColMajorBackendType S_colmajor_;
  mutable ColMajorBackendType M_ofield_colmajor_;
  mutable ColMajorBackendType S_schur_ofield_colmajor_;
  // Linear solvers and linear operators needed for solvers
  std::shared_ptr<DirectSolverType> stokes_solver_;
  mutable std::shared_ptr<SolverType> ofield_solver_;
  mutable std::shared_ptr<DirectSolverType> ofield_mass_matrix_solver_;
  XT::LA::IdentityPreconditioner<PfieldPhiMuMatrixOperatorType> identity_prec_;
  mutable std::shared_ptr<OfieldSchurSolverType> ofield_schur_solver_;
  VectorType ofield_tmp_vec_;
  // Indices for restricted operator in DEIM context
  std::vector<std::vector<size_t>> ofield_deim_input_dofs_;
  std::vector<size_t> Pnat_deim_input_dofs_begin_;
  // output dofs that were computed by the DEIM algorithm
  std::vector<std::shared_ptr<std::vector<size_t>>> ofield_deim_output_dofs_;
  std::vector<std::vector<size_t>> ofield_deim_unique_output_dofs_;
  // output dofs for the respective variables, shifted to range [0, size_phi)
  std::vector<std::vector<size_t>> P_deim_output_dofs_;
  std::vector<std::vector<size_t>> Pnat_deim_output_dofs_;
  // Entities that we have to walk over to calculate values at output dofs
  std::vector<std::vector<E>> ofield_deim_entities_;
  // DiscreteFunctions and vectors to be used in (nonlinear) calculations where a source vector is provided

  // Sparsity pattern of one block of phase field system matrix
  const XT::LA::SparsityPatternDefault pfield_submatrix_pattern_;
  // Phase field system matrix S = (M/dt+D E 0; G H J; A 0 C)
  // MatrixType S_pfield_;
  // MatrixViewType S_pfield_00_;
  // MatrixViewType S_pfield_01_;
  // MatrixViewType S_pfield_10_;
  // MatrixViewType S_pfield_11_;
  // MatrixViewType S_pfield_12_;
  // MatrixViewType S_pfield_20_;
  // MatrixViewType S_pfield_22_;
  MatrixType M_pfield_;
  MatrixType D_pfield_;
  MatrixType M_ell_pfield_;
  MatrixType M_nonlin_pfield_;
  MatrixType G_pfield_;
  mutable ColMajorBackendType M_pfield_colmajor_;
  // Pfield solvers and linear operators
  mutable std::shared_ptr<SolverType> pfield_solver_;
  mutable std::shared_ptr<DirectSolverType> pfield_mass_matrix_solver_;
  PfieldPhiMuMatrixOperatorType pfield_phimu_matrixop_;
  PfieldMatrixLinearPartOperator<VectorType, MatrixType, PhiDirichletConstraintsType, CellModelSolver>
      pfield_jac_linear_op_;
  mutable std::shared_ptr<PfieldPhiMuSolverType> pfield_phimu_solver_;
  // Matrix operators for phasefield matrices
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, 1>> S_pfield_00_op_;
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, 1>> S_pfield_10_op_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> D_pfield_op_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> G_pfield_op_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> M_nonlin_pfield_op_;
  // Phase field rhs vector (g h f)
  VectorType pfield_rhs_vector_;
  VectorViewType pfield_g_vector_;
  VectorViewType pfield_h_vector_;
  VectorViewType pfield_f_vector_;
  EigenVectorType phi_tmp_eigen_;
  EigenVectorType phi_tmp_eigen2_;
  std::vector<EigenVectorType> pfield_tmp_phimu_;
  std::vector<EigenVectorViewType> pfield_tmp_phimu_phi_;
  std::vector<EigenVectorViewType> pfield_tmp_phimu_mu_;
  EigenVectorType pfield_tmp_rhs_phimu_;
  EigenVectorViewType pfield_tmp_rhs_phimu_phi_;
  EigenVectorViewType pfield_tmp_rhs_phimu_mu_;
  // Indices for restricted operator in DEIM context
  std::vector<std::vector<size_t>> pfield_deim_input_dofs_;
  // phinat_deim_input_dofs_begin_[cell] contains index of first phinat input dof in pfield_deim_input_dofs_[cell]
  // vector
  std::vector<size_t> phinat_deim_input_dofs_begin_;
  std::vector<size_t> mu_deim_input_dofs_begin_;
  // output dofs that were computed by the DEIM algorithm
  std::vector<std::shared_ptr<std::vector<size_t>>> pfield_deim_output_dofs_;
  std::vector<std::vector<size_t>> pfield_deim_unique_output_dofs_;
  // output dofs for the respective variables, shifted to range [0, size_phi)
  std::vector<std::vector<size_t>> phi_deim_output_dofs_;
  std::vector<std::vector<size_t>> phinat_deim_output_dofs_;
  std::vector<std::vector<size_t>> mu_deim_output_dofs_;
  // Indices that are both in phinat and mu output dofs
  std::vector<std::vector<size_t>> both_mu_and_phi_deim_output_dofs_;
  // Entities that we have to walk over to calculate values at output dofs
  std::vector<std::vector<E>> pfield_deim_entities_;
  // DiscreteFunctions and vectors to be used in (nonlinear) calculations where a source vector is provided
  VectorType pfield_tmp_vec_;
  VectorType pfield_tmp_vec2_;
  VectorType phi_tmp_vec_;
  VectorType u_tmp_vec_;
  mutable VectorDiscreteFunctionType u_tmp_;
  mutable std::vector<VectorDiscreteFunctionType> P_tmp_;
  mutable std::vector<VectorDiscreteFunctionType> Pnat_tmp_;
  mutable std::vector<DiscreteFunctionType> phi_tmp_;
  mutable std::vector<DiscreteFunctionType> phinat_tmp_;
  mutable std::vector<DiscreteFunctionType> mu_tmp_;
  mutable std::shared_ptr<PerThreadVectorLocalFunc> u_tmp_local_;
  mutable std::shared_ptr<PerThreadVectorLocalFuncs> P_tmp_local_;
  mutable std::shared_ptr<PerThreadVectorLocalFuncs> Pnat_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> phi_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> phinat_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> mu_tmp_local_;
  double dt_;
};

#endif // DUNE_GDT_EXAMPLES_CELLMODEL_HH
