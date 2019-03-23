// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_INTERNAL_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_INTERNAL_HH

#include <boost/multi_array.hpp>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "slopes.hh"

namespace Dune {
namespace GDT {
namespace internal {


template <class MatImp>
XT::Common::Configuration hyperbolic_default_eigensolver_options()
{
  using EigenSolverOptionsType = typename XT::LA::EigenSolverOptions<MatImp>;
  using MatrixInverterOptionsType = typename XT::LA::MatrixInverterOptions<MatImp>;
  XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
  //  XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
  eigensolver_options["assert_eigendecomposition"] = "1e-6";
  eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
  eigensolver_options["disable_checks"] =
#ifdef NDEBUG
      "true";
#else
      "false";
#endif
  XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
  matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
  matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
  eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
  return eigensolver_options;
} // ... hyperbolic_default_eigensolver_options()


// Wrapper for thread-safe and consistent handling of the different jacobians (Usual matrices vs. block matrices in the
// partial moment case)
template <class AnalyticalFluxType, class MatrixImp, class VectorImp>
class EigenvectorWrapperBase
{
public:
  using MatrixType = MatrixImp;
  using VectorType = VectorImp;
  static constexpr size_t dimDomain = AnalyticalFluxType::r;
  static constexpr size_t dimRange = AnalyticalFluxType::rC;
  using RangeFieldType = typename AnalyticalFluxType::R;
  using DomainType = FieldVector<RangeFieldType, dimDomain>;
  using E = typename AnalyticalFluxType::E;
  using LocalFluxType = typename AnalyticalFluxType::LocalFunctionType;
  using FluxDomainType = XT::Common::FieldVector<typename AnalyticalFluxType::D, dimRange>;

  EigenvectorWrapperBase(const AnalyticalFluxType& analytical_flux, const bool flux_is_affine)
    : analytical_flux_(analytical_flux)
    , local_flux_(analytical_flux_.local_function())
    , flux_is_affine_(flux_is_affine)
    , computed_(false)
  {}

  virtual ~EigenvectorWrapperBase() {}

  virtual void compute_eigenvectors(const E& entity,
                                    const DomainType& x_local,
                                    const VectorType& u,
                                    const XT::Common::Parameter& param)
  {
    if (!computed_ || !flux_is_affine_) {
      compute_eigenvectors_impl(entity, x_local, u, param);
      computed_ = true;
    }
  }

  virtual void compute_eigenvectors_impl(const E& entity,
                                         const DomainType& x_local,
                                         const VectorType& u,
                                         const XT::Common::Parameter& param) = 0;

  virtual void apply_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const = 0;

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const = 0;

  virtual const MatrixType& eigenvectors(const size_t dd) const = 0;

protected:
  const AnalyticalFluxType& analytical_flux_;
  const std::unique_ptr<LocalFluxType> local_flux_;
  const bool flux_is_affine_;
  bool computed_;
}; // class EigenvectorWrapperBase<...>


template <class AnalyticalFluxType,
          class MatrixType =
              FieldMatrix<typename AnalyticalFluxType::R, AnalyticalFluxType::rC, AnalyticalFluxType::rC>,
          class VectorType = FieldVector<typename AnalyticalFluxType::RangeFieldType, AnalyticalFluxType::rC>>
class EigenvectorWrapper : public EigenvectorWrapperBase<AnalyticalFluxType, MatrixType, VectorType>
{
  using BaseType = EigenvectorWrapperBase<AnalyticalFluxType, MatrixType, VectorType>;

protected:
  using V = XT::Common::VectorAbstraction<VectorType>;
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using EigenSolverType = typename XT::LA::EigenSolver<MatrixType>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::RangeFieldType;
  using JacobianType = XT::Common::FieldVector<MatrixType, dimDomain>;

  EigenvectorWrapper(const AnalyticalFluxType& analytical_flux, const bool flux_is_affine)
    : BaseType(analytical_flux, flux_is_affine)
    , work_(1)
    , scale_(dimRange)
    , rconde_(dimRange)
    , rcondv_(dimRange)
    , iwork_(2 * dimRange - 2)
    , jacobian_(std::make_unique<JacobianType>())
    , eigenvectors_(std::make_unique<JacobianType>())
    , eigenvectors_rcond_(1.)
    , eigenvalues_(std::vector<RangeFieldType>(dimRange))
    , QR_(std::make_unique<JacobianType>())
    , tau_(V::create(dimRange))
  {
#if HAVE_MKL || HAVE_LAPACKE
    int ilo, ihi;
    double norm;
    if (M::storage_layout == XT::Common::StorageLayout::dense_row_major) {
      // get optimal working size in work[0] (requested by lwork = -1)
      int info = XT::Common::Lapacke::dgeevx_work(XT::Common::Lapacke::row_major(),
                                                  /*both diagonally scale and permute*/ 'B',
                                                  /*do_not_compute_left_eigenvectors:*/ 'N',
                                                  /*compute_right_eigenvectors:*/ 'V',
                                                  /*do not compute condition numbers*/ 'N',
                                                  static_cast<int>(dimRange),
                                                  M::data((*jacobian_)[0]),
                                                  static_cast<int>(dimRange),
                                                  &(eigenvalues_[0][0]),
                                                  &(dummy_complex_eigenvalues_[0]),
                                                  nullptr,
                                                  static_cast<int>(dimRange),
                                                  M::data((*eigenvectors_)[0]),
                                                  static_cast<int>(dimRange),
                                                  &ilo,
                                                  &ihi,
                                                  scale_.data(),
                                                  &norm,
                                                  rconde_.data(),
                                                  rcondv_.data(),
                                                  work_.data(),
                                                  -1,
                                                  iwork_.data());
      if (info != 0)
        DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
      work_.resize(static_cast<size_t>(work_[0] + 0.5));
    }
#endif
  }

  virtual void compute_eigenvectors_impl(const E& entity,
                                         const DomainType& x_local,
                                         const VectorType& u,
                                         const XT::Common::Parameter& param) override final
  {
    local_flux_->bind(entity);
    *jacobian_ = local_flux_->jacobian(x_local, u, param);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      try {
        if (false) {
          ;
#if HAVE_MKL || HAVE_LAPACKE
        } else if (M::storage_layout == XT::Common::StorageLayout::dense_row_major) {
          int ilo, ihi;
          double norm;
          int info = XT::Common::Lapacke::dgeevx_work(XT::Common::Lapacke::row_major(),
                                                      /*both diagonally scale and permute*/ 'B',
                                                      /*do not compute left eigenvectors*/ 'N',
                                                      /*compute right eigenvectors*/ 'V',
                                                      /*do not compute condition numbers*/ 'N',
                                                      static_cast<int>(dimRange),
                                                      M::data((*jacobian_)[dd]),
                                                      static_cast<int>(dimRange),
                                                      eigenvalues_[dd].data(),
                                                      &(dummy_complex_eigenvalues_[0]),
                                                      nullptr,
                                                      static_cast<int>(dimRange),
                                                      M::data((*eigenvectors_)[dd]),
                                                      static_cast<int>(dimRange),
                                                      &ilo,
                                                      &ihi,
                                                      scale_.data(),
                                                      &norm,
                                                      rconde_.data(),
                                                      rcondv_.data(),
                                                      work_.data(),
                                                      static_cast<int>(work_.size()),
                                                      iwork_.data());
          if (info != 0)
            DUNE_THROW(Dune::MathError, "The lapack backend reported '" << info << "'!");
#endif // HAVE_MKL || HAVE_LAPACKE
        } else {
          static auto eigensolver_options = hyperbolic_default_eigensolver_options<MatrixType>();
          const auto eigensolver = EigenSolverType((*jacobian_)[dd], &eigensolver_options);
          (*eigenvectors_)[dd] = eigensolver.real_eigenvectors();
          eigenvalues_[dd] = eigensolver.real_eigenvalues();
        }
        (*QR_)[dd] = (*eigenvectors_)[dd];
        XT::LA::qr((*QR_)[dd], tau_[dd], permutations_[dd]);
#if HAVE_MKL || HAVE_LAPACKE
        int info = XT::Common::Lapacke::dtrcon(XT::Common::Lapacke::row_major(),
                                               '1',
                                               'U',
                                               'N',
                                               static_cast<int>(dimRange),
                                               M::data((*QR_)[dd]),
                                               static_cast<int>(dimRange),
                                               &eigenvectors_rcond_);
        if (info || eigenvectors_rcond_ < 1e-5)
          DUNE_THROW(Dune::MathError, "Eigenvector condition too high!");
#endif
      } catch (const Dune::MathError&) {
        // use scalar limiters, i.e. eigenvectors matrix is eye-matrix.
        XT::LA::eye_matrix((*eigenvectors_)[dd]);
        std::fill(eigenvalues_[dd].begin(), eigenvalues_[dd].end(), 1.);
        (*QR_)[dd] = (*eigenvectors_)[dd];
        XT::LA::qr((*QR_)[dd], tau_[dd], permutations_[dd]);
      }
    } // dd
    // we do not need jacobian_ anymore if the flux is affine, so save memory
    if (flux_is_affine_)
      jacobian_ = nullptr;
  }

  virtual void apply_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const override final
  {
    (*eigenvectors_)[dd].mv(u, ret);
  }

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const override final
  {
    thread_local VectorType work = V::create(dimRange);
    XT::LA::solve_qr_factorized((*QR_)[dd], tau_[dd], permutations_[dd], ret, u, &work);
  }

  virtual const MatrixType& eigenvectors(const size_t dd) const override final
  {
    return (*eigenvectors_)[dd];
  }

protected:
  using BaseType::computed_;
  using BaseType::flux_is_affine_;
  using BaseType::local_flux_;
  std::vector<RangeFieldType> work_, scale_, rconde_, rcondv_;
  std::vector<int> iwork_;
  std::unique_ptr<JacobianType> jacobian_;
  std::unique_ptr<JacobianType> eigenvectors_;
  RangeFieldType eigenvectors_rcond_;
  FieldVector<std::vector<RangeFieldType>, dimDomain> eigenvalues_;
  FieldVector<RangeFieldType, dimRange> dummy_complex_eigenvalues_;
  std::unique_ptr<JacobianType> QR_;
  FieldVector<VectorType, dimDomain> tau_;
  FieldVector<FieldVector<int, dimRange>, dimDomain> permutations_;
}; // class EigenvectorWrapper<...>


template <class AnalyticalFluxType, size_t block_size = (AnalyticalFluxType::r == 1) ? 2 : 4>
class BlockedEigenvectorWrapper
  : public EigenvectorWrapperBase<
        AnalyticalFluxType,
        XT::Common::BlockedFieldMatrix<typename AnalyticalFluxType::R, AnalyticalFluxType::rC / block_size, block_size>,
        XT::Common::BlockedFieldVector<typename AnalyticalFluxType::R, AnalyticalFluxType::rC / block_size, block_size>>
{
  using BaseType = EigenvectorWrapperBase<
      AnalyticalFluxType,
      XT::Common::BlockedFieldMatrix<typename AnalyticalFluxType::R, AnalyticalFluxType::rC / block_size, block_size>,
      XT::Common::BlockedFieldVector<typename AnalyticalFluxType::R, AnalyticalFluxType::rC / block_size, block_size>>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::FluxDomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::VectorType;
  using JacobianType = XT::Common::FieldVector<MatrixType, dimDomain>;
  static constexpr size_t num_blocks = VectorType::num_blocks;
  static_assert(dimRange % block_size == 0, "dimRange has to be a multiple of block_size");
  using BlockedIntVectorType = XT::Common::BlockedFieldVector<int, num_blocks, block_size>;
  using LocalVectorType = typename VectorType::BlockType;
  using LocalMatrixType = typename MatrixType::BlockType;
  using EigenSolverType = typename XT::LA::EigenSolver<LocalMatrixType>;
  using LocalM = typename XT::Common::MatrixAbstraction<LocalMatrixType>;
  using LocalV = typename XT::Common::VectorAbstraction<LocalVectorType>;
  using NonblockedJacobianType =
      XT::Common::FieldVector<XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>;

  BlockedEigenvectorWrapper(const AnalyticalFluxType& analytical_flux, const bool flux_is_affine)
    : BaseType(analytical_flux, flux_is_affine)
    , jacobian_(std::make_unique<JacobianType>())
    , nonblocked_jacobian_(std::make_unique<NonblockedJacobianType>())
  {
    std::fill_n(&(eigenvalues_[0][0]), dimDomain * num_blocks, std::vector<double>(block_size, 0.));
  }

  virtual void compute_eigenvectors_impl(const E& entity,
                                         const DomainType& x_local,
                                         const VectorType& u,
                                         const XT::Common::Parameter& param) override final
  {
    local_flux_->bind(entity);
    const FluxDomainType nonblocked_u = u.operator FluxDomainType();
    *jacobian_ = local_flux_->jacobian(x_local, nonblocked_u, param);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        try {
          if (block_size == 2) {
            const auto& jac = (*jacobian_)[dd].block(jj);
            const auto trace = jac[0][0] + jac[1][1];
            const auto det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
            const auto sqrt_val = std::sqrt(0.25 * trace * trace - det);
            auto& eigvals = eigenvalues_[dd][jj];
            eigvals[0] = 0.5 * trace + sqrt_val;
            eigvals[1] = 0.5 * trace - sqrt_val;
            auto& eigvecs = eigenvectors_[dd].block(jj);
            if (std::abs(jac[1][0]) > std::abs(jac[0][1])) {
              if (XT::Common::FloatCmp::ne(jac[1][0], 0.)) {
                eigvecs[0][0] = eigvals[0] - jac[1][1];
                eigvecs[0][1] = eigvals[1] - jac[1][1];
                eigvecs[1][0] = eigvecs[1][1] = jac[1][0];
              } else {
                eigvecs[0][0] = eigvecs[1][1] = 1.;
                eigvecs[0][1] = eigvecs[1][0] = 0.;
              }
            } else {
              if (XT::Common::FloatCmp::ne(jac[0][1], 0.)) {
                eigvecs[1][0] = eigvals[0] - jac[0][0];
                eigvecs[1][1] = eigvals[1] - jac[0][0];
                eigvecs[0][0] = eigvecs[0][1] = jac[0][1];
              } else {
                eigvecs[0][0] = eigvecs[1][1] = 1.;
                eigvecs[0][1] = eigvecs[1][0] = 0.;
              }
            }
            // normalize such that the eigenvectors have norm 1
            for (size_t col = 0; col < 2; ++col) {
              RangeFieldType two_norm = 0;
              two_norm = std::sqrt(std::pow(eigvecs[0][col], 2) + std::pow(eigvecs[1][col], 2));
              for (size_t row = 0; row < 2; ++row)
                eigvecs[row][col] /= two_norm;
            }
          } else {
            // For the small matrices (usually 4x4) used here it causes a lot of overhead to call into LAPACK every
            // time, so we just use our own eigensolver most of the time. Occasionally, however, our eigensolver fails
            // where the LAPACK eigensolver succeeds (due to a superior shifting strategy), so in these cases we call
            // LAPACK.
            try {
              XT::LA::internal::fmatrix_compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(
                  (*jacobian_)[dd].block(jj), eigenvalues_[dd][jj], eigenvectors_[dd].block(jj));
            } catch (const Dune::MathError&) {
              // Our own eigensolver failed, try the default one instead (Lapacke, Eigen or Numpy, if none of these is
              // available, we solve again using our own eigensolver, which will throw the error again.
              static auto eigensolver_options = hyperbolic_default_eigensolver_options<LocalMatrixType>();
              const auto eigensolver = EigenSolverType((*jacobian_)[dd].block(jj), &eigensolver_options);
              eigenvectors_[dd].block(jj) = eigensolver.real_eigenvectors();
              eigenvalues_[dd][jj] = eigensolver.real_eigenvalues();
            }
          } // else (block_size == 2)
          QR_[dd].block(jj) = eigenvectors_[dd].block(jj);
          XT::LA::qr(QR_[dd].block(jj), tau_[dd].block(jj), permutations_[dd].block(jj));
#if HAVE_MKL || HAVE_LAPACKE
          int info = XT::Common::Lapacke::dtrcon(XT::Common::Lapacke::row_major(),
                                                 '1',
                                                 'U',
                                                 'N',
                                                 static_cast<int>(block_size),
                                                 LocalM::data(QR_[dd].block(jj)),
                                                 static_cast<int>(block_size),
                                                 &eigenvectors_rcond_);
          if (info || eigenvectors_rcond_ < 1e-5)
            DUNE_THROW(Dune::MathError, "Eigenvector condition too high!");
#endif
        } catch (const Dune::MathError&) {
          // use scalar limiters, i.e. eigenvectors matrix is eye-matrix.
          XT::LA::eye_matrix(eigenvectors_[dd].block(jj));
          std::fill(eigenvalues_[dd][jj].begin(), eigenvalues_[dd][jj].end(), 1.);
          QR_[dd].block(jj) = eigenvectors_[dd].block(jj);
          XT::LA::qr(QR_[dd].block(jj), tau_[dd].block(jj), permutations_[dd].block(jj));
        }
      } // jj
    } // dd
    if (flux_is_affine_) {
      jacobian_ = nullptr;
      nonblocked_jacobian_ = nullptr;
    }
  }

  virtual void apply_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const override final
  {
    eigenvectors_[dd].mv(u, ret);
  }

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& u, VectorType& ret) const override final
  {
    LocalVectorType work;
    for (size_t jj = 0; jj < num_blocks; ++jj)
      XT::LA::solve_qr_factorized(
          QR_[dd].block(jj), tau_[dd].block(jj), permutations_[dd].block(jj), ret.block(jj), u.block(jj), &work);
  }

  virtual const MatrixType& eigenvectors(const size_t dd) const override final
  {
    return eigenvectors_[dd];
  }

protected:
  using BaseType::flux_is_affine_;
  using BaseType::local_flux_;
  std::unique_ptr<JacobianType> jacobian_;
  std::unique_ptr<NonblockedJacobianType> nonblocked_jacobian_;
  FieldVector<FieldVector<std::vector<double>, num_blocks>, dimDomain> eigenvalues_;
  FieldVector<MatrixType, dimDomain> eigenvectors_;
  FieldVector<MatrixType, dimDomain> QR_;
  FieldVector<VectorType, dimDomain> tau_;
  RangeFieldType eigenvectors_rcond_;
  FieldVector<BlockedIntVectorType, dimDomain> permutations_;
}; // BlockedEigenvectorWrapper<...>


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_INTERNAL_HH
