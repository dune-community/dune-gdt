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
} // ... create_eigensolver_opts()


// Wrapper for thread-safe and consistent handling of the different jacobians (Usual matrices vs. block matrices in the
// partial moment case)
template <class AnalyticalFluxType, class MatrixImp, class VectorImp>
class JacobianWrapperBase
{
public:
  using MatrixType = MatrixImp;
  using VectorType = VectorImp;
  static constexpr size_t dimDomain = AnalyticalFluxType::dimDomain;
  static constexpr size_t dimRange = AnalyticalFluxType::dimRange;
  using DomainType = typename AnalyticalFluxType::DomainType;
  using RangeFieldType = typename AnalyticalFluxType::RangeFieldType;
  using EntityType = typename AnalyticalFluxType::EntityType;
  using JacobianType = FieldVector<MatrixType, dimDomain>;
  using StateRangeType = typename AnalyticalFluxType::StateRangeType;

  JacobianWrapperBase()
    : computed_(false)
  {
  }

  virtual ~JacobianWrapperBase()
  {
  }

  virtual void get_jacobian(const size_t dd,
                            const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const VectorType& u,
                            const XT::Common::Parameter& param) = 0;

  virtual void get_jacobian(const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const VectorType& u,
                            const XT::Common::Parameter& param) = 0;

  virtual void compute(const size_t dd) = 0;

  virtual void compute()
  {
    for (size_t dd = 0; dd < dimDomain; ++dd)
      compute(dd);
  }

  virtual bool computed(const size_t dd) const
  {
    return computed_[dd];
  }

  virtual bool computed() const
  {
    for (size_t dd = 0; dd < dimDomain; ++dd)
      if (!computed(dd))
        return false;
    return true;
  }

  virtual std::unique_ptr<JacobianType>& jacobian()
  {
    return jacobian_;
  }

  virtual const std::unique_ptr<JacobianType>& jacobian() const
  {
    return jacobian_;
  }

  virtual MatrixType& jacobian(const size_t dd)
  {
    return (*jacobian_)[dd];
  }

  virtual const MatrixType& jacobian(const size_t dd) const
  {
    return (*jacobian_)[dd];
  }

  virtual void apply_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const = 0;

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const = 0;

  virtual const MatrixType& eigenvectors(const size_t dd) const = 0;

protected:
  std::unique_ptr<JacobianType> jacobian_;
  FieldVector<bool, dimDomain> computed_;
}; // class JacobianWrapperBase<...>


template <class AnalyticalFluxType,
          class MatrixType = FieldMatrix<typename AnalyticalFluxType::RangeFieldType,
                                         AnalyticalFluxType::dimRange,
                                         AnalyticalFluxType::dimRange>,
          class VectorType = FieldVector<typename AnalyticalFluxType::RangeFieldType, AnalyticalFluxType::dimRange>>
class JacobianWrapper : public JacobianWrapperBase<AnalyticalFluxType, MatrixType, VectorType>
{
  using BaseType = JacobianWrapperBase<AnalyticalFluxType, MatrixType, VectorType>;

protected:
  using V = XT::Common::VectorAbstraction<VectorType>;
  using M = XT::Common::MatrixAbstraction<MatrixType>;
  using EigenSolverType = typename XT::LA::EigenSolver<MatrixType>;
  using typename BaseType::RangeFieldType;

public:
  static constexpr size_t dimDomain = AnalyticalFluxType::dimDomain;
  static constexpr size_t dimRange = AnalyticalFluxType::dimRange;
  using typename BaseType::JacobianType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::StateRangeType;

  using BaseType::jacobian;

  JacobianWrapper()
    : work_(1)
    , scale_(dimRange)
    , rconde_(dimRange)
    , rcondv_(dimRange)
    , iwork_(2 * dimRange - 2)
    , eigenvalues_(std::vector<RangeFieldType>(dimRange))
    , tau_(V::create(dimRange))
  {
    jacobian() = std::make_unique<JacobianType>(eigenvectors_);
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
                                                  M::data(jacobian(0)),
                                                  static_cast<int>(dimRange),
                                                  &(eigenvalues_[0][0]),
                                                  &(dummy_complex_eigenvalues_[0]),
                                                  nullptr,
                                                  static_cast<int>(dimRange),
                                                  M::data(eigenvectors_[0]),
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

  virtual void get_jacobian(const size_t dd,
                            const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const StateRangeType& u,
                            const XT::Common::Parameter& param) override final
  {
    return analytical_flux.local_function(entity)->partial_u_col(dd, x_in_inside_coords, u, jacobian(dd), param);
  }

  virtual void get_jacobian(const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const StateRangeType& u,
                            const XT::Common::Parameter& param) override final
  {
    analytical_flux.local_function(entity)->partial_u(x_in_inside_coords, u, *jacobian(), param);
  }

  using BaseType::compute;

  virtual void compute(const size_t dd) override
  {
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
                                                  M::data(jacobian(dd)),
                                                  static_cast<int>(dimRange),
                                                  eigenvalues_[dd].data(),
                                                  &(dummy_complex_eigenvalues_[0]),
                                                  nullptr,
                                                  static_cast<int>(dimRange),
                                                  M::data(eigenvectors_[dd]),
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
        DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
#endif // HAVE_MKL || HAVE_LAPACKE
    } else {
      static auto eigensolver_options = hyperbolic_default_eigensolver_options<MatrixType>();
      const auto eigensolver = EigenSolverType(jacobian(dd), &eigensolver_options);
      eigenvectors_[dd] = eigensolver.real_eigenvectors();
      eigenvalues_[dd] = eigensolver.real_eigenvalues();
    }
    QR_[dd] = eigenvectors_[dd];
    XT::LA::qr(QR_[dd], tau_[dd], permutations_[dd]);
    computed_[dd] = true;
  }

  virtual void apply_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const override final
  {
    eigenvectors_[dd].mv(x, ret);
  }

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const override final
  {
    thread_local VectorType work = V::create(dimRange);
    XT::LA::solve_qr_factorized(QR_[dd], tau_[dd], permutations_[dd], ret, x, &work);
  }

  virtual const MatrixType& eigenvectors(const size_t dd) const override final
  {
    return eigenvectors_[dd];
  }

protected:
  std::vector<RangeFieldType> work_, scale_, rconde_, rcondv_;
  std::vector<int> iwork_;
  using BaseType::computed_;
  JacobianType eigenvectors_;
  FieldVector<std::vector<RangeFieldType>, dimDomain> eigenvalues_;
  FieldVector<RangeFieldType, dimRange> dummy_complex_eigenvalues_;
  JacobianType QR_;
  FieldVector<VectorType, dimDomain> tau_;
  FieldVector<FieldVector<int, dimRange>, dimDomain> permutations_;
}; // class JacobianWrapper<...>


template <class AnalyticalFluxType, size_t block_size = (AnalyticalFluxType::dimDomain == 1) ? 2 : 4>
class BlockedJacobianWrapper
    : public JacobianWrapperBase<AnalyticalFluxType,
                                 XT::Common::BlockedFieldMatrix<typename AnalyticalFluxType::RangeFieldType,
                                                                AnalyticalFluxType::dimRange / block_size,
                                                                block_size>,
                                 XT::Common::BlockedFieldVector<typename AnalyticalFluxType::RangeFieldType,
                                                                AnalyticalFluxType::dimRange / block_size,
                                                                block_size>>
{
  using BaseType = JacobianWrapperBase<AnalyticalFluxType,
                                       XT::Common::BlockedFieldMatrix<typename AnalyticalFluxType::RangeFieldType,
                                                                      AnalyticalFluxType::dimRange / block_size,
                                                                      block_size>,
                                       XT::Common::BlockedFieldVector<typename AnalyticalFluxType::RangeFieldType,
                                                                      AnalyticalFluxType::dimRange / block_size,
                                                                      block_size>>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::RangeFieldType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::JacobianType;
  using typename BaseType::StateRangeType;
  static constexpr size_t num_blocks = VectorType::num_blocks;
  static_assert(dimRange % block_size == 0, "dimRange has to be a multiple of block_size");
  using BlockedIntVectorType = XT::Common::BlockedFieldVector<int, num_blocks, block_size>;
  using LocalVectorType = typename VectorType::BlockType;
  using LocalMatrixType = typename MatrixType::BlockType;
  using EigenSolverType = typename XT::LA::EigenSolver<LocalMatrixType>;
  using LocalM = typename XT::Common::MatrixAbstraction<LocalMatrixType>;
  using LocalV = typename XT::Common::VectorAbstraction<LocalVectorType>;

  BlockedJacobianWrapper()
  {
    std::fill_n(&(eigenvalues_[0][0]), dimDomain * num_blocks, std::vector<double>(block_size, 0.));
    jacobian() = std::make_unique<JacobianType>();
    nonblocked_jacobians_ = std::make_unique<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>>();
  }

  using BaseType::jacobian;
  using BaseType::compute;

  virtual void compute(const size_t dd) override final
  {
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      if (block_size == 2) {
        const auto& jac = jacobian(dd).block(jj);
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
        XT::LA::internal::fmatrix_compute_real_eigenvalues_and_real_right_eigenvectors_using_qr(
            jacobian(dd).block(jj), eigenvalues_[dd][jj], eigenvectors_[dd].block(jj));
      } // else (block_size == 2)
      QR_[dd].block(jj) = eigenvectors_[dd].block(jj);
      XT::LA::qr(QR_[dd].block(jj), tau_[dd].block(jj), permutations_[dd].block(jj));
    } // jj
    computed_[dd] = true;
  }

  virtual void apply_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const override final
  {
    eigenvectors_[dd].mv(x, ret);
  }

  virtual void get_jacobian(const size_t dd,
                            const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const VectorType& u,
                            const XT::Common::Parameter& param) override final
  {
    const auto local_func = analytical_flux.local_function(entity);
    local_func->partial_u_col(dd, x_in_inside_coords, u, (*nonblocked_jacobians_)[dd], param);
    jacobian(dd) = (*nonblocked_jacobians_)[dd];
  }

  virtual void get_jacobian(const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const VectorType& u,
                            const XT::Common::Parameter& param) override final
  {
    const auto local_func = analytical_flux.local_function(entity);
    local_func->partial_u(x_in_inside_coords, u, *nonblocked_jacobians_, param);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      jacobian(dd) = (*nonblocked_jacobians_)[dd];
  }

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const override final
  {
    LocalVectorType work;
    for (size_t jj = 0; jj < num_blocks; ++jj)
      XT::LA::solve_qr_factorized(
          QR_[dd].block(jj), tau_[dd].block(jj), permutations_[dd].block(jj), ret.block(jj), x.block(jj), &work);
  }

  virtual const MatrixType& eigenvectors(const size_t dd) const override final
  {
    return eigenvectors_[dd];
  }

protected:
  std::unique_ptr<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>> nonblocked_jacobians_;
  using BaseType::computed_;
  FieldVector<FieldVector<std::vector<double>, num_blocks>, dimDomain> eigenvalues_;
  FieldVector<MatrixType, dimDomain> eigenvectors_;
  FieldVector<MatrixType, dimDomain> QR_;
  FieldVector<VectorType, dimDomain> tau_;
  FieldVector<BlockedIntVectorType, dimDomain> permutations_;
}; // BlockedJacobianWrapper<...>


// Helper to iterate over all possible multiindices in a boost::multi_array
template <class MultiArrayType>
class MultiIndexIterator
{
public:
  static constexpr size_t dim = MultiArrayType::dimensionality;
  static constexpr size_t array_dim = dim == 0 ? 1 : dim;
  using IndicesType = std::array<size_t, array_dim>;

  MultiIndexIterator(const MultiArrayType& multi_array, const IndicesType& indices)
    : multi_array_(multi_array)
    , indices_(indices)
  {
  }

  IndicesType& operator*()
  {
    return indices_;
  }

  const IndicesType& operator*() const
  {
    return indices_;
  }

  MultiIndexIterator& operator++()
  {
    size_t ii = 0;
    for (; ii < array_dim; ++ii) {
      indices_[ii]++;
      if (indices_[ii] < (dim == 0 ? 1 : multi_array_.shape()[ii]))
        break;
      indices_[ii] = 0;
    }
    if (ii == array_dim) {
      if (dim == 0)
        indices_[0] = 1;
      else
        std::copy_n(multi_array_.shape(), array_dim, indices_.begin());
    }
    return *this;
  } // ... operator++()

  MultiIndexIterator operator++(int)
  {
    MultiIndexIterator ret = *this;
    this->operator++();
    return ret;
  } // ... operator++(int)

  bool operator==(const MultiIndexIterator& other) const
  {
    return indices_ == other.indices_;
  }

  bool operator!=(const MultiIndexIterator& other) const
  {
    return indices_ != other.indices_;
  }

private:
  const MultiArrayType& multi_array_;
  IndicesType indices_;
}; // class MultiIndexIterator<...>

template <class MultiArrayType>
class MultiIndexProvider
{
public:
  using IteratorType = MultiIndexIterator<MultiArrayType>;
  using IndicesType = typename IteratorType::IndicesType;

  MultiIndexProvider(const MultiArrayType& multi_array)
    : multi_array_(multi_array)
  {
  }

  IteratorType begin()
  {
    static const IndicesType zero_indices = []() {
      IndicesType ret;
      ret.fill(0);
      return ret;
    }();
    return IteratorType(multi_array_, zero_indices);
  }

  IteratorType end()
  {
    IndicesType indices;
    if (MultiArrayType::dimensionality == 0)
      indices[0] = 1; // for dimension 0, the end iterator has index 1
    else
      std::copy_n(multi_array_.shape(), MultiArrayType::dimensionality, indices.begin());
    return IteratorType(multi_array_, indices);
  }

private:
  const MultiArrayType& multi_array_;
}; // class MultiIndexProvider<...>

// Helper functor to build indices.
template <typename RangeArrayType, size_t num_ranges, size_t dimension>
struct IndicesBuilder
{
  static boost::detail::multi_array::index_gen<num_ranges, dimension> build(const RangeArrayType& ranges)
  {
    return boost::detail::multi_array::index_gen<num_ranges, dimension>(
        IndicesBuilder<RangeArrayType, num_ranges - 1, dimension>::build(ranges), ranges[num_ranges - 1]);
  }
};

// Helper functor specialization to terminate recursion.
template <typename RangeArrayType, size_t dimension>
struct IndicesBuilder<RangeArrayType, 0, dimension>
{
  static boost::detail::multi_array::index_gen<0, dimension> build(const RangeArrayType& /*ranges*/)
  {
    return boost::detail::multi_array::index_gen<0, dimension>();
  }
};

template <class RangeType, size_t dimDomain>
class Slice : public boost::multi_array<RangeType, dimDomain - 1>
{
};

template <class RangeType>
class Slice<RangeType, 1>
{
public:
  template <class MultiIndexType>
  RangeType& operator()(const MultiIndexType&)
  {
    return value_;
  }

  template <class IndicesType>
  void resize(const IndicesType&)
  {
  }

private:
  RangeType value_;
};


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_INTERNAL_HH
