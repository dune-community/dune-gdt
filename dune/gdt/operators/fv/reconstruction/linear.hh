// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH

#include <boost/multi_array.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/lapacke.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "../quadrature.hh"
#include "reconstructed_function.hh"
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
  //        XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
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
                            const StateRangeType& u,
                            const XT::Common::Parameter& param) = 0;

  virtual void get_jacobian(const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const StateRangeType& u,
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
      work_.resize(work_[0]);
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
                                 FieldVector<FieldMatrix<typename AnalyticalFluxType::RangeFieldType,
                                                         block_size,
                                                         block_size>,
                                             AnalyticalFluxType::dimRange / block_size>,
                                 typename AnalyticalFluxType::StateRangeType>
{
  using BaseType =
      JacobianWrapperBase<AnalyticalFluxType,
                          FieldVector<FieldMatrix<typename AnalyticalFluxType::RangeFieldType, block_size, block_size>,
                                      AnalyticalFluxType::dimRange / block_size>,
                          typename AnalyticalFluxType::StateRangeType>;

public:
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  static constexpr size_t num_blocks = dimRange / block_size;
  static_assert(dimRange % block_size == 0, "dimRange has to be a multiple of block_size");
  using LocalRangeType = FieldVector<RangeFieldType, block_size>;
  using LocalMatrixType = FieldMatrix<RangeFieldType, block_size, block_size>;
  using MatrixType = FieldVector<LocalMatrixType, num_blocks>;
  using EigenSolverType = typename XT::LA::EigenSolver<LocalMatrixType>;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::JacobianType;
  using typename BaseType::StateRangeType;
  using typename BaseType::VectorType;

  BlockedJacobianWrapper()
    : work_(1)
  {
    jacobian() = std::make_unique<JacobianType>();
    nonblocked_jacobians_ = std::make_unique<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>>();
#if HAVE_MKL || HAVE_LAPACKE
    // get optimal working size in work[0] (requested by lwork = -1)
    int info = XT::Common::Lapacke::dgeev_work(XT::Common::Lapacke::row_major(),
                                               /*do_not_compute_left_eigenvectors: */ 'N',
                                               /*compute_right_eigenvectors: */ 'V',
                                               static_cast<int>(block_size),
                                               &(jacobian(0)[0][0][0]),
                                               static_cast<int>(block_size),
                                               &(eigenvalues_[0][0][0]),
                                               &(dummy_complex_eigenvalues_[0]),
                                               nullptr,
                                               static_cast<int>(block_size),
                                               &(jacobian(0)[0][0][0]),
                                               static_cast<int>(block_size),
                                               work_.data(),
                                               -1);
    if (info != 0)
      DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
    work_.resize(work_[0]);
#endif
  }

  using BaseType::jacobian;
  using BaseType::compute;

  virtual void compute(const size_t dd) override final
  {
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      if (block_size == 2) {
        const auto& jac = jacobian(dd)[jj];
        const auto trace = jac[0][0] + jac[1][1];
        const auto det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        const auto sqrt_val = std::sqrt(0.25 * trace * trace - det);
        const auto eigval1 = 0.5 * trace + sqrt_val;
        const auto eigval2 = 0.5 * trace - sqrt_val;
        if (std::abs(jac[1][0]) > std::abs(jac[0][1])) {
          if (XT::Common::FloatCmp::ne(jac[1][0], 0.)) {
            eigenvectors_[dd][jj][0][0] = eigval1 - jac[1][1];
            eigenvectors_[dd][jj][0][1] = eigval2 - jac[1][1];
            eigenvectors_[dd][jj][1][0] = eigenvectors_[dd][jj][1][1] = jac[1][0];
          } else {
            eigenvectors_[dd][jj][0][0] = eigenvectors_[dd][jj][1][1] = 1.;
            eigenvectors_[dd][jj][0][1] = eigenvectors_[dd][jj][1][0] = 0.;
          }
        } else {
          if (XT::Common::FloatCmp::ne(jac[0][1], 0.)) {
            eigenvectors_[dd][jj][1][0] = eigval1 - jac[0][0];
            eigenvectors_[dd][jj][1][1] = eigval2 - jac[0][0];
            eigenvectors_[dd][jj][0][0] = eigenvectors_[dd][jj][0][1] = jac[0][1];
          } else {
            eigenvectors_[dd][jj][0][0] = eigenvectors_[dd][jj][1][1] = 1.;
            eigenvectors_[dd][jj][0][1] = eigenvectors_[dd][jj][1][0] = 0.;
          }
        }
      } else {
#if HAVE_MKL || HAVE_LAPACKE
        int info = XT::Common::Lapacke::dgeev_work(XT::Common::Lapacke::row_major(),
                                                   /*do_not_compute_left_eigenvectors: */ 'N',
                                                   /*compute_right_eigenvectors: */ 'V',
                                                   static_cast<int>(block_size),
                                                   &(jacobian(dd)[jj][0][0]),
                                                   static_cast<int>(block_size),
                                                   &(eigenvalues_[jj][0][0]),
                                                   &(dummy_complex_eigenvalues_[0]),
                                                   nullptr,
                                                   static_cast<int>(block_size),
                                                   &(eigenvectors_[dd][jj][0][0]),
                                                   static_cast<int>(block_size),
                                                   work_.data(),
                                                   static_cast<int>(work_.size()));
        if (info != 0)
          DUNE_THROW(Dune::XT::LA::Exceptions::eigen_solver_failed, "The lapack backend reported '" << info << "'!");
#else // HAVE_MKL || HAVE_LAPACKE
        static XT::Common::Configuration eigensolver_options = hyperbolic_default_eigensolver_options<MatrixType>();
        const auto eigensolver = EigenSolverType(jacobian(dd)[jj], &eigensolver_options);
        eigenvectors_[dd][jj] = eigensolver.real_eigenvectors();
#endif // HAVE_MKL || HAVE_LAPACKE
      } // else (block_size == 2)
      QR_[dd][jj] = eigenvectors_[dd][jj];
      XT::LA::qr(QR_[dd][jj], tau_[dd][jj], permutations_[dd][jj]);
    } // jj
    computed_[dd] = true;
  }

  virtual void apply_eigenvectors(const size_t dd, const StateRangeType& x, StateRangeType& ret) const override final
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          ret[offset + ll] += eigenvectors_[dd][jj][ll][mm] * x[offset + mm];
    } // jj
  }

  virtual void get_jacobian(const size_t dd,
                            const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const StateRangeType& u,
                            const XT::Common::Parameter& param) override final
  {
    const auto local_func = analytical_flux.local_function(entity);
    local_func->partial_u_col(dd, x_in_inside_coords, u, (*nonblocked_jacobians_)[dd], param);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = jj * block_size;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          jacobian(dd)[jj][ll][mm] = (*nonblocked_jacobians_)[dd][offset + ll][offset + mm];
    } // jj
  }

  virtual void get_jacobian(const EntityType& entity,
                            const AnalyticalFluxType& analytical_flux,
                            const DomainType& x_in_inside_coords,
                            const StateRangeType& u,
                            const XT::Common::Parameter& param) override final
  {
    const auto local_func = analytical_flux.local_function(entity);
    local_func->partial_u(x_in_inside_coords, u, *nonblocked_jacobians_, param);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = jj * block_size;
        for (size_t ll = 0; ll < block_size; ++ll)
          for (size_t mm = 0; mm < block_size; ++mm)
            jacobian(dd)[jj][ll][mm] = (*nonblocked_jacobians_)[dd][offset + ll][offset + mm];
      } // jj
    } // dd
  }

  virtual void apply_inverse_eigenvectors(const size_t dd, const VectorType& x, VectorType& ret) const override final
  {
    LocalRangeType work;
    LocalRangeType tmp_ret, tmp_x;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (size_t ll = 0; ll < block_size; ++ll)
        tmp_x[ll] = x[jj * block_size + ll];
      XT::LA::solve_qr_factorized(QR_[dd][jj], tau_[dd][jj], permutations_[dd][jj], tmp_ret, tmp_x, &work);
      for (size_t ll = 0; ll < block_size; ++ll)
        ret[jj * block_size + ll] = tmp_ret[ll];
    }
  }

  virtual const MatrixType& eigenvectors(const size_t dd) const override final
  {
    return eigenvectors_[dd];
  }

protected:
  std::vector<RangeFieldType> work_;
  std::unique_ptr<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>> nonblocked_jacobians_;
  using BaseType::computed_;
  FieldVector<FieldVector<FieldVector<RangeFieldType, block_size>, num_blocks>, dimDomain> eigenvalues_;
  FieldVector<RangeFieldType, dimRange> dummy_complex_eigenvalues_;
  FieldVector<MatrixType, dimDomain> eigenvectors_;
  FieldVector<MatrixType, dimDomain> QR_;
  FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain> tau_;
  FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain> permutations_;
}; // BlockedJacobianWrapper<...>


} // namespace internal


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

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, class JacobianWrapperType>
class LocalLinearReconstructionOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t axis_size = 3;
  using EntityType = typename GridLayerType::template Codim<0>::Entity;
  using IndexSetType = typename GridLayerType::IndexSet;
  using DomainType = typename BoundaryValueType::DomainType;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using RangeType = typename BoundaryValueType::RangeType;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using Quadrature1dType = Dune::QuadratureRule<DomainFieldType, 1>;
  using IntersectionType = typename GridLayerType::Intersection;
  using IntersectionVectorType = FieldVector<IntersectionType, 2 * dimDomain>;
  using IntersectionLocalCoordType = typename IntersectionType::Geometry::LocalCoordinate;
  using AnalyticalFluxLocalfunctionType = typename AnalyticalFluxType::LocalfunctionType;
  using StateRangeType = typename AnalyticalFluxLocalfunctionType::StateRangeType;
  using BoundaryInfoType = typename XT::Grid::BoundaryInfo<IntersectionType>;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using StencilType = boost::multi_array<boost::optional<RangeType>, dimDomain>;
  using MultiArrayType = boost::multi_array<RangeType, dimDomain>;
  using SliceType = Slice<RangeType, dimDomain>;
  using CoordsType = std::array<size_t, dimDomain>;
  using MatrixType = typename JacobianWrapperType::MatrixType;
  using SlopeType = SlopeBase<RangeType, MatrixType>;

public:
  explicit LocalLinearReconstructionOperator(const std::vector<RangeType>& source_values,
                                             const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const SlopeType& slope,
                                             const GridLayerType& grid_layer,
                                             const XT::Common::Parameter& param,
                                             const Quadrature1dType& quadrature,
                                             ReconstructedFunctionType& reconstructed_function,
                                             XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , slope_(slope)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_(quadrature)
    , reconstructed_function_(reconstructed_function)
    , jacobian_wrapper_(jacobian_wrapper)
  {
    param_.set("boundary", {0.});
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    static const CoordsType stencil_sizes = []() {
      CoordsType ret;
      const auto ax_size = axis_size; // avoid linker error
      ret.fill(ax_size);
      return ret;
    }();
    thread_local StencilType stencil(stencil_sizes);
    bool valid = fill_stencil(stencil, entity);
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!valid)
      return;
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    const auto entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_function_.values()[entity_index];

    // get jacobian
    auto& jac = *jacobian_wrapper_;
    if (!jac.computed() || !analytical_flux_.is_affine()) {
      const auto& u_entity = source_values_[entity_index];
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      jac.get_jacobian(entity, analytical_flux_, x_in_inside_coords, u_entity, param_);
      jac.compute();
      if (analytical_flux_.is_affine())
        jac.jacobian() = nullptr;
    }

    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (quadrature_.size() == 1) {
        // no need to reconstruct in all directions, as we are only regarding the center of the face, which will always
        // have the same value assigned, independent of the slope in the other directions
        std::array<size_t, dimDomain> indices;
        indices.fill(1);
        FieldVector<RangeType, axis_size> stencil_1d, stencil_1d_char;
        FieldVector<RangeType, 2> reconstructed_values;
        for (size_t ii = 0; ii < axis_size; ++ii) { // transform to characteristic variables
          indices[dd] = ii;
          stencil_1d[ii] = *stencil(indices);
          jac.apply_inverse_eigenvectors(dd, stencil_1d[ii], stencil_1d_char[ii]);
        }
        // perform the actual reconstruction
        linear_reconstruction_1d(stencil_1d, stencil_1d_char, reconstructed_values, jac.eigenvectors(dd));

        // convert back to non-characteristic variables
        auto tmp_value = reconstructed_values[0];
        jac.apply_eigenvectors(dd, tmp_value, reconstructed_values[0]);
        tmp_value = reconstructed_values[1];
        jac.apply_eigenvectors(dd, tmp_value, reconstructed_values[1]);

        // store reconstructed values
        reconstructed_values_map.emplace(intersections[2 * dd].geometryInInside().center(), reconstructed_values[0]);
        reconstructed_values_map.emplace(intersections[2 * dd + 1].geometryInInside().center(),
                                         reconstructed_values[1]);
      } else {
        thread_local MultiArrayType reconstructed_values(stencil_sizes);
        thread_local auto tmp_multiarray = reconstructed_values;
        tmp_multiarray.resize(stencil_sizes);

        // Transform values on stencil to characteristic variables of the current coordinate direction dd
        for (size_t ii = 0; ii < stencil.num_elements(); ++ii)
          jac.apply_inverse_eigenvectors(dd, *stencil.data()[ii], tmp_multiarray.data()[ii]);

        size_t curr_dir = dd;
        size_t last_dir = curr_dir;
        RangeType tmp_value;
        CoordsType current_sizes;
        for (size_t dir = 0; dir < dimDomain; ++dir) {
          curr_dir = (dd + dir) % dimDomain;
          // Transform to characteristic variables of the current reconstruction direction.
          if (dir > 0) {
            std::copy_n(reconstructed_values.shape(), dimDomain, current_sizes.begin());
            tmp_multiarray.resize(current_sizes);
            tmp_multiarray = reconstructed_values;
            std::for_each(
                tmp_multiarray.data(), tmp_multiarray.data() + tmp_multiarray.num_elements(), [&](RangeType& value) {
                  jac.apply_eigenvectors(last_dir, value, tmp_value);
                  jac.apply_inverse_eigenvectors(curr_dir, tmp_value, value);
                });
          } // if (dir > 0)
          // perform the actual reconstruction
          const auto& curr_quadrature = dir > 0 ? quadrature_ : get_left_right_quadrature();
          linear_reconstruction(
              curr_dir, curr_quadrature, tmp_multiarray, reconstructed_values, jac.eigenvectors(curr_dir));
          last_dir = curr_dir;
        } // dir
        // convert back to non-characteristic variables
        std::for_each(reconstructed_values.data(),
                      reconstructed_values.data() + reconstructed_values.num_elements(),
                      [&](RangeType& value) {
                        tmp_value = value;
                        jac.apply_eigenvectors(last_dir, tmp_value, value);
                      });

        // Convert coordinates on face to local entity coordinates and store reconstructed values
        MultiIndexProvider<MultiArrayType> multi_indices(reconstructed_values);
        for (const auto& multi_index : multi_indices) {
          IntersectionLocalCoordType quadrature_point;
          for (size_t ii = 0; ii < dimDomain; ++ii)
            if (ii != dd)
              quadrature_point[ii < dd ? ii : ii - 1] = quadrature_[multi_index[ii]].position();
          reconstructed_values_map.emplace(
              intersections[2 * dd + multi_index[dd]].geometryInInside().global(quadrature_point),
              reconstructed_values(multi_index));
        } // multi_indices
      }
    } // dd
  } // void apply_local(...)

  // stencil is in original coordinates, stencil_char in characteristic coordinates
  // returned reconstructed values are in characteristic coordinates
  void linear_reconstruction_1d(const FieldVector<RangeType, axis_size>& stencil,
                                const FieldVector<RangeType, axis_size>& stencil_char,
                                FieldVector<RangeType, 2>& reconstructed_values,
                                const MatrixType& eigenvectors)
  {
    const auto& u_entity_char = stencil_char[1];
    const auto slope_char = slope_.get(stencil, stencil_char, eigenvectors) * 0.5;
    reconstructed_values[0] = u_entity_char - slope_char;
    reconstructed_values[1] = u_entity_char + slope_char;
  } // void linear_reconstruction_1d(...)

  // This is only used in several dimensions if we use an interface quadrature different from the midpoint quadrature
  // Todo: adapt to new slopes
  void linear_reconstruction(const size_t dir,
                             const Quadrature1dType& quadrature,
                             const MultiArrayType& stencil,
                             MultiArrayType& reconstructed_values,
                             const MatrixType& /*eigenvectors*/)
  {
    DUNE_THROW(Dune::NotImplemented, "This needs to be adapted to the new slopes");
    // resize the reconstructed_values array
    CoordsType new_shape;
    std::copy_n(stencil.shape(), dimDomain, new_shape.begin());
    new_shape[dir] = quadrature.size();
    reconstructed_values.resize(new_shape);

    // get array_views corresponding to the left, center and right values
    using IndexRangeType = typename MultiArrayType::index_range;
    using RangeArrayType = XT::Common::FieldVector<IndexRangeType, dimDomain>;
    using IndicesBuilderType = IndicesBuilder<RangeArrayType, dimDomain, dimDomain - 1>;

    assert(stencil.shape()[dir] == 3);
    RangeArrayType ranges;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, stencil.shape()[ii]);
    ranges[dir] = IndexRangeType(0);
    const auto u_left = stencil[IndicesBuilderType::build(ranges)];
    ranges[dir] = IndexRangeType(1);
    const auto u_entity = stencil[IndicesBuilderType::build(ranges)];
    ranges[dir] = IndexRangeType(2);
    const auto u_right = stencil[IndicesBuilderType::build(ranges)];

    // calculate slopes
    thread_local SliceType slope;
    std::array<size_t, dimDomain - 1> slope_sizes;
    std::copy_n(u_entity.shape(), dimDomain - 1, slope_sizes.begin());
    slope.resize(slope_sizes);
    MultiIndexProvider<decltype(u_left)> multi_indices(u_left);
    // Todo: adapt to new slopes
    //    for (const auto& multi_index : multi_indices) {
    // slope(multi_index) = slope_limiter_.get(u_left(multi_index), u_entity(multi_index),
    // u_right(multi_index), eigenvectors);
    //    }

    // calculate reconstructed values
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ranges[ii] = IndexRangeType(0, new_shape[ii]);
    for (size_t ii = 0; ii < quadrature.size(); ++ii) {
      ranges[dir] = IndexRangeType(ii);
      auto reconstructed_ii = reconstructed_values[IndicesBuilderType::build(ranges)];
      for (const auto& multi_index : multi_indices) {
        reconstructed_ii(multi_index) = slope(multi_index);
        reconstructed_ii(multi_index) *= (quadrature[ii].position() - 0.5);
        reconstructed_ii(multi_index) += u_entity(multi_index);
      }
    } // ii (quadrature.size())
  } // void linear_reconstruction(...)

private:
  CoordsType center(const StencilType& stencil) const
  {
    CoordsType ret;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      assert(stencil.shape()[ii] % 2 && "Center not well-defined if one of the axis_sizes is even!");
      ret[ii] = stencil.shape()[ii] / 2;
    }
    return ret;
  }

  bool fill_stencil(StencilType& stencil, const EntityType& entity)
  {
    const int dir = -2;
    auto coords = center(stencil);
    std::fill_n(stencil.data(), stencil.num_elements(), boost::none);
    return fill_impl(stencil, entity, dir, coords);
  } // void fill(...)

private:
  bool fill_impl(StencilType& stencil, const EntityType& entity, const int dir, const CoordsType& coords)
  {
    bool ret = true;
    const auto entity_index = grid_layer_.indexSet().index(entity);
    stencil(coords) = source_values_[entity_index];
    std::vector<int> boundary_dirs;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const auto new_dir = intersection.indexInInside();
      if (direction_allowed(dir, new_dir) && !end_of_stencil(stencil, new_dir, coords)) {
        auto new_coords = coords;
        if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
          boundary_dirs.push_back(new_dir);
          auto boundary_value = boundary_values_.local_function(entity)->evaluate(
              intersection, entity.geometry().local(intersection.geometry().center()), source_values_[entity_index]);
          while (!end_of_stencil(stencil, new_dir, new_coords)) {
            next_coords_in_dir(new_dir, new_coords);
            stencil(new_coords) = boundary_value;
          }
        } else if (intersection.neighbor()) { // inner and periodic intersections
          const auto& outside = intersection.outside();
          next_coords_in_dir(new_dir, new_coords);
          ret = ret && fill_impl(stencil, outside, new_dir, new_coords);
        } else if (!intersection.neighbor() && !intersection.boundary()) { // processor boundary
          return false;
        }
      } // if (!end_of_stencil(...))
    } // intersections

    assert(boundary_dirs.size() <= dimDomain);
    if (boundary_dirs.size() > 1) {
      auto new_coords = coords;
      next_coords_in_dir(boundary_dirs[0], new_coords);
      const auto& boundary_value = stencil(new_coords);

      std::for_each(stencil.data(),
                    stencil.data() + stencil.num_elements(),
                    [&boundary_value](boost::optional<RangeType>& value) {
                      if (!value)
                        value = boundary_value;
                    });
    } // if (boundary_dirs.size() > 1)
    return ret;
  }

  //  get next coords in direction dir (increase or decrease coords in that direction)
  static void next_coords_in_dir(const int dir, CoordsType& coords)
  {
    dir % 2 ? coords[dir / 2]++ : coords[dir / 2]--;
  }

  // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
  // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
  // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
  // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
  // without emitting new iterators).
  static bool direction_allowed(const int dir, const int new_dir)
  {
    return new_dir == dir || new_dir / 2 > dir / 2;
  }

  bool end_of_stencil(const StencilType& stencil, const int dir, const CoordsType& coords)
  {
    return coords[dir / 2] == stencil.shape()[dir / 2] - 1 || coords[dir / 2] == 0;
  }


  // quadrature rule containing left and right interface points
  static const Quadrature1dType& get_left_right_quadrature()
  {
    static const Quadrature1dType ret = left_right_quadrature();
    return ret;
  }

  static Quadrature1dType left_right_quadrature()
  {
    Quadrature1dType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SlopeType& slope_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const Quadrature1dType& quadrature_;
  ReconstructedFunctionType& reconstructed_function_;
  XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper_;
}; // class LocalLinearReconstructionOperator


template <class AnalyticalFluxImp, class BoundaryValueImp, class JacobianWrapperImp, class Traits>
class LinearReconstructionOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, class JacobianWrapperImp>
struct LinearReconstructionOperatorTraits
{
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  using JacobianWrapperType = JacobianWrapperImp;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using RangeFieldType = typename BoundaryValueType::DomainFieldType;
  using FieldType = DomainFieldType;
  using JacobianType = NoJacobian;
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  using ProductQuadratureType = QuadratureRule<DomainFieldType, dimDomain - 1>;
  using Quadrature1dType = Dune::QuadratureRule<DomainFieldType, 1>;
  using derived_type = LinearReconstructionOperator<AnalyticalFluxType,
                                                    BoundaryValueType,
                                                    JacobianWrapperType,
                                                    LinearReconstructionOperatorTraits>;
};


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class JacobianWrapperImp = internal::JacobianWrapper<AnalyticalFluxImp,
                                                               FieldMatrix<typename BoundaryValueImp::RangeFieldType,
                                                                           BoundaryValueImp::dimRange,
                                                                           BoundaryValueImp::dimRange>,
                                                               FieldVector<typename BoundaryValueImp::RangeFieldType,
                                                                           BoundaryValueImp::dimRange>>,
          class Traits =
              internal::LinearReconstructionOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, JacobianWrapperImp>>
class LinearReconstructionOperator : public OperatorInterface<Traits>
{
public:
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BoundaryValueType = typename Traits::BoundaryValueType;
  using JacobianWrapperType = typename Traits::JacobianWrapperType;
  using Quadrature1dType = typename Traits::Quadrature1dType;
  using ProductQuadratureType = typename Traits::ProductQuadratureType;
  using DomainFieldType = typename Traits::DomainFieldType;
  using RangeFieldType = typename Traits::RangeFieldType;
  using MatrixType = typename JacobianWrapperType::MatrixType;
  using SlopeType = SlopeBase<typename BoundaryValueType::RangeType, MatrixType>;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static constexpr size_t dimRange = Traits::dimRange;

  LinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const SlopeType& slope,
                               const Quadrature1dType& quadrature_1d = default_1d_quadrature<DomainFieldType>(1))
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , slope_(slope)
    , quadrature_1d_(quadrature_1d)
    , product_quadrature_(product_quadrature_on_intersection<DomainFieldType, dimDomain>(quadrature_1d))
  {
  }

  const ProductQuadratureType& quadrature() const
  {
    return product_quadrature_;
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    static_assert(is_reconstructed_localizable_function<RangeType>::value,
                  "RangeType has to be derived from ReconstructedLocalizableFunction!");
    // evaluate cell averages
    const auto& grid_layer = source.space().grid_layer();
    const auto& index_set = grid_layer.indexSet();
    std::vector<typename SourceType::RangeType> source_values(index_set.size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto entity_index = index_set.index(entity);
      const auto local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator =
        LocalLinearReconstructionOperator<AnalyticalFluxType,
                                          BoundaryValueType,
                                          typename SourceType::SpaceType::GridLayerType,
                                          JacobianWrapperType>(source_values,
                                                               analytical_flux_,
                                                               boundary_values_,
                                                               slope_,
                                                               grid_layer,
                                                               param,
                                                               quadrature_1d_,
                                                               range,
                                                               jacobian_wrapper_);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SlopeType& slope_;
  const Quadrature1dType quadrature_1d_;
  const ProductQuadratureType product_quadrature_;
  mutable XT::Common::PerThreadValue<JacobianWrapperType> jacobian_wrapper_;
}; // class LinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
