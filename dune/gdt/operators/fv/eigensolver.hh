// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
#define DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH

#include <config.h>

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

namespace Dune {
namespace GDT {


// TODO: remove Eigen-dependency and use generic eigenvalue solver
#if HAVE_EIGEN

template <class LocalFluxFunctionImp, bool selfadjoint = false>
class EigenSolver
{
  static const size_t stateDimRange = LocalFluxFunctionImp::StateType::dimRange;
  static const size_t dimRange = LocalFluxFunctionImp::dimRange;
  static const size_t dimRangeCols = LocalFluxFunctionImp::dimRangeCols;
  static_assert(stateDimRange == dimRange,
                "Jacobian has to be a square matrix, i.e. stateDimRange and dimRange have to be equal!");
  typedef typename LocalFluxFunctionImp::RangeFieldType RangeFieldType;
  typedef typename LocalFluxFunctionImp::DomainType DomainType;
  typedef typename LocalFluxFunctionImp::StateRangeType StateRangeType;
  typedef typename XT::LA::EigenDenseVector<RangeFieldType> EigenVectorType;
  typedef typename XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef FieldVector<EigenVectorType, dimRangeCols> EigenValuesType;
  typedef FieldVector<EigenMatrixType, dimRangeCols> EigenVectorsType;
  typedef typename EigenMatrixType::BackendType EigenMatrixBackendType;
  typedef typename std::conditional<selfadjoint,
                                    ::Eigen::SelfAdjointEigenSolver<EigenMatrixBackendType>,
                                    ::Eigen::EigenSolver<EigenMatrixBackendType>>::type EigenSolverType;

public:
  EigenSolver(const LocalFluxFunctionImp& local_flux_function,
              const DomainType& x_local,
              const StateRangeType& u,
              const XT::Common::Parameter& param)
  {
    const auto partial_u = XT::Functions::JacobianRangeTypeConverter<stateDimRange, dimRange, dimRangeCols>::convert(
        local_flux_function.partial_u(x_local, u, param));
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      const auto partial_u_eigen =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<EigenMatrixType, dimRange, stateDimRange>::convert(
              partial_u[ii]);
      EigenSolverType eigen_solver(partial_u_eigen.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto& eigenvalues_eigen = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto& eigenvectors_eigen = eigen_solver.eigenvectors(); // <- this should be an Eigen vector of std::complex
      assert(size_t(eigenvalues_eigen.size()) == dimRange);
      assert(XT::Common::FloatCmp::eq(EigenVectorType(eigenvalues_eigen.imag()), EigenVectorType(dimRange, 0.)));
      assert(XT::Common::FloatCmp::eq(EigenMatrixType(eigenvectors_eigen.imag()),
                                      EigenMatrixType(dimRange, dimRange, 0.)));
      eigenvalues_[ii] = EigenVectorType(eigenvalues_eigen.real());
      eigenvectors_[ii] = EigenMatrixType(eigenvectors_eigen.real());
    }
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType eigenvectors_inverse() const
  {
    EigenVectorsType ret;
    for (size_t ii = 0; ii < dimRangeCols; ++ii)
      ret[ii] = EigenMatrixType(eigenvectors_[ii].backend().inverse());
    return ret;
  }

private:
  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
}; // class EigenSolver<...>

#else // HAVE_EIGEN

template <class LocalFluxFunctionImp, bool selfadjoint = false>
class EigenSolver
{
  static_assert(AlwaysFalse<LocalFluxFunctionImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
