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

#ifndef DUNE_GDT_OPERATORS_FV_HH
#define DUNE_GDT_OPERATORS_FV_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/local/fluxes/laxfriedrichs.hh>
#include <dune/gdt/local/operators/fv.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/base.hh>

#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include "interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {


enum class NumericalFluxes
{
  godunov,
  godunov_with_reconstruction,
  laxfriedrichs,
  laxfriedrichs_with_reconstruction,
  local_laxfriedrichs,
  local_laxfriedrichs_with_reconstruction
};

// forwards
template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionLaxFriedrichsOperator;

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionGodunovOperator;

template <class RHSEvaluationImp>
class AdvectionRHSOperator;


namespace internal {


// TODO: add static assert once type of BoundaryValueFunctionImp is decided
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp>
class AdvectionTraitsBase
{
  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  //  static_assert(Stuff::is_???< BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef typename AnalyticalFluxType::DomainFieldType FieldType;
}; // class AdvectionTraitsBase


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter>
class AdvectionLaxFriedrichsOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxImp,
                                         BoundaryValueFunctionImp,
                                         LocalizableFunctionImp,
                                         slope_limiter>
      derived_type;

}; // class AdvectionLaxFriedrichsOperatorTraits

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionGodunovOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp>
{
public:
  typedef AdvectionGodunovOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
}; // class AdvectionGodunovOperatorTraits


template <class RHSEvaluationImp>
class AdvectionRHSOperatorTraits
{
public:
  typedef AdvectionRHSOperator<RHSEvaluationImp> derived_type;
  typedef RHSEvaluationImp RHSEvaluationType;
  typedef typename RHSEvaluationImp::DomainFieldType FieldType;
}; // class AdvectionRHSOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueFunctionImp,
          class SourceImp,
          class RangeImp>
class AdvectionLocalizableDefault
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp> BaseType;

  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  static_assert(is_local_numerical_coupling_flux<NumericalCouplingFluxImp>::value,
                "NumericalCouplingFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_boundary_flux<NumericalBoundaryFluxImp>::value,
                "NumericalBoundaryFluxImp has to be derived from LocalNumericalBoundaryFluxInterface!");
  //  static_assert(std::is_base_of< ???, BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
  static_assert(is_discrete_function<SourceImp>::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalCouplingFvOperator<NumericalCouplingFluxType> LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFvOperator<NumericalBoundaryFluxType> LocalBoundaryOperatorType;

  template <class... LocalOperatorArgTypes>
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const BoundaryValueFunctionType& boundary_values,
                              const SourceType& src,
                              RangeType& rng,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(rng.space().grid_view(), src, rng)
    , local_operator_(analytical_flux, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
    , local_boundary_operator_(
          analytical_flux, boundary_values, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
  {
    this->append(local_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->append(local_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GridViewType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridViewType>());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


template <class SourceImp, class RangeImp, class BoundaryValueFunctionImp, class MatrixImp, SlopeLimiters slope_limiter>
class LinearReconstructionLocalizable
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp> BaseType;
  typedef LinearReconstructionLocalizable<SourceImp, RangeImp, BoundaryValueFunctionImp, MatrixImp, slope_limiter>
      ThisType;

public:
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef MatrixImp MatrixType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalReconstructionFvOperator<MatrixType, BoundaryValueFunctionType, slope_limiter>
      LocalOperatorType;

  LinearReconstructionLocalizable(const SourceType& src,
                                  RangeType& rng,
                                  const MatrixType& eigenvectors,
                                  const MatrixType& eigenvectors_inverse,
                                  const BoundaryValueFunctionType& boundary_values)
    : BaseType(rng.space().grid_view(), src, rng)
    , local_operator_(eigenvectors, eigenvectors_inverse, boundary_values)
    , source_(src)
    , range_(rng)
  {
    this->append(local_operator_);
  }

private:
  const LocalOperatorType local_operator_;
  const SourceType& source_;
  RangeType& range_;
}; // class LinearReconstructionLocalizable


// TODO: remove eigen dependency of LocalGodunovNumericalCouplingFlux/LocalGodunovNumericalBoundaryFlux
#if HAVE_EIGEN

namespace internal {


template <size_t domainDim, size_t rangeDim, class MatrixType, class EigenMatrixType, class AnalyticalFluxType>
struct EigenvectorInitializer
{
  static void initialize(const AnalyticalFluxType& /*analytical_flux*/,
                         const bool /*is_linear*/,
                         const bool use_linear_reconstruction,
                         MatrixType& /*eigenvectors*/,
                         MatrixType& /*eigenvectors_inverse*/)
  {
    if (use_linear_reconstruction) {
      DUNE_THROW(Dune::NotImplemented, "Linear reconstruction is only implemented in 1D!");
    }
  }
}; // struct EigenvectorInitializer<...>

template <size_t rangeDim, class MatrixType, class EigenMatrixType, class AnalyticalFluxType>
struct EigenvectorInitializer<1, rangeDim, MatrixType, EigenMatrixType, AnalyticalFluxType>
{
  static void initialize(const AnalyticalFluxType& analytical_flux,
                         const bool is_linear,
                         const bool use_linear_reconstruction,
                         MatrixType& eigenvectors,
                         MatrixType& eigenvectors_inverse)
  {
    if (use_linear_reconstruction) {
      assert(is_linear && "Linear reconstruction is only implemented for linear analytical fluxes!");
      // calculate matrix of eigenvectors of A, where A is the jacobian of the linear analytical flux, i.e. u_t + A*u_x
      // = 0. As the analytical flux is linear, the jacobian A is constant, so it is enough to evaluate at 0.
      ::Eigen::EigenSolver<typename EigenMatrixType::BackendType> eigen_solver(
          Dune::XT::Common::from_string<EigenMatrixType>(
              Dune::XT::Common::to_string(analytical_flux.jacobian(typename AnalyticalFluxType::RangeType(0))))
              .backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigen_eigenvectors = eigen_solver.eigenvectors();
#ifndef NDEBUG
      for (size_t ii = 0; ii < rangeDim; ++ii)
        for (size_t jj = 0; jj < rangeDim; ++jj)
          assert(eigen_eigenvectors(ii, jj).imag() < 1e-15);
#endif
      eigenvectors = Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.real())));
      eigenvectors_inverse = Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.inverse().real())));
    }
  }
}; // struct EigenvectorInitializer<1, ...>

template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          SlopeLimiters slope_limiter>
struct AdvectionOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueFunctionType,
            class MatrixType,
            class SourceType,
            class RangeType,
            class... LocalOperatorArgTypes>
  static void apply(const AnalyticalFluxType& analytical_flux,
                    const BoundaryValueFunctionType& boundary_values,
                    const MatrixType& eigenvectors,
                    const MatrixType& eigenvectors_inverse,
                    const SourceType& source,
                    RangeType& range,
                    const double time,
                    const bool use_linear_reconstruction,
                    LocalOperatorArgTypes&&... local_operator_args)
  {
    const auto current_boundary_values = boundary_values.evaluate_at_time(time);
    if (use_linear_reconstruction) {
      typedef DunePdelabDgProductSpaceWrapper<typename SourceType::SpaceType::GridViewType,
                                              1, // polOrder
                                              RangeFieldType,
                                              dimRange,
                                              dimRangeCols>
          DGSpaceType;
      typedef DiscreteFunction<DGSpaceType, typename SourceType::VectorType> ReconstructedDiscreteFunctionType;
      const auto dg_space = Dune::XT::Common::make_unique<DGSpaceType>(range.space().grid_view());
      const auto reconstruction =
          Dune::XT::Common::make_unique<ReconstructedDiscreteFunctionType>(*dg_space, "reconstructed");
      LinearReconstructionLocalizable<SourceType,
                                      ReconstructedDiscreteFunctionType,
                                      typename BoundaryValueFunctionType::ExpressionFunctionType,
                                      MatrixType,
                                      slope_limiter>
          reconstruction_operator(
              source, *reconstruction, eigenvectors, eigenvectors_inverse, *current_boundary_values);
      reconstruction_operator.apply();
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::ExpressionFunctionType,
                                  ReconstructedDiscreteFunctionType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               *current_boundary_values,
                               *reconstruction,
                               range,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply();
    } else {
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::ExpressionFunctionType,
                                  SourceType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               *current_boundary_values,
                               source,
                               range,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply();
    }
  }
}; // struct AdvectionOperatorApplier


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter>
class AdvectionLaxFriedrichsOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                                                         BoundaryValueFunctionImp,
                                                                                         LocalizableFunctionImp,
                                                                                         slope_limiter>>
{
public:
  typedef internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                         BoundaryValueFunctionImp,
                                                         LocalizableFunctionImp,
                                                         slope_limiter>
      Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = AnalyticalFluxType::dimRangeCols;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;

private:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  typedef typename Dune::GDT::LocalLaxFriedrichsNumericalCouplingFlux<AnalyticalFluxType,
                                                                      LocalizableFunctionType,
                                                                      dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      LocalLaxFriedrichsDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                       typename BoundaryValueFunctionType::TimeIndependentFunctionType,
                                                       LocalizableFunctionType,
                                                       dimDomain>
          NumericalBoundaryFluxType;

  AdvectionLaxFriedrichsOperator(const AnalyticalFluxType& analytical_flux,
                                 const BoundaryValueFunctionType& boundary_values,
                                 const LocalizableFunctionType& dx,
                                 const double dt,
                                 const bool is_linear = false,
                                 const bool use_linear_reconstruction = false,
                                 const bool use_local_laxfriedrichs_flux = false,
                                 const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
    , use_local_laxfriedrichs_flux_(use_local_laxfriedrichs_flux)
    , entity_geometries_equal_(entity_geometries_equal)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, is_linear_, use_linear_reconstruction_, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             source,
                                                             range,
                                                             time,
                                                             use_linear_reconstruction_,
                                                             dx_,
                                                             dt_,
                                                             is_linear_,
                                                             use_local_laxfriedrichs_flux_,
                                                             entity_geometries_equal_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_linear_reconstruction_;
  const bool use_local_laxfriedrichs_flux_;
  const bool entity_geometries_equal_;
  MatrixType eigenvectors_;
  MatrixType eigenvectors_inverse_;
}; // class AdvectionLaxFriedrichsOperator


// TODO: 0 boundary by default, so no need to specify boundary conditions for periodic grid views
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionGodunovOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                                   BoundaryValueFunctionImp,
                                                                                   slope_limiter>>
{
public:
  typedef internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = AnalyticalFluxType::dimRangeCols;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;

private:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  typedef typename Dune::GDT::LocalGodunovNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      LocalGodunovNumericalBoundaryFlux<AnalyticalFluxType,
                                        typename BoundaryValueFunctionType::TimeIndependentFunctionType,
                                        dimDomain>
          NumericalBoundaryFluxType;

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueFunctionType& boundary_values,
                           const bool is_linear = false,
                           const bool use_linear_reconstruction = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, is_linear_, use_linear_reconstruction_, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             source,
                                                             range,
                                                             time,
                                                             use_linear_reconstruction_,
                                                             is_linear_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const bool is_linear_;
  const bool use_linear_reconstruction_;
  MatrixType eigenvectors_;
  MatrixType eigenvectors_inverse_;
}; // class AdvectionGodunovOperator

#else // HAVE_EIGEN

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter>
class AdvectionLaxFriedrichsOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionGodunovOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN


template <class RHSEvaluationImp>
class AdvectionRHSOperator : public Dune::GDT::OperatorInterface<internal::AdvectionRHSOperatorTraits<RHSEvaluationImp>>
{
  static_assert(is_rhs_evaluation<RHSEvaluationImp>::value, "RHSEvaluationImp has to be derived from RHSInterface!");

public:
  typedef internal::AdvectionRHSOperatorTraits<RHSEvaluationImp> Traits;
  typedef typename Traits::RHSEvaluationType RHSEvaluationType;
  typedef LocalRhsFvOperator<RHSEvaluationType> LocalOperatorType;

  AdvectionRHSOperator(const RHSEvaluationType& rhs_evaluation)
    : local_operator_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const double /*time*/ = 0.0) const
  {
    LocalizableOperatorBase<typename RangeType::SpaceType::GridViewType, SourceType, RangeType> localizable_operator(
        range.space().grid_view(), source, range);
    localizable_operator.append(local_operator_);
    localizable_operator.apply();
  }

private:
  const LocalOperatorType local_operator_;
}; // class AdvectionRHSOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_HH
