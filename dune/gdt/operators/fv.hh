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

#include "config.h"

#include <memory>
#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/base.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/local/fluxes/laxfriedrichs.hh>
#include <dune/gdt/local/fluxes/kinetic.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/fv.hh>
#include <dune/gdt/local/operators/fv.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include "interfaces.hh"
#include "base.hh"
#include "matrix_exponential.cpp"

namespace Dune {
namespace GDT {


enum class NumericalFluxes
{
  godunov,
  godunov_with_reconstruction,
  laxfriedrichs,
  laxfriedrichs_with_reconstruction,
  local_laxfriedrichs,
  local_laxfriedrichs_with_reconstruction,
  kinetic
};

// forwards
template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionLaxFriedrichsOperator;

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder = 2,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionLaxFriedrichsWENOOperator;

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder = 2,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionGodunovWENOOperator;


template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionGodunovOperator;

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod>
class AdvectionKineticOperator;

template <class RHSEvaluationImp>
class AdvectionRHSOperator;


namespace internal {


// TODO: add static assert once type of BoundaryValueFunctionImp is decided
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionTraitsBase
{
  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  //  static_assert(Stuff::is_???< BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
public:
  static const SlopeLimiters slope_limiter = slope_limiter_type;
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = AnalyticalFluxType::dimRangeCols;
  typedef typename AnalyticalFluxType::DomainFieldType FieldType;
  typedef typename AnalyticalFluxType::FluxJacobianRangeType JacobianType;
}; // class AdvectionTraitsBase


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_limiter_type>
class AdvectionLaxFriedrichsOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxType,
                                                                      LocalizableFunctionType,
                                                                      dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                       typename BoundaryValueFunctionType::NonparametricType,
                                                       LocalizableFunctionType,
                                                       dimDomain>
          NumericalBoundaryFluxType;
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxImp,
                                         BoundaryValueFunctionImp,
                                         LocalizableFunctionImp,
                                         slope_limiter>
      derived_type;
}; // class AdvectionLaxFriedrichsOperatorTraits

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder,
          SlopeLimiters slope_limiter_type>
class AdvectionLaxFriedrichsWENOOperatorTraits : public AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                                                             BoundaryValueFunctionImp,
                                                                                             LocalizableFunctionImp,
                                                                                             slope_limiter_type>
{
public:
  typedef AdvectionLaxFriedrichsWENOOperator<AnalyticalFluxImp,
                                             BoundaryValueFunctionImp,
                                             LocalizableFunctionImp,
                                             GridViewType,
                                             Polyhedron_3,
                                             polOrder,
                                             slope_limiter_type>
      derived_type;
}; // class AdvectionLaxFriedrichsWENOOperatorTraits

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionGodunovOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovLocalNumericalBoundaryFlux<AnalyticalFluxType,
                                                                typename BoundaryValueFunctionType::NonparametricType,
                                                                dimDomain>
      NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
}; // class AdvectionGodunovOperatorTraits

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder,
          SlopeLimiters slope_limiter_type>
class AdvectionGodunovWENOOperatorTraits
    : public AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
public:
  typedef AdvectionGodunovWENOOperator<AnalyticalFluxImp,
                                       BoundaryValueFunctionImp,
                                       GridViewType,
                                       Polyhedron_3,
                                       polOrder,
                                       slope_limiter_type>
      derived_type;
}; // class AdvectionGodunovWENOOperatorTraits


template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
class AdvectionKineticOperatorTraits
    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

public:
  using BaseType::slope_limiter;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::KineticLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::KineticLocalNumericalBoundaryFlux<AnalyticalFluxType,
                                                                typename BoundaryValueFunctionType::NonparametricType,
                                                                dimDomain>
      NumericalBoundaryFluxType;
  typedef AdvectionKineticOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
}; // class AdvectionKineticOperatorTraits


template <class RHSEvaluationImp>
class AdvectionRHSOperatorTraits
{
public:
  typedef AdvectionRHSOperator<RHSEvaluationImp> derived_type;
  typedef RHSEvaluationImp RHSEvaluationType;
  typedef typename RHSEvaluationImp::DomainFieldType FieldType;
  typedef typename RHSEvaluationImp::RhsJacobianRangeType JacobianType;
}; // class AdvectionRHSOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueFunctionImp,
          class SourceImp,
          class RangeImp>
class AdvectionLocalizableDefault
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp> BaseType;

  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  static_assert(is_local_numerical_coupling_flux<NumericalCouplingFluxImp>::value,
                "NumericalCouplingFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_boundary_flux<NumericalBoundaryFluxImp>::value,
                "NumericalBoundaryFluxImp has to be derived from LocalNumericalBoundaryFluxInterface!");
  //  static_assert(std::is_base_of< ???, BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
  //  static_assert(is_discrete_function<SourceImp>::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridLayerType GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename Dune::GDT::LocalCouplingFvOperator<NumericalCouplingFluxType> LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFvOperator<NumericalBoundaryFluxType> LocalBoundaryOperatorType;

  template <class QuadratureRuleType, class... LocalOperatorArgTypes>
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const std::shared_ptr<BoundaryValueFunctionType>& boundary_values,
                              const SourceType& source,
                              RangeType& range,
                              const QuadratureRuleType& quadrature_rule,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(range.space().grid_view(), source, range)
    , local_operator_(quadrature_rule, analytical_flux, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
    , local_boundary_operator_(quadrature_rule,
                               analytical_flux,
                               boundary_values,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...)
  {
    this->append(local_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridLayerType>());
    this->append(local_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GridLayerType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridLayerType>());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


template <class SourceImp, class RangeImp, class BoundaryValueFunctionImp, class MatrixImp, SlopeLimiters slope_limiter>
class LinearReconstructionLocalizable
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp> BaseType;
  typedef LinearReconstructionLocalizable<SourceImp, RangeImp, BoundaryValueFunctionImp, MatrixImp, slope_limiter>
      ThisType;

public:
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef MatrixImp MatrixType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridLayerType GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename Dune::GDT::LocalReconstructionFvOperator<MatrixType, BoundaryValueFunctionType, slope_limiter>
      LocalOperatorType;

  LinearReconstructionLocalizable(const SourceType& source,
                                  RangeType& range,
                                  const MatrixType& eigenvectors,
                                  const MatrixType& eigenvectors_inverse,
                                  const std::shared_ptr<BoundaryValueFunctionType>& boundary_values)
    : BaseType(range.space().grid_view(), source, range)
    , local_operator_(eigenvectors, eigenvectors_inverse, boundary_values)
    , source_(source)
    , range_(range)
  {
    this->append(local_operator_);
  }

private:
  const LocalOperatorType local_operator_;
  const SourceType& source_;
  RangeType& range_;
}; // class LinearReconstructionLocalizable


// TODO: remove eigen dependency of GodunovLocalNumericalCouplingFlux/GodunovLocalNumericalBoundaryFlux
#if HAVE_EIGEN

namespace internal {


template <size_t domainDim, size_t rangeDim, class MatrixType, class EigenMatrixType, class AnalyticalFluxType>
struct EigenvectorInitializer
{
  static void initialize(const AnalyticalFluxType& /*analytical_flux*/,
                         const bool /*flux_is_linear*/,
                         const bool use_linear_reconstruction,
                         std::shared_ptr<MatrixType>& /*eigenvectors*/,
                         std::shared_ptr<MatrixType>& /*eigenvectors_inverse*/)
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
                         const bool flux_is_linear,
                         const bool use_linear_reconstruction,
                         std::shared_ptr<MatrixType>& eigenvectors,
                         std::shared_ptr<MatrixType>& eigenvectors_inverse)
  {
    if (use_linear_reconstruction) {
      assert(flux_is_linear && "Linear reconstruction is only implemented for linear analytical fluxes!");
      // calculate matrix of eigenvectors of A, where A is the jacobian of the linear analytical flux, i.e. u_t + A*u_x
      // = 0. As the analytical flux is linear, the jacobian A is constant, so it is enough to evaluate at 0.
      ::Eigen::EigenSolver<typename EigenMatrixType::BackendType> eigen_solver(
          Dune::XT::Common::from_string<EigenMatrixType>(
              Dune::XT::Common::to_string(
                  analytical_flux.jacobian(typename AnalyticalFluxType::RangeType(0),
                                           typename AnalyticalFluxType::EntityType{},
                                           typename AnalyticalFluxType::EntityType::Geometry::LocalCoordinate(0))))
              .backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigen_eigenvectors = eigen_solver.eigenvectors();
#ifndef NDEBUG
      for (size_t ii = 0; ii < rangeDim; ++ii)
        for (size_t jj = 0; jj < rangeDim; ++jj)
          assert(eigen_eigenvectors(ii, jj).imag() < 1e-15);
#endif
      eigenvectors = std::make_shared<MatrixType>(Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.real()))));
      eigenvectors_inverse = std::make_shared<MatrixType>(Dune::XT::Common::from_string<MatrixType>(
          Dune::XT::Common::to_string(EigenMatrixType(eigen_eigenvectors.inverse().real()))));
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
                    const SourceType& source,
                    RangeType& range,
                    const XT::Common::Parameter& param,
                    const bool use_linear_reconstruction,
                    const std::shared_ptr<MatrixType>& eigenvectors,
                    const std::shared_ptr<MatrixType>& eigenvectors_inverse,
                    LocalOperatorArgTypes&&... local_operator_args)
  {
    const std::shared_ptr<typename BoundaryValueFunctionType::NonparametricType> current_boundary_values =
        boundary_values.with_parameter(param);
    typedef typename SourceType::SpaceType::GridViewType GridViewType;
    const GridViewType& grid_view = source.space().grid_view();
    const auto quadrature_rule = Dune::QuadratureRules<RangeFieldType, GridViewType::dimension - 1>::rule(
        grid_view.ibegin(*(grid_view.template begin<0>()))->geometry().type(), 0);
    if (use_linear_reconstruction) {
      typedef DunePdelabDgProductSpaceWrapper<typename SourceType::SpaceType::GridLayerType,
                                              1, // polOrder
                                              RangeFieldType,
                                              dimRange,
                                              dimRangeCols>
          DGSpaceType;
      typedef DiscreteFunction<DGSpaceType, typename SourceType::VectorType> ReconstructedDiscreteFunctionType;
      const auto dg_space = Dune::XT::Common::make_unique<const DGSpaceType>(range.space().grid_layer());
      const auto reconstruction =
          Dune::XT::Common::make_unique<ReconstructedDiscreteFunctionType>(*dg_space, "reconstructed");
      LinearReconstructionLocalizable<SourceType,
                                      ReconstructedDiscreteFunctionType,
                                      typename BoundaryValueFunctionType::NonparametricType,
                                      MatrixType,
                                      slope_limiter>
          reconstruction_operator(
              source, *reconstruction, *eigenvectors, *eigenvectors_inverse, current_boundary_values);
      reconstruction_operator.apply();
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  ReconstructedDiscreteFunctionType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               *reconstruction,
                               range,
                               quadrature_rule,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply();
    } else {
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  SourceType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               source,
                               range,
                               quadrature_rule,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply(true);
    }
  }
}; // struct AdvectionOperatorApplier


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template <class GridViewImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class ReconstructedLocalizableFunction
    : public XT::Functions::LocalizableFunctionInterface<typename GridViewImp::template Codim<0>::Entity,
                                                         DomainFieldImp,
                                                         domainDim,
                                                         RangeFieldImp,
                                                         rangeDim,
                                                         rangeDimCols>
{
  typedef XT::Functions::LocalizableFunctionInterface<typename GridViewImp::template Codim<0>::Entity,
                                                      DomainFieldImp,
                                                      domainDim,
                                                      RangeFieldImp,
                                                      rangeDim,
                                                      rangeDimCols>
      BaseType;

public:
  static const constexpr size_t dimDomain = BaseType::dimDomain;
  static const constexpr size_t dimRange = BaseType::dimRange;
  static const constexpr size_t dimRangeCols = BaseType::dimRangeCols;

  typedef GridViewImp GridViewType;
  typedef typename GridViewType::IndexSet IndexSetType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename GridViewType::template Codim<0>::Geometry::LocalCoordinate LocalCoordinateType;

private:
  class ReconstructedLocalfunction
      : public XT::Functions::
            LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
  {
    typedef typename XT::Functions::
        LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
            BaseType;

  public:
    ReconstructedLocalfunction(const EntityType& entity,
                               const std::map<LocalCoordinateType, RangeType, internal::FieldVectorLess>& values)
      : BaseType(entity)
      , values_(values)
    {
    }

    virtual size_t order() const
    {
      DUNE_THROW(Dune::InvalidStateException, "This function can't be integrated!");
      return 2;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const
    {
      try {
        ret = values_.at(xx);
      } catch (const std::out_of_range& e) {
        DUNE_THROW(Dune::RangeError,
                   "There are no values for xx = " << XT::Common::to_string(xx) << " in this function!");
      }
    }

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const std::map<LocalCoordinateType, RangeType, internal::FieldVectorLess>& values_;
  };

public:
  typedef ReconstructedLocalfunction LocalfunctionType;
  typedef XT::Functions::
      LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
          LocalfunctionInterfaceType;

  static const bool available = true;

  ReconstructedLocalizableFunction(
      const GridViewType& grid_view,
      const std::vector<std::map<LocalCoordinateType, RangeType, internal::FieldVectorLess>>& reconstructed_values)
    : index_set_(grid_view.indexSet())
    , reconstructed_values_(reconstructed_values)
  {
  }

  virtual std::unique_ptr<LocalfunctionInterfaceType> local_function(const EntityType& entity) const
  {
    return std::make_unique<LocalfunctionType>(entity, reconstructed_values_[index_set_.index(entity)]);
  }

  static std::string static_id()
  {
    return "reconstructed localizable function";
  }

  virtual std::string type() const
  {
    return "reconstructed localizable function";
  }

  virtual std::string name() const
  {
    return "reconstructed localizable function";
  }

private:
  const IndexSetType& index_set_;
  const std::vector<std::map<LocalCoordinateType, RangeType, internal::FieldVectorLess>>& reconstructed_values_;
}; // class ReconstructedLocalizableFunction


template <class BoundaryValueType, class GridViewType, class SourceType, class VectorType, size_t dimDomain>
struct CellAveragesGetter
{
  static std::vector<std::vector<std::vector<VectorType>>>
  get_cell_averages(const BoundaryValueType& boundary_values,
                    const GridViewType& grid_view,
                    const FieldVector<size_t, dimDomain>& grid_sizes,
                    const std::vector<FieldVector<size_t, dimDomain>>& entity_indices,
                    const SourceType& source)
  {
    static const size_t precision = 15;
    std::vector<std::vector<std::vector<VectorType>>> cell_averages(
        grid_sizes[0] + 2,
        std::vector<std::vector<VectorType>>(grid_sizes[1] + 2, std::vector<VectorType>(grid_sizes[2] + 2)));
    for (const auto& entity : Dune::elements(grid_view)) {
      const auto& local_source = source.local_function(entity);
      const auto& indices = entity_indices[grid_view.indexSet().index(entity)];
      cell_averages[indices[0]][indices[1]][indices[2]] = XT::Common::from_string<VectorType>(XT::Common::to_string(
          local_source->evaluate(entity.geometry().local(entity.geometry().center())), precision));
      if (entity.hasBoundaryIntersections()) {
        const auto local_boundary_values = boundary_values.local_function(entity);
        for (const auto& intersection : Dune::intersections(grid_view, entity)) {
          if (intersection.boundary()) {
            const auto intersection_center = intersection.geometry().local(intersection.geometry().center());
            const auto& n = intersection.unitOuterNormal(intersection_center);
            size_t direction(0);
            bool right(false);
            for (size_t dd = 0; dd < dimDomain; ++dd) {
              if (XT::Common::FloatCmp::eq(n[dd], 1.)) {
                direction = dd;
                right = true;
                break;
              } else if (XT::Common::FloatCmp::eq(n[dd], -1.)) {
                direction = dd;
                right = false;
                break;
              }
            } // dd
            if (direction == 0)
              cell_averages[grid_sizes[0] + right][indices[1]][indices[2]] =
                  XT::Common::from_string<VectorType>(XT::Common::to_string(
                      local_boundary_values->evaluate(intersection.geometryInInside().global(intersection_center)),
                      precision));
            else if (direction == 1)
              cell_averages[indices[0]][grid_sizes[1] + right][indices[2]] =
                  XT::Common::from_string<VectorType>(XT::Common::to_string(
                      local_boundary_values->evaluate(intersection.geometryInInside().global(intersection_center)),
                      precision));
            else
              cell_averages[indices[0]][indices[1]][grid_sizes[2] + right] =
                  XT::Common::from_string<VectorType>(XT::Common::to_string(
                      local_boundary_values->evaluate(intersection.geometryInInside().global(intersection_center)),
                      precision));
          } // intersection.boundary()
        } // intersections
      } // entity.hasBoundaryIntersections()
    } // entities

    // TODO: proper handling of boundary conditions
    // for now, just copy everything from adjacent boundary cells
    for (size_t ii = 0; ii < 2; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        for (size_t kk = 0; kk < 2; ++kk)
          cell_averages[grid_sizes[0] + ii][grid_sizes[1] + jj][grid_sizes[2] + kk] =
              cell_averages[grid_sizes[0] + ii][jj == 0 ? 0 : grid_sizes[1] - 1][kk == 0 ? 0 : grid_sizes[2] - 1];

    for (size_t ii = 0; ii < 2; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        for (size_t kk = 0; kk < grid_sizes[2]; ++kk)
          cell_averages[grid_sizes[0] + ii][grid_sizes[1] + jj][kk] =
              cell_averages[grid_sizes[0] + ii][jj == 0 ? 0 : grid_sizes[1] - 1][kk];

    for (size_t ii = 0; ii < grid_sizes[0]; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        for (size_t kk = 0; kk < 2; ++kk)
          cell_averages[ii][grid_sizes[1] + jj][grid_sizes[2] + kk] =
              cell_averages[ii][grid_sizes[1] + jj][kk == 0 ? 0 : grid_sizes[2] - 1];

    for (size_t ii = 0; ii < 2; ++ii)
      for (size_t jj = 0; jj < grid_sizes[1]; ++jj)
        for (size_t kk = 0; kk < 2; ++kk)
          cell_averages[grid_sizes[0] + ii][jj][grid_sizes[2] + kk] =
              cell_averages[grid_sizes[0] + ii][jj][kk == 0 ? 0 : grid_sizes[2] - 1];

    return cell_averages;
  } // ... get_cell_averages(...)
}; // class CellAverageGetter<...>

template <class BoundaryValueType, class GridViewType, class SourceType, class VectorType>
struct CellAveragesGetter<BoundaryValueType, GridViewType, SourceType, VectorType, 1>
{
  static std::vector<VectorType> get_cell_averages(const BoundaryValueType& boundary_values,
                                                   const GridViewType& grid_view,
                                                   const FieldVector<size_t, 1>& grid_sizes,
                                                   const std::vector<FieldVector<size_t, 1>>& entity_indices,
                                                   const SourceType& source)
  {
    static const size_t precision = 15;
    std::vector<VectorType> cell_averages(grid_sizes[0] + 2);
    for (const auto& entity : Dune::elements(grid_view)) {
      const auto& local_source = source.local_function(entity);
      const auto& indices = entity_indices[grid_view.indexSet().index(entity)];
      const auto& value = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
      cell_averages[indices[0]] = XT::Common::from_string<VectorType>(XT::Common::to_string(value, 15));
      if (entity.hasBoundaryIntersections()) {
        const auto local_boundary_values = boundary_values.local_function(entity);
        for (const auto& intersection : Dune::intersections(grid_view, entity)) {
          if (intersection.boundary()) {
            const auto intersection_center = intersection.geometry().local(intersection.geometry().center());
            bool right =
                XT::Common::FloatCmp::eq(intersection.unitOuterNormal(intersection_center)[0], 1.) ? true : false;
            cell_averages[grid_sizes[0] + right] = XT::Common::from_string<VectorType>(XT::Common::to_string(
                local_boundary_values->evaluate(intersection.geometryInInside().global(intersection_center)),
                precision));
          } // intersection.boundary()
        } // intersections
      } // entity.hasBoundaryIntersections()
    } // entities
    return cell_averages;
  }
};

template <class FieldType, size_t dimDomain>
struct QuadratureRuleGetter
{
  static Dune::QuadratureRule<FieldType, dimDomain - 1>
  get_quadrature(const FieldVector<Dune::QuadratureRule<FieldType, 1>, dimDomain>& quadrature_rules)
  {
    Dune::QuadratureRule<FieldType, dimDomain - 1> ret;
    const auto& quadrature_rule = quadrature_rules[0];
    for (size_t ii = 0; ii < quadrature_rule.size(); ++ii)
      for (size_t jj = 0; jj < quadrature_rule.size(); ++jj)
        ret.push_back(Dune::QuadraturePoint<FieldType, dimDomain - 1>(
            {quadrature_rule[ii].position()[0], quadrature_rule[jj].position()[0]},
            quadrature_rule[ii].weight() * quadrature_rule[jj].weight()));
    return ret;
  }
};

template <class FieldType>
struct QuadratureRuleGetter<FieldType, 1>
{
  static Dune::QuadratureRule<FieldType, 0>
  get_quadrature(const FieldVector<Dune::QuadratureRule<FieldType, 1>, 1>& /*quadrature_rules*/)
  {
    Dune::QuadratureRule<FieldType, 0> ret;
    ret.push_back(Dune::QuadraturePoint<FieldType, 0>(FieldVector<FieldType, 0>(0), 1));
    return ret;
  }
};


template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          class RangeFieldType,
          class Polyhedron_3,
          size_t dimDomain,
          size_t dimRange,
          size_t dimRangeCols,
          size_t polOrder,
          SlopeLimiters slope_limiter>
struct AdvectionWENOOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueFunctionType,
            class SourceType,
            class RangeType,
            class... LocalOperatorArgTypes>
  static void apply(const AnalyticalFluxType& analytical_flux,
                    const BoundaryValueFunctionType& boundary_values,
                    const SourceType& source,
                    RangeType& range,
                    const XT::Common::Parameter& param,
                    const bool use_reconstruction,
                    const size_t quadrature_order,
                    const std::vector<FieldVector<size_t, dimDomain>>& entity_indices,
                    const FieldVector<size_t, dimDomain>& grid_sizes,
                    const Polyhedron_3& poly,
                    LocalOperatorArgTypes&&... local_operator_args)
  {
    typedef double FieldType;
    typedef XT::LA::EigenDenseVector<FieldType> EigenVectorType;
    const std::shared_ptr<typename BoundaryValueFunctionType::NonparametricType> current_boundary_values =
        boundary_values.with_parameter(param);
    if (use_reconstruction) {
      typedef typename SourceType::SpaceType::GridViewType GridViewType;
      const GridViewType& grid_view = source.space().grid_view();
      // collect cell averages in array
      const auto cell_averages = CellAveragesGetter<typename BoundaryValueFunctionType::NonparametricType,
                                                    typename SourceType::SpaceType::GridViewType,
                                                    SourceType,
                                                    EigenVectorType,
                                                    dimDomain>::get_cell_averages(*current_boundary_values,
                                                                                  grid_view,
                                                                                  grid_sizes,
                                                                                  entity_indices,
                                                                                  source);

      // get 1D quadrature rules
      const auto quadrature_rule = Dune::QuadratureRules<FieldType, 1>::rule(
          grid_view.begin<dimDomain - 1>()->geometry().type(), quadrature_order, QuadratureType::Enum::GaussLegendre);
      FieldVector<Dune::QuadratureRule<FieldType, 1>, dimDomain> quadrature_rules;
      std::fill(quadrature_rules.begin(), quadrature_rules.end(), quadrature_rule);

      // do reconstruction
      std::vector<std::map<typename GridViewType::template Codim<0>::Geometry::LocalCoordinate,
                           typename AnalyticalFluxType::RangeType,
                           internal::FieldVectorLess>>
          reconstructed_values(grid_view.size(0));

      auto local_reconstruction_operator = LocalWENOReconstructionFvOperator<GridViewType,
                                                                             AnalyticalFluxType,
                                                                             dimDomain,
                                                                             dimRange,
                                                                             polOrder,
                                                                             slope_limiter>(
          grid_view, analytical_flux, param, cell_averages, entity_indices, quadrature_rules, reconstructed_values);
      auto walker = XT::Grid::Walker<GridViewType>(grid_view);
      walker.append(local_reconstruction_operator);
      walker.walk(true);

      // do limiting for realizability in M_N models
      auto local_realizability_limiter = LocalRealizabilityLimiter<GridViewType, Polyhedron_3, dimDomain, dimRange>(
          grid_view, poly, cell_averages, entity_indices, reconstructed_values);
      walker.clear();
      walker.append(local_realizability_limiter);
      walker.walk(false);

      typedef ReconstructedLocalizableFunction<GridViewType, FieldType, dimDomain, RangeFieldType, dimRange>
          ReconstructedLocalizableFunctionType;
      const ReconstructedLocalizableFunctionType reconstructed_function(grid_view, reconstructed_values);

      // get quadrature rule on intersection from 1d quadratures;
      const auto intersection_quadrature_rule =
          QuadratureRuleGetter<FieldType, dimDomain>::get_quadrature(quadrature_rules);

      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  ReconstructedLocalizableFunctionType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               reconstructed_function,
                               range,
                               intersection_quadrature_rule,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply(true);
    } else {
      typedef typename SourceType::SpaceType::GridViewType GridViewType;
      const GridViewType& grid_view = source.space().grid_view();
      const auto intersection_quadrature_rule =
          Dune::QuadratureRules<RangeFieldType, GridViewType::dimension - 1>::rule(
              grid_view.ibegin(*(grid_view.template begin<0>()))->geometry().type(), 0);
      AdvectionLocalizableDefault<AnalyticalFluxType,
                                  NumericalCouplingFluxType,
                                  NumericalBoundaryFluxType,
                                  typename BoundaryValueFunctionType::NonparametricType,
                                  SourceType,
                                  RangeType>
          localizable_operator(analytical_flux,
                               current_boundary_values,
                               source,
                               range,
                               intersection_quadrature_rule,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...);
      localizable_operator.apply(true);
    }
  }
}; // struct AdvectionWENOOperatorApplier

template <class AnalyticalFluxType,
          class BoundaryValueFunctionType,
          class SourceType,
          class RangeType,
          class MatrixType,
          class Tuple,
          class FunctionType,
          size_t... Is>
constexpr auto apply_impl(const AnalyticalFluxType& analytical_flux,
                          const BoundaryValueFunctionType& boundary_values,
                          const SourceType& source,
                          RangeType& range,
                          const XT::Common::Parameter& param,
                          const bool use_linear_reconstruction,
                          std::shared_ptr<MatrixType> eigenvectors,
                          std::shared_ptr<MatrixType> eigenvectors_inverse,
                          Tuple t,
                          FunctionType f,
                          std::index_sequence<Is...>)
{
  return f(analytical_flux,
           boundary_values,
           source,
           range,
           param,
           use_linear_reconstruction,
           eigenvectors,
           eigenvectors_inverse,
           get<Is>(t)...);
}

template <class AnalyticalFluxType,
          class BoundaryValueFunctionType,
          class SourceType,
          class RangeType,
          class MatrixType,
          class Tuple,
          class FunctionType>
constexpr auto apply(const AnalyticalFluxType& analytical_flux,
                     const BoundaryValueFunctionType& boundary_values,
                     const SourceType& source,
                     RangeType& range,
                     const XT::Common::Parameter& param,
                     const bool use_linear_reconstruction,
                     std::shared_ptr<MatrixType> eigenvectors,
                     std::shared_ptr<MatrixType> eigenvectors_inverse,
                     Tuple t,
                     FunctionType f)
{
  return apply_impl(analytical_flux,
                    boundary_values,
                    source,
                    range,
                    param,
                    use_linear_reconstruction,
                    eigenvectors,
                    eigenvectors_inverse,
                    t,
                    f,
                    std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}


} // namespace internal

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          SlopeLimiters slope_lim>
class AdvectionLaxFriedrichsOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                                                         BoundaryValueFunctionImp,
                                                                                         LocalizableFunctionImp,
                                                                                         slope_lim>>
{
public:
  typedef internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                         BoundaryValueFunctionImp,
                                                         LocalizableFunctionImp,
                                                         slope_lim>
      Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionLaxFriedrichsOperator(const AnalyticalFluxType& analytical_flux,
                                 const BoundaryValueFunctionType& boundary_values,
                                 const LocalizableFunctionType& dx,
                                 const bool flux_is_linear = false,
                                 const bool use_linear_reconstruction = false,
                                 const bool use_local_laxfriedrichs_flux = false,
                                 const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , flux_is_linear_(flux_is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
    , use_local_laxfriedrichs_flux_(use_local_laxfriedrichs_flux)
    , entity_geometries_equal_(entity_geometries_equal)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, flux_is_linear, use_linear_reconstruction, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             use_linear_reconstruction_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             dx_,
                                                             param,
                                                             flux_is_linear_,
                                                             use_local_laxfriedrichs_flux_,
                                                             entity_geometries_equal_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const bool flux_is_linear_;
  const bool use_linear_reconstruction_;
  const bool use_local_laxfriedrichs_flux_;
  const bool entity_geometries_equal_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
}; // class AdvectionLaxFriedrichsOperator


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder,
          SlopeLimiters slope_lim>
class AdvectionLaxFriedrichsWENOOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionLaxFriedrichsWENOOperatorTraits<AnalyticalFluxImp,
                                                                                             BoundaryValueFunctionImp,
                                                                                             LocalizableFunctionImp,
                                                                                             GridViewType,
                                                                                             Polyhedron_3,
                                                                                             polOrder,
                                                                                             slope_lim>>
{
public:
  typedef internal::AdvectionLaxFriedrichsWENOOperatorTraits<AnalyticalFluxImp,
                                                             BoundaryValueFunctionImp,
                                                             LocalizableFunctionImp,
                                                             GridViewType,
                                                             Polyhedron_3,
                                                             polOrder,
                                                             slope_lim>
      Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionLaxFriedrichsWENOOperator(const AnalyticalFluxType& analytical_flux,
                                     const BoundaryValueFunctionType& boundary_values,
                                     const LocalizableFunctionType& dx,
                                     const GridViewType& grid_view,
                                     const FieldVector<size_t, dimDomain> grid_sizes,
                                     const Polyhedron_3& poly,
                                     const bool flux_is_linear = false,
                                     const bool use_reconstruction = false,
                                     const size_t quadrature_order = 2,
                                     const bool use_local_laxfriedrichs_flux = false,
                                     const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , grid_sizes_(grid_sizes)
    , poly_(poly)
    , flux_is_linear_(flux_is_linear)
    , use_reconstruction_(use_reconstruction)
    , quadrature_order_(quadrature_order)
    , use_local_laxfriedrichs_flux_(use_local_laxfriedrichs_flux)
    , entity_geometries_equal_(entity_geometries_equal)
    , entity_indices_(grid_view.size(0))
  {
    FieldVector<size_t, dimDomain> indices;
    for (const auto& entity : Dune::elements(grid_view)) {
      const auto& index = grid_view.indexSet().index(entity);
      const auto indices_array = entity.seed().impl().coord();
      for (size_t dd = 0; dd < dimDomain; ++dd)
        indices[dd] = indices_array[dd];
      entity_indices_[index] = indices;
    }
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
  {
    internal::AdvectionWENOOperatorApplier<NumericalCouplingFluxType,
                                           NumericalBoundaryFluxType,
                                           RangeFieldType,
                                           Polyhedron_3,
                                           dimDomain,
                                           dimRange,
                                           dimRangeCols,
                                           polOrder,
                                           slope_limiter>::apply(analytical_flux_,
                                                                 boundary_values_,
                                                                 source,
                                                                 range,
                                                                 param,
                                                                 use_reconstruction_,
                                                                 quadrature_order_,
                                                                 entity_indices_,
                                                                 grid_sizes_,
                                                                 poly_,
                                                                 dx_,
                                                                 param,
                                                                 flux_is_linear_,
                                                                 use_local_laxfriedrichs_flux_,
                                                                 entity_geometries_equal_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const FieldVector<size_t, dimDomain> grid_sizes_;
  const Polyhedron_3& poly_;
  const bool flux_is_linear_;
  const bool use_reconstruction_;
  const size_t quadrature_order_;
  const bool use_local_laxfriedrichs_flux_;
  const bool entity_geometries_equal_;
  std::vector<FieldVector<size_t, dimDomain>> entity_indices_;
}; // class AdvectionLaxFriedrichsWENOOperator


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class GridViewType,
          class Polyhedron_3,
          size_t polOrder,
          SlopeLimiters slope_lim>
class AdvectionGodunovWENOOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovWENOOperatorTraits<AnalyticalFluxImp,
                                                                                       BoundaryValueFunctionImp,
                                                                                       GridViewType,
                                                                                       Polyhedron_3,
                                                                                       polOrder,
                                                                                       slope_lim>>
{
public:
  typedef internal::AdvectionGodunovWENOOperatorTraits<AnalyticalFluxImp,
                                                       BoundaryValueFunctionImp,
                                                       GridViewType,
                                                       Polyhedron_3,
                                                       polOrder,
                                                       slope_lim>
      Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionGodunovWENOOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueFunctionType& boundary_values,
                               const GridViewType& grid_view,
                               const FieldVector<size_t, dimDomain> grid_sizes,
                               const Polyhedron_3& poly,
                               const bool flux_is_linear = false,
                               const bool use_reconstruction = false,
                               const size_t quadrature_order = 2)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_sizes_(grid_sizes)
    , poly_(poly)
    , flux_is_linear_(flux_is_linear)
    , use_reconstruction_(use_reconstruction)
    , quadrature_order_(quadrature_order)
    , entity_indices_(grid_view.size(0))
  {
    FieldVector<size_t, dimDomain> indices;
    for (const auto& entity : Dune::elements(grid_view)) {
      const auto& index = grid_view.indexSet().index(entity);
      const auto indices_array = entity.seed().impl().coord();
      for (size_t dd = 0; dd < dimDomain; ++dd)
        indices[dd] = indices_array[dd];
      entity_indices_[index] = indices;
    }
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
  {
    internal::AdvectionWENOOperatorApplier<NumericalCouplingFluxType,
                                           NumericalBoundaryFluxType,
                                           RangeFieldType,
                                           Polyhedron_3,
                                           dimDomain,
                                           dimRange,
                                           dimRangeCols,
                                           polOrder,
                                           slope_limiter>::apply(analytical_flux_,
                                                                 boundary_values_,
                                                                 source,
                                                                 range,
                                                                 param,
                                                                 use_reconstruction_,
                                                                 quadrature_order_,
                                                                 entity_indices_,
                                                                 grid_sizes_,
                                                                 poly_,
                                                                 param,
                                                                 flux_is_linear_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const FieldVector<size_t, dimDomain> grid_sizes_;
  const Polyhedron_3& poly_;
  const bool flux_is_linear_;
  const bool use_reconstruction_;
  const size_t quadrature_order_;
  std::vector<FieldVector<size_t, dimDomain>> entity_indices_;
}; // class AdvectionGodunovWENOOperator

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
class AdvectionKineticOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp,
                                                                                   BoundaryValueFunctionImp,
                                                                                   slope_lim>>
{
  typedef typename internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim>
      Traits;

public:
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux, const BoundaryValueFunctionType& boundary_values)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             false,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             param);
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
}; // class AdvectionKineticOperator


// TODO: 0 boundary by default, so no need to specify boundary conditions for periodic grid views
template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
class AdvectionGodunovOperator
    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                                   BoundaryValueFunctionImp,
                                                                                   slope_lim>>
{
public:
  typedef internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim> Traits;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

protected:
  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

public:
  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueFunctionType& boundary_values,
                           const bool flux_is_linear = false,
                           const bool use_linear_reconstruction = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , flux_is_linear_(flux_is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
  {
    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType, AnalyticalFluxType>::initialize(
        analytical_flux_, flux_is_linear, use_linear_reconstruction, eigenvectors_, eigenvectors_inverse_);
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param = {}) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       RangeFieldType,
                                       dimRange,
                                       dimRangeCols,
                                       slope_limiter>::apply(analytical_flux_,
                                                             boundary_values_,
                                                             source,
                                                             range,
                                                             param,
                                                             use_linear_reconstruction_,
                                                             eigenvectors_,
                                                             eigenvectors_inverse_,
                                                             flux_is_linear_);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const bool flux_is_linear_;
  const bool use_linear_reconstruction_;
  std::shared_ptr<MatrixType> eigenvectors_;
  std::shared_ptr<MatrixType> eigenvectors_inverse_;
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
class AdvectionKineticOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter>
class AdvectionGodunovOperator
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN

template <class GridViewType, class MatrixType, class DiscreteFunctionType>
class MatrixSolveFunctor : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename XT::Grid::Functor::Codim0<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixSolveFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& rhs,
                     DiscreteFunctionType& solution)
    : matrices_(matrices)
    , rhs_(rhs)
    , solution_(solution)
  {
  }

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = rhs_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_solution(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_rhs(local_solution.size(), 0.);
    const auto& rhs_vector = rhs_.vector();
    auto& solution_vector = solution_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      local_rhs[ii] = rhs_vector.get_entry(global_indices[ii]);
    // solve
    matrices_[rhs_.space().grid_view().indexSet().index(entity)].solve(local_solution, local_rhs);
    // write solution
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      solution_vector.set_entry(global_indices[ii], local_solution[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& rhs_;
  DiscreteFunctionType& solution_;
};

template <class DiscreteFunctionType, class RhsEvaluationType>
class MatrixExponentialFunctor
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridViewType>
{
  typedef typename XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridViewType> BaseType;
  typedef typename RhsEvaluationType::RangeFieldType FieldType;
  static const size_t dimRange = RhsEvaluationType::dimRange;

public:
  using typename BaseType::EntityType;

  MatrixExponentialFunctor(DiscreteFunctionType& solution,
                           const double t,
                           const double dt,
                           const RhsEvaluationType& rhs_evaluation)
    : solution_(solution)
    , t_(t)
    , dt_(dt)
    , rhs_evaluation_(rhs_evaluation)
  {
  }

  // solves d_t u(t) = A u(t) + b locally on each entity
  virtual void apply_local(const EntityType& entity)
  {
    auto solution_local = solution_.local_discrete_function(entity);

    // get u
    const auto center = entity.geometry().local(entity.geometry().center());
    const auto u0 = solution_local->evaluate(center);

    // get A and b
    auto zero = u0;
    zero *= 0.;
    const FieldMatrix<FieldType, dimRange, dimRange> A = rhs_evaluation_.jacobian(u0, entity, center);
    const auto b = rhs_evaluation_.evaluate(zero, entity, center);

    // calculate matrix exponential exp(A*dt)
    auto Adt = A;
    Adt *= dt_;
    // get pointer to the underlying array of the FieldMatrix
    double* Adt_array = &(Adt[0][0]);

    const double* exp_Adt_array = r8mat_expm1(dimRange, Adt_array);

    FieldMatrix<FieldType, dimRange, dimRange> exp_Adt;
    std::copy_n(exp_Adt_array, dimRange * dimRange, &(exp_Adt[0][0]));
    delete[] exp_Adt_array;

    // get inverse of A
    auto A_inverse = A;
    A_inverse.invert();

    // calculate solution u = exp(A dt) (u0 + A^{-1} (exp(A dt) - I) b)
    FieldVector<FieldType, dimRange> u(0), ret;
    auto exp_Adt_m_I = exp_Adt;
    for (size_t ii = 0; ii < dimRange; ++ii)
      exp_Adt_m_I[ii][ii] -= 1.;
    exp_Adt_m_I.leftmultiply(A_inverse);
    exp_Adt_m_I.mv(b, u);
    u += u0;
    exp_Adt.mv(u, ret);

    // write to return vector
    auto& local_vector = solution_local->vector();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_vector.set(ii, ret[ii]);
  }

private:
  DiscreteFunctionType& solution_;
  const double t_;
  const double dt_;
  const RhsEvaluationType& rhs_evaluation_;
};


template <class GridViewType, class MatrixType, class DiscreteFunctionType>
class MatrixApplyFunctor : public XT::Grid::Functor::Codim0<GridViewType>
{
  typedef typename XT::Grid::Functor::Codim0<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixApplyFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& vector,
                     DiscreteFunctionType& result)
    : matrices_(matrices)
    , vector_(vector)
    , result_(result)
  {
  }

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = vector_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_vector(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_result(local_vector.size(), 0.);
    const auto& vector_vector = vector_.vector();
    auto& result_vector = result_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      local_vector[ii] = vector_vector.get_entry(global_indices[ii]);
    // solve
    matrices_[vector_.space().grid_view().indexSet().index(entity)].mv(local_vector, local_result);

    // write solution
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      result_vector.set_entry(global_indices[ii], local_result[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& vector_;
  DiscreteFunctionType& result_;
};

template <class RHSEvaluationImp>
class AdvectionRHSOperator : public Dune::GDT::OperatorInterface<internal::AdvectionRHSOperatorTraits<RHSEvaluationImp>>
{
  static_assert(is_rhs_evaluation<RHSEvaluationImp>::value, "RHSEvaluationImp has to be derived from RHSInterface!");

public:
  typedef internal::AdvectionRHSOperatorTraits<RHSEvaluationImp> Traits;
  typedef typename Traits::RHSEvaluationType RHSEvaluationType;
  using BaseType = typename Dune::GDT::OperatorInterface<internal::AdvectionRHSOperatorTraits<RHSEvaluationImp>>;

  AdvectionRHSOperator(const RHSEvaluationType& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& /*param*/) const
  {
    std::fill(range.vector().begin(), range.vector().end(), 0);
    LocalVolumeIntegralFunctional<LocalFvRhsIntegrand<RHSEvaluationType, SourceType>> local_functional(rhs_evaluation_,
                                                                                                       source);
    VectorFunctionalBase<typename RangeType::VectorType,
                         typename RangeType::SpaceType,
                         typename RangeType::SpaceType::GridViewType,
                         typename RangeType::DomainFieldType>
        functional_assembler(range.vector(), range.space());
    functional_assembler.append(local_functional);
    functional_assembler.assemble(true);
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixTraits>
  void assemble_jacobian(XT::LA::MatrixInterface<MatrixTraits, typename SourceType::RangeFieldType>& jac,
                         const SourceType& source,
                         const XT::Common::Parameter& /*param*/ = {}) const
  {
    typedef LocalVolumeIntegralOperator<LocalFvRhsJacobianIntegrand<RHSEvaluationType, SourceType>> LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source);
    LocalVolumeTwoFormAssembler<LocalOperatorType> local_assembler(local_operator);
    SystemAssembler<typename SourceType::SpaceType> assembler(source.space());
    assembler.append(local_assembler, jac);
    assembler.assemble(true);
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixType>
  void assemble_newton_matrix(std::vector<MatrixType>& newton_matrices,
                              const SourceType& source,
                              const XT::Common::Parameter& param) const
  {
    typedef LocalVolumeIntegralOperator<LocalFvRhsNewtonIntegrand<RHSEvaluationType, SourceType>> LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source, param);
    LocalVolumeTwoFormAssembler<LocalOperatorType> local_assembler(local_operator);
    SystemAssembler<typename SourceType::SpaceType> assembler(source.space());
    assembler.append(local_assembler, newton_matrices);
    assembler.assemble(false);
  }

  // solves with matrix exponential on each entity
  template <class SourceType>
  void apply_matrix_exponential(SourceType& solution,
                                const double t,
                                const double dt,
                                const XT::Common::Parameter& /*param*/ = {}) const
  {
    MatrixExponentialFunctor<SourceType, RHSEvaluationType> functor(solution, t, dt, rhs_evaluation_);
    SystemAssembler<typename SourceType::SpaceType> assembler(solution.space());
    assembler.append(functor);
    assembler.assemble(false);
  }

  // solves with local jacobian on each entity
  template <class SourceType, class RangeType, class MatrixType>
  void solve(const std::vector<MatrixType>& newton_matrices,
             const SourceType& rhs,
             RangeType& solution,
             const XT::Common::Parameter& /*param*/ = {}) const
  {
    MatrixSolveFunctor<typename SourceType::SpaceType::GridViewType, MatrixType, SourceType> functor(
        newton_matrices, rhs, solution);
    SystemAssembler<typename SourceType::SpaceType> assembler(rhs.space());
    assembler.append(functor);
    assembler.assemble(false);
  }

  // applies local jacobian on each entity
  template <class SourceType, class RangeType, class MatrixType>
  void mv(const std::vector<MatrixType>& newton_matrices,
          const SourceType& vector,
          RangeType& result,
          const XT::Common::Parameter& /*param*/ = {}) const
  {
    MatrixApplyFunctor<typename SourceType::SpaceType::GridViewType, MatrixType, SourceType> functor(
        newton_matrices, vector, result);
    SystemAssembler<typename SourceType::SpaceType> assembler(vector.space());
    assembler.append(functor);
    assembler.assemble(false);
  }

private:
  const RHSEvaluationType& rhs_evaluation_;
}; // class AdvectionRHSOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_HH
