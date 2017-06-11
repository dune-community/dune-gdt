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

#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/basisfunctions.hh>

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
  local_laxfriedrichs_with_reconstruction,
  kinetic
};

template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class Traits>
class AdvectionLaxFriedrichsOperator;

// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder = 1,
//          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
// class AdvectionGodunovWENOOperator;


// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter =
// SlopeLimiters::minmod>
// class AdvectionGodunovOperator;

// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder = 1,
//          SlopeLimiters slope_limiter = SlopeLimiters::minmod>
// class AdvectionKineticWENOOperator;

// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter =
// SlopeLimiters::minmod>
// class AdvectionKineticOperator;

template <class RhsEvaluationImp>
class AdvectionRhsOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp>
class AdvectionTraitsBase
{
public:
  static const size_t polOrder = reconstructionOrder;
  static const SlopeLimiters slope_limiter = slope_lim;
  static const bool realizability_limiting = realizability_lim;
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueType;
  typedef BasisFunctionImp BasisFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename BoundaryValueFunctionImp::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueFunctionImp::RangeFieldType RangeFieldType;
  typedef RangeFieldType FieldType;
  typedef typename BoundaryValueFunctionImp::DomainType DomainType;
  typedef typename AnalyticalFluxType::JacobianWrtURangeType JacobianType;
}; // class AdvectionTraitsBase


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp>
class AdvectionLaxFriedrichsOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                        BoundaryValueFunctionImp,
                                                                        reconstructionOrder,
                                                                        slope_lim,
                                                                        realizability_lim,
                                                                        BasisFunctionImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueFunctionImp,
                              reconstructionOrder,
                              slope_lim,
                              realizability_lim,
                              BasisFunctionImp>
      BaseType;

public:
  using BaseType::polOrder;
  using BaseType::slope_limiter;
  using BaseType::realizability_limiting;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::BasisFunctionType;
  using BaseType::dimDomain;
  typedef typename Dune::GDT::LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxType,
                                                                      LocalizableFunctionType,
                                                                      dimDomain>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                                               BoundaryValueType,
                                                                               LocalizableFunctionType,
                                                                               dimDomain>
      NumericalBoundaryFluxType;
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxImp,
                                         BoundaryValueFunctionImp,
                                         LocalizableFunctionImp,
                                         polOrder,
                                         slope_limiter,
                                         realizability_limiting,
                                         BasisFunctionType,
                                         AdvectionLaxFriedrichsOperatorTraits>
      derived_type;
}; // class AdvectionLaxFriedrichsOperatorTraits

// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
// class AdvectionGodunovOperatorTraits
//    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
//{
//  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

// public:
//  using BaseType::slope_limiter;
//  using typename BaseType::AnalyticalFluxType;
//  using typename BaseType::BoundaryValueType;
//  using BaseType::dimDomain;
//  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
//      NumericalCouplingFluxType;
//  typedef
//      typename Dune::GDT::GodunovLocalNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType, dimDomain>
//          NumericalBoundaryFluxType;
//  typedef AdvectionGodunovOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
//}; // class AdvectionGodunovOperatorTraits

// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder,
//          SlopeLimiters slope_limiter_type>
// class AdvectionGodunovWENOOperatorTraits
//    : public AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
//{
// public:
//  typedef AdvectionGodunovWENOOperator<AnalyticalFluxImp,
//                                       BoundaryValueFunctionImp,
//                                       GridLayerType,
//                                       basis_function_type,
//                                       polOrder,
//                                       slope_limiter_type>
//      derived_type;
//}; // class AdvectionGodunovWENOOperatorTraits


// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter_type>
// class AdvectionKineticOperatorTraits
//    : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
//{
//  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type> BaseType;

// public:
//  using BaseType::slope_limiter;
//  using typename BaseType::AnalyticalFluxType;
//  using typename BaseType::BoundaryValueType;
//  using BaseType::dimDomain;
//  typedef typename Dune::GDT::KineticLocalNumericalCouplingFlux<AnalyticalFluxType, dimDomain>
//      NumericalCouplingFluxType;
//  typedef
//      typename Dune::GDT::KineticLocalNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType, dimDomain>
//          NumericalBoundaryFluxType;
//  typedef AdvectionKineticOperator<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter> derived_type;
//}; // class AdvectionKineticOperatorTraits

// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder,
//          SlopeLimiters slope_limiter_type>
// class AdvectionKineticWENOOperatorTraits
//    : public AdvectionKineticOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter_type>
//{
// public:
//  typedef AdvectionKineticWENOOperator<AnalyticalFluxImp,
//                                       BoundaryValueFunctionImp,
//                                       GridLayerType,
//                                       basis_function_type,
//                                       polOrder,
//                                       slope_limiter_type>
//      derived_type;
//}; // class AdvectionKineticWENOOperatorTraits


template <class RhsEvaluationImp>
class AdvectionRhsOperatorTraits
{
public:
  typedef AdvectionRhsOperator<RhsEvaluationImp> derived_type;
  typedef RhsEvaluationImp RhsEvaluationType;
  typedef typename RhsEvaluationImp::DomainFieldType FieldType;
  typedef typename RhsEvaluationImp::JacobianWrtURangeType JacobianType;
}; // class AdvectionRhsOperatorTraits


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

  //  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
  //                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
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
  typedef BoundaryValueFunctionImp BoundaryValueType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridLayerType GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename Dune::GDT::LocalCouplingFvOperator<NumericalCouplingFluxType> LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFvOperator<NumericalBoundaryFluxType> LocalBoundaryOperatorType;

  template <class QuadratureRuleType, class... LocalOperatorArgTypes>
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const BoundaryValueType& boundary_values,
                              const SourceType& source,
                              RangeType& range,
                              const XT::Common::Parameter& param,
                              const QuadratureRuleType& quadrature_rule,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(range.space().grid_layer(), source, range)
    , local_operator_(
          quadrature_rule, analytical_flux, param, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
    , local_boundary_operator_(quadrature_rule,
                               analytical_flux,
                               boundary_values,
                               param,
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


// TODO: remove eigen dependency of GodunovLocalNumericalCouplingFlux/GodunovLocalNumericalBoundaryFlux
#if HAVE_EIGEN


namespace internal {


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

  typedef GridViewImp GridLayerType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate LocalCoordinateType;

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

    virtual void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& /*param*/) const
    {
      try {
        ret = values_.at(xx);
      } catch (const std::out_of_range& e) {
        DUNE_THROW(Dune::RangeError,
                   "There are no values for xx = " << XT::Common::to_string(xx) << " in this function!");
      }
    }

    virtual void
    jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/, const XT::Common::Parameter& /*param*/) const
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
      const GridLayerType& grid_layer,
      const std::vector<std::map<LocalCoordinateType, RangeType, internal::FieldVectorLess>>& reconstructed_values)
    : index_set_(grid_layer.indexSet())
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


template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          bool realizability_limiting>
struct AdvectionOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class SourceType,
            class RangeType,
            class RangeFieldType,
            class DomainFieldType,
            class BasisFunctionType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        const Dune::QuadratureRule<DomainFieldType, 1> intersection_quadrature_1d,
        const Dune::QuadratureRule<DomainFieldType, BoundaryValueType::dimDomain - 1> intersection_quadrature,
        const Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain>& quadrature,
        const RangeFieldType epsilon,
        const std::shared_ptr<const BasisFunctionType> basis_functions,
        LocalOperatorArgTypes&&... local_operator_args)

  {
    typedef XT::LA::EigenDenseVector<RangeFieldType> EigenVectorType;
    typedef typename SourceType::SpaceType::GridLayerType GridLayerType;
    typedef typename BoundaryValueType::DomainType DomainType;
    static const size_t dimDomain = BoundaryValueType::dimDomain;
    static const size_t dimRange = BoundaryValueType::dimRange;
    const GridLayerType& grid_layer = source.space().grid_layer();

    // evaluate cell averages as EigenVectorType
    std::vector<EigenVectorType> source_values(grid_layer.indexSet().size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto& entity_index = grid_layer.indexSet().index(entity);
      const auto& local_source = source.local_function(entity);
      source_values[entity_index] = XT::LA::internal::FieldVectorToLaVector<EigenVectorType, dimRange>::convert(
          local_source->evaluate(entity.geometry().local(entity.geometry().center())));
    }

    // do reconstruction
    std::vector<std::map<DomainType, typename BoundaryValueType::RangeType, internal::FieldVectorLess>>
        reconstructed_values(grid_layer.size(0));

    auto local_reconstruction_operator =
        LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>(
            source_values,
            analytical_flux,
            boundary_values,
            grid_layer,
            param,
            intersection_quadrature_1d,
            reconstructed_values);
    auto walker = XT::Grid::Walker<GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);

    if (realizability_limiting) {
      assert(basis_functions);
      // do limiting for realizability in M_N models
      auto local_realizability_limiter = LocalRealizabilityLimiter<SourceType, BasisFunctionType, dimDomain, dimRange>(
          source, reconstructed_values, *basis_functions, quadrature, epsilon);
      walker.clear();
      walker.append(local_realizability_limiter);
      walker.walk(true);
    }

    typedef ReconstructedLocalizableFunction<GridLayerType, RangeFieldType, dimDomain, RangeFieldType, dimRange>
        ReconstructedLocalizableFunctionType;
    const ReconstructedLocalizableFunctionType reconstructed_function(grid_layer, reconstructed_values);

    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                ReconstructedLocalizableFunctionType,
                                RangeType>
        localizable_operator(analytical_flux,
                             boundary_values,
                             reconstructed_function,
                             range,
                             param,
                             intersection_quadrature,
                             std::forward<LocalOperatorArgTypes>(local_operator_args)...);
    localizable_operator.apply(true);
  }
}; // struct AdvectionOperatorApplier<..., polOrder>0,...>

template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          SlopeLimiters slope_limiter,
          bool realizability_limiting>
struct AdvectionOperatorApplier<NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                0,
                                slope_limiter,
                                realizability_limiting>
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class SourceType,
            class RangeType,
            class RangeFieldType,
            class DomainFieldType,
            class BasisFunctionType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        const Dune::QuadratureRule<DomainFieldType, 1> /*intersection_quadrature_1d*/,
        const Dune::QuadratureRule<DomainFieldType, BoundaryValueType::dimDomain - 1> intersection_quadrature,
        const Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain> /*quadrature*/,
        const RangeFieldType /*epsilon*/,
        const std::shared_ptr<const BasisFunctionType> /*basis_functions*/,
        LocalOperatorArgTypes&&... local_operator_args)
  {
    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                SourceType,
                                RangeType>
        localizable_operator(analytical_flux,
                             boundary_values,
                             source,
                             range,
                             param,
                             intersection_quadrature,
                             std::forward<LocalOperatorArgTypes>(local_operator_args)...);
    localizable_operator.apply(true);
  }
}; // struct AdvectionOperatorApplier<..., polOrder=0,...>


} // namespace internal


template <class Traits>
class AdvectionOperatorBase
{
public:
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const size_t polOrder = Traits::polOrder;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  static const bool realizability_limiting = Traits::realizability_limiting;
  typedef typename Traits::BasisFunctionType BasisFunctionType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

  typedef Dune::QuadratureRule<DomainFieldType, 1> Intersection1dQuadratureType;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain - 1> IntersectionQuadratureType;
  typedef Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain> QuadratureType;

public:
  AdvectionOperatorBase(const AnalyticalFluxType& analytical_flux, const BoundaryValueType& boundary_values)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , intersection_1d_quadrature_(helper<>::default_1d_quadrature())
    , intersection_quadrature_(helper2<>::get_quadrature(intersection_1d_quadrature_))
    , epsilon_(1e-14)
  {
  }

  template <class SourceType, class RangeType, class... Args>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param, Args&&... args) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       polOrder,
                                       slope_limiter,
                                       realizability_limiting>::apply(analytical_flux_,
                                                                      boundary_values_,
                                                                      source,
                                                                      range,
                                                                      param,
                                                                      intersection_1d_quadrature_,
                                                                      intersection_quadrature_,
                                                                      quadrature_,
                                                                      epsilon_,
                                                                      basis_functions_,
                                                                      std::forward<Args>(args)...);
  }

  void set_1d_quadrature(const Intersection1dQuadratureType& quadrature)
  {
    intersection_1d_quadrature_ = quadrature;
    intersection_quadrature_ = helper2<>::get_quadrature(quadrature);
  }

  void set_quadrature(const QuadratureType& quadrature)
  {
    quadrature_ = quadrature;
  }

  void set_epsilon(const RangeFieldType& epsilon)
  {
    epsilon_ = epsilon;
  }

  void set_basisfunctions(const std::shared_ptr<const BasisFunctionType> basis_functions)
  {
    basis_functions_ = basis_functions;
  }

private:
  template <size_t reconstructionOrder = polOrder, class anything = void>
  struct helper
  {
    static Intersection1dQuadratureType default_1d_quadrature()
    {
      return Dune::QuadratureRules<DomainFieldType, 1>::rule(Dune::GeometryType(Dune::GeometryType::BasicType::cube, 1),
                                                             2 * polOrder);
    }
  };

  template <class anything>
  struct helper<1, anything>
  {
    static Intersection1dQuadratureType default_1d_quadrature()
    {
      // get 1D quadrature rules
      Intersection1dQuadratureType quadrature;
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. - 1. / std::sqrt(3)), 0.5));
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. + 1. / std::sqrt(3)), 0.5));
      return quadrature;
    }
  };

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper2;

  template <class anything>
  struct helper2<1, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& /*quadrature_1d*/)
    {
      Dune::QuadratureRule<DomainFieldType, dimDomain - 1> ret;
      ret.push_back(Dune::QuadraturePoint<DomainFieldType, 0>(FieldVector<DomainFieldType, 0>(0), 1));
      return ret;
    }
  };

  template <class anything>
  struct helper2<2, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& quadrature_1d)
    {
      return quadrature_1d;
    }
  };

  template <class anything>
  struct helper2<3, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& quadrature_1d)
    {
      Dune::QuadratureRule<DomainFieldType, dimDomain - 1> ret;
      for (size_t ii = 0; ii < quadrature_1d.size(); ++ii)
        for (size_t jj = 0; jj < quadrature_1d.size(); ++jj)
          ret.push_back(Dune::QuadraturePoint<DomainFieldType, dimDomain - 1>(
              {quadrature_1d[ii].position()[0], quadrature_1d[jj].position()[0]},
              quadrature_1d[ii].weight() * quadrature_1d[jj].weight()));
      return ret;
    }
  };

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  Intersection1dQuadratureType intersection_1d_quadrature_;
  IntersectionQuadratureType intersection_quadrature_;
  QuadratureType quadrature_;
  std::shared_ptr<const BasisFunctionType> basis_functions_;
  RangeFieldType epsilon_;
}; // class AdvectionOperatorBase<...>


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          bool realizability_lim = false,
          class BasisFunctionImp = Hyperbolic::Problems::HatFunctions<typename BoundaryValueFunctionImp::DomainFieldImp,
                                                                      BoundaryValueFunctionImp::dimDomain,
                                                                      typename BoundaryValueFunctionImp::RangeFieldType,
                                                                      BoundaryValueFunctionImp::dimRange,
                                                                      BoundaryValueFunctionImp::dimRangeCols>,
          class Traits = internal::AdvectionLaxFriedrichsOperatorTraits<AnalyticalFluxImp,
                                                                        BoundaryValueFunctionImp,
                                                                        LocalizableFunctionImp,
                                                                        polOrder,
                                                                        slope_lim,
                                                                        realizability_lim,
                                                                        BasisFunctionImp>>
class AdvectionLaxFriedrichsOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionLaxFriedrichsOperator(const AnalyticalFluxType& analytical_flux,
                                 const BoundaryValueType& boundary_values,
                                 const LocalizableFunctionType& dx,
                                 const bool use_local_laxfriedrichs_flux = false,
                                 const bool is_linear = false,
                                 const DomainType lambda = DomainType(0))
    : BaseType(analytical_flux, boundary_values)
    , dx_(dx)
    , use_local_laxfriedrichs_flux_(use_local_laxfriedrichs_flux)
    , is_linear_(is_linear)
    , lambda_(lambda)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, use_local_laxfriedrichs_flux_, is_linear_, lambda_);
  }

private:
  const LocalizableFunctionType& dx_;
  const bool use_local_laxfriedrichs_flux_;
  const bool is_linear_;
  const DomainType lambda_;
}; // class AdvectionLaxFriedrichsOperator


// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder,
//          SlopeLimiters slope_lim>
// class AdvectionGodunovWENOOperator
//    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovWENOOperatorTraits<AnalyticalFluxImp,
//                                                                                       BoundaryValueFunctionImp,
//                                                                                       GridLayerType,
//                                                                                       basis_function_type,
//                                                                                       polOrder,
//                                                                                       slope_lim>>
//{
// public:
//  typedef internal::AdvectionGodunovWENOOperatorTraits<AnalyticalFluxImp,
//                                                       BoundaryValueFunctionImp,
//                                                       GridLayerType,
//                                                       basis_function_type,
//                                                       polOrder,
//                                                       slope_lim>
//      Traits;
//  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
//  typedef typename Traits::BoundaryValueType BoundaryValueType;
//  static const size_t dimDomain = Traits::dimDomain;
//  static const size_t dimRange = Traits::dimRange;
//  static const size_t dimRangeCols = Traits::dimRangeCols;
//  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
//  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
//  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
//  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

// protected:
//  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
//  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

// public:
//  AdvectionGodunovWENOOperator(
//      const AnalyticalFluxType& analytical_flux,
//      const BoundaryValueType& boundary_values,
//      const GridLayerType& grid_layer,
//      const FieldVector<size_t, dimDomain> grid_sizes,
//      const std::vector<std::pair<FieldVector<RangeFieldType, dimRange>, RangeFieldType>>& plane_coefficients,
//      const bool flux_is_linear = false,
//      const bool use_reconstruction = false,
//      const FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain> quadrature_rules =
//          FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain>())
//    : analytical_flux_(analytical_flux)
//    , boundary_values_(boundary_values)
//    , grid_sizes_(grid_sizes)
//    , plane_coefficients_(plane_coefficients)
//    , flux_is_linear_(flux_is_linear)
//    , use_reconstruction_(use_reconstruction)
//    , quadrature_rules_(quadrature_rules)
//    , entity_indices_(grid_layer.size(0))
//  {
//    FieldVector<size_t, dimDomain> indices;
//    for (const auto& entity : Dune::elements(grid_layer)) {
//      const auto& index = grid_layer.indexSet().index(entity);
//      const auto indices_array = entity.seed().impl().coord();
//      for (size_t dd = 0; dd < dimDomain; ++dd)
//        indices[dd] = indices_array[dd];
//      entity_indices_[index] = indices;
//    }
//    if (quadrature_rules_[0].empty())
//      quadrature_rules_ = default_quadrature_rules(grid_layer);
//  }

//  template <class SourceType, class RangeType>
//  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
//  {

//  }

// private:
//  const AnalyticalFluxType& analytical_flux_;
//  const BoundaryValueType& boundary_values_;
//  const FieldVector<size_t, dimDomain> grid_sizes_;
//  const std::vector<std::pair<FieldVector<RangeFieldType, dimRange>, RangeFieldType>>& plane_coefficients_;
//  const bool flux_is_linear_;
//  const bool use_reconstruction_;
//  FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain> quadrature_rules_;
//  std::vector<FieldVector<size_t, dimDomain>> entity_indices_;
//}; // class AdvectionGodunovWENOOperator

// template <class AnalyticalFluxImp,
//          class BoundaryValueFunctionImp,
//          class GridLayerType,
//          BasisFunction basis_function_type,
//          size_t polOrder,
//          SlopeLimiters slope_lim>
// class AdvectionKineticWENOOperator
//    : public Dune::GDT::OperatorInterface<internal::AdvectionKineticWENOOperatorTraits<AnalyticalFluxImp,
//                                                                                       BoundaryValueFunctionImp,
//                                                                                       GridLayerType,
//                                                                                       basis_function_type,
//                                                                                       polOrder,
//                                                                                       slope_lim>>
//{
// public:
//  typedef internal::AdvectionKineticWENOOperatorTraits<AnalyticalFluxImp,
//                                                       BoundaryValueFunctionImp,
//                                                       GridLayerType,
//                                                       basis_function_type,
//                                                       polOrder,
//                                                       slope_lim>
//      Traits;
//  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
//  typedef typename Traits::BoundaryValueType BoundaryValueType;
//  static const size_t dimDomain = Traits::dimDomain;
//  static const size_t dimRange = Traits::dimRange;
//  static const size_t dimRangeCols = Traits::dimRangeCols;
//  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
//  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
//  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
//  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

// protected:
//  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
//  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

// public:
//  AdvectionKineticWENOOperator(
//      const AnalyticalFluxType& analytical_flux,
//      const BoundaryValueType& boundary_values,
//      const GridLayerType& grid_layer,
//      const FieldVector<size_t, dimDomain> grid_sizes,
//      const std::vector<std::pair<FieldVector<RangeFieldType, dimRange>, RangeFieldType>>& plane_coefficients,
//      const bool flux_is_linear = false,
//      const bool use_reconstruction = false,
//      const FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain> quadrature_rules =
//          FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain>(),
//      const RangeFieldType epsilon = 1e-10)
//    : analytical_flux_(analytical_flux)
//    , boundary_values_(boundary_values)
//    , grid_sizes_(grid_sizes)
//    , plane_coefficients_(plane_coefficients)
//    , flux_is_linear_(flux_is_linear)
//    , use_reconstruction_(use_reconstruction)
//    , quadrature_rules_(quadrature_rules)
//    , epsilon_(epsilon)
//    , entity_indices_(grid_layer.size(0))
//  {
//    FieldVector<size_t, dimDomain> indices;
//    for (const auto& entity : Dune::elements(grid_layer)) {
//      const auto& index = grid_layer.indexSet().index(entity);
//      const auto indices_array = entity.seed().impl().coord();
//      for (size_t dd = 0; dd < dimDomain; ++dd)
//        indices[dd] = indices_array[dd];
//      entity_indices_[index] = indices;
//    }
//    if (quadrature_rules_[0].empty())
//      quadrature_rules_ = default_quadrature_rules(grid_layer);
//  }

//  template <class SourceType, class RangeType>
//  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
//  {
//    internal::AdvectionWENOOperatorApplier<NumericalCouplingFluxType,
//                                           NumericalBoundaryFluxType,
//                                           RangeFieldType,
//                                           basis_function_type,
//                                           dimDomain,
//                                           dimRange,
//                                           dimRangeCols,
//                                           polOrder,
//                                           slope_limiter>::apply(analytical_flux_,
//                                                                 boundary_values_,
//                                                                 source,
//                                                                 range,
//                                                                 param,
//                                                                 use_reconstruction_,
//                                                                 quadrature_rules_,
//                                                                 epsilon_,
//                                                                 entity_indices_,
//                                                                 grid_sizes_,
//                                                                 plane_coefficients_,
//                                                                 param);
//  }

// private:
//  const AnalyticalFluxType& analytical_flux_;
//  const BoundaryValueType& boundary_values_;
//  const FieldVector<size_t, dimDomain> grid_sizes_;
//  const std::vector<std::pair<FieldVector<RangeFieldType, dimRange>, RangeFieldType>>& plane_coefficients_;
//  const bool flux_is_linear_;
//  const bool use_reconstruction_;
//  FieldVector<Dune::QuadratureRule<RangeFieldType, 1>, dimDomain> quadrature_rules_;
//  const RangeFieldType epsilon_;
//  std::vector<FieldVector<size_t, dimDomain>> entity_indices_;
//}; // class AdvectionKineticWENOOperator


// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
// class AdvectionKineticOperator
//    : public Dune::GDT::OperatorInterface<internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp,
//                                                                                   BoundaryValueFunctionImp,
//                                                                                   slope_lim>>
//{
//  typedef typename internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim>
//      Traits;

// public:
//  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
//  typedef typename Traits::BoundaryValueType BoundaryValueType;
//  static const size_t dimDomain = Traits::dimDomain;
//  static const size_t dimRange = Traits::dimRange;
//  static const size_t dimRangeCols = Traits::dimRangeCols;
//  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
//  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
//  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
//  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

// protected:
//  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
//  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

// public:
//  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux, const BoundaryValueType&
//  boundary_values)
//    : analytical_flux_(analytical_flux)
//    , boundary_values_(boundary_values)
//  {
//  }

//  template <class SourceType, class RangeType>
//  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param) const
//  {
//    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
//                                       NumericalBoundaryFluxType,
//                                       RangeFieldType,
//                                       dimRange,
//                                       dimRangeCols,
//                                       slope_limiter>::apply(analytical_flux_,
//                                                             boundary_values_,
//                                                             source,
//                                                             range,
//                                                             param,
//                                                             false,
//                                                             eigenvectors_,
//                                                             eigenvectors_inverse_,
//                                                             param);
//  }

//  const AnalyticalFluxType& analytical_flux_;
//  const BoundaryValueType& boundary_values_;
//  std::shared_ptr<MatrixType> eigenvectors_;
//  std::shared_ptr<MatrixType> eigenvectors_inverse_;
//}; // class AdvectionKineticOperator


//// TODO: 0 boundary by default, so no need to specify boundary conditions for periodic grid views
// template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_lim>
// class AdvectionGodunovOperator
//    : public Dune::GDT::OperatorInterface<internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
//                                                                                   BoundaryValueFunctionImp,
//                                                                                   slope_lim>>
//{
// public:
//  typedef internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueFunctionImp, slope_lim> Traits;
//  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
//  typedef typename Traits::BoundaryValueType BoundaryValueType;
//  static const size_t dimDomain = Traits::dimDomain;
//  static const size_t dimRange = Traits::dimRange;
//  static const size_t dimRangeCols = Traits::dimRangeCols;
//  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
//  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
//  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
//  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;

// protected:
//  typedef typename Dune::XT::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
//  typedef typename Dune::XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

// public:
//  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
//                           const BoundaryValueType& boundary_values,
//                           const bool flux_is_linear = false,
//                           const bool use_linear_reconstruction = false)
//    : analytical_flux_(analytical_flux)
//    , boundary_values_(boundary_values)
//    , flux_is_linear_(flux_is_linear)
//    , use_linear_reconstruction_(use_linear_reconstruction)
//  {
//    //    internal::EigenvectorInitializer<dimDomain, dimRange, MatrixType, EigenMatrixType,
//    //    AnalyticalFluxType>::initialize(
//    //        analytical_flux_, flux_is_linear, use_linear_reconstruction, eigenvectors_, eigenvectors_inverse_);
//  }

//  template <class SourceType, class RangeType>
//  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter param = {}) const
//  {
//    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
//                                       NumericalBoundaryFluxType,
//                                       RangeFieldType,
//                                       dimRange,
//                                       dimRangeCols,
//                                       slope_limiter>::apply(analytical_flux_,
//                                                             boundary_values_,
//                                                             source,
//                                                             range,
//                                                             param,
//                                                             use_linear_reconstruction_,
//                                                             eigenvectors_,
//                                                             eigenvectors_inverse_,
//                                                             param,
//                                                             flux_is_linear_);
//  }

// private:
//  const AnalyticalFluxType& analytical_flux_;
//  const BoundaryValueType& boundary_values_;
//  const bool flux_is_linear_;
//  const bool use_linear_reconstruction_;
//  std::shared_ptr<MatrixType> eigenvectors_;
//  std::shared_ptr<MatrixType> eigenvectors_inverse_;
//}; // class AdvectionGodunovOperator

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

template <class GridLayerType, class MatrixType, class DiscreteFunctionType>
class MatrixSolveFunctor : public XT::Grid::Functor::Codim0<GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<GridLayerType> BaseType;

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
    matrices_[rhs_.space().grid_layer().indexSet().index(entity)].solve(local_solution, local_rhs);
    // write solution
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      solution_vector.set_entry(global_indices[ii], local_solution[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& rhs_;
  DiscreteFunctionType& solution_;
};

template <class GridLayerType, class MatrixType, class DiscreteFunctionType>
class MatrixApplyFunctor : public XT::Grid::Functor::Codim0<GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<GridLayerType> BaseType;

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
    matrices_[vector_.space().grid_layer().indexSet().index(entity)].mv(local_vector, local_result);

    // write solution
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      result_vector.set_entry(global_indices[ii], local_result[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& vector_;
  DiscreteFunctionType& result_;
};

template <class RhsEvaluationImp>
class AdvectionRhsOperator : public Dune::GDT::OperatorInterface<internal::AdvectionRhsOperatorTraits<RhsEvaluationImp>>
{
  //  static_assert(is_rhs_evaluation<RhsEvaluationImp>::value, "RhsEvaluationImp has to be derived from
  //  RhsInterface!");

public:
  typedef internal::AdvectionRhsOperatorTraits<RhsEvaluationImp> Traits;
  typedef typename Traits::RhsEvaluationType RhsEvaluationType;
  using BaseType = typename Dune::GDT::OperatorInterface<internal::AdvectionRhsOperatorTraits<RhsEvaluationImp>>;

  AdvectionRhsOperator(const RhsEvaluationType& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& /*param*/) const
  {
    std::fill(range.vector().begin(), range.vector().end(), 0);
    LocalVolumeIntegralFunctional<LocalFvRhsIntegrand<RhsEvaluationType, SourceType>,
                                  typename RangeType::SpaceType::BaseFunctionSetType>
        local_functional(rhs_evaluation_, source);
    VectorFunctionalBase<typename RangeType::VectorType,
                         typename RangeType::SpaceType,
                         typename RangeType::SpaceType::GridLayerType,
                         typename RangeType::DomainFieldType>
        functional_assembler(range.vector(), range.space());
    functional_assembler.append(local_functional);
    functional_assembler.assemble(true);
  }

  const RhsEvaluationType& evaluation() const
  {
    return rhs_evaluation_;
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixTraits>
  void assemble_jacobian(XT::LA::MatrixInterface<MatrixTraits, typename SourceType::RangeFieldType>& jac,
                         const SourceType& source,
                         const XT::Common::Parameter& /*param*/ = {}) const
  {
    typedef typename SourceType::SpaceType SpaceType;
    typedef typename SpaceType::BaseFunctionSetType BasisType;
    typedef LocalVolumeIntegralOperator<LocalFvRhsJacobianIntegrand<RhsEvaluationType, SourceType>, BasisType>
        LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source);
    SystemAssembler<SpaceType> assembler(source.space());
    assembler.append(local_operator, jac);
    assembler.assemble(true);
  }

  //  // assembles jacobian (jacobian is assumed to be zero initially)
  //  template <class SourceType, class MatrixType>
  //  void assemble_newton_matrix(std::vector<MatrixType>& newton_matrices,
  //                              const SourceType& source,
  //                              const XT::Common::Parameter& param) const
  //  {
  //    typedef LocalVolumeIntegralOperator<LocalFvRhsNewtonIntegrand<RhsEvaluationType, SourceType>> LocalOperatorType;
  //    LocalOperatorType local_operator(rhs_evaluation_, source, param);
  //    LocalVolumeTwoFormAssembler<LocalOperatorType> local_assembler(local_operator);
  //    SystemAssembler<typename SourceType::SpaceType> assembler(source.space());
  //    assembler.append(local_assembler, newton_matrices);
  //    assembler.assemble(false);
  //  }

  //  // solves with local jacobian on each entity
  //  template <class SourceType, class RangeType, class MatrixType>
  //  void solve(const std::vector<MatrixType>& newton_matrices,
  //             const SourceType& rhs,
  //             RangeType& solution,
  //             const XT::Common::Parameter& /*param*/ = {}) const
  //  {
  //    MatrixSolveFunctor<typename SourceType::SpaceType::GridLayerType, MatrixType, SourceType> functor(
  //        newton_matrices, rhs, solution);
  //    SystemAssembler<typename SourceType::SpaceType> assembler(rhs.space());
  //    assembler.append(functor);
  //    assembler.assemble(false);
  //  }

  //  // applies local jacobian on each entity
  //  template <class SourceType, class RangeType, class MatrixType>
  //  void mv(const std::vector<MatrixType>& newton_matrices,
  //          const SourceType& vector,
  //          RangeType& result,
  //          const XT::Common::Parameter& /*param*/ = {}) const
  //  {
  //    MatrixApplyFunctor<typename SourceType::SpaceType::GridLayerType, MatrixType, SourceType> functor(
  //        newton_matrices, vector, result);
  //    SystemAssembler<typename SourceType::SpaceType> assembler(vector.space());
  //    assembler.append(functor);
  //    assembler.assemble(false);
  //  }

private:
  const RhsEvaluationType& rhs_evaluation_;
}; // class AdvectionRhsOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_HH
