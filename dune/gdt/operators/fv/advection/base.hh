// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_BASE_HH
#define DUNE_GDT_OPERATORS_FV_BASE_HH

#include <utility> // for std::forward

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/grid/walker/apply-on.hh>

#include <dune/gdt/local/operators/fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/fv/boundary.hh>
#include <dune/gdt/operators/fv/quadrature.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp>
class AdvectionTraitsBase
{
  static_assert(XT::Functions::is_localizable_flux_function<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from LocalizableFluxFunctionInterface!");
  static_assert(is_localizable_boundary_value<BoundaryValueImp>::value,
                "BoundaryValueImp has to be derived from LocalizableBoundaryValueInterface!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueImp BoundaryValueType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef RangeFieldType FieldType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename AnalyticalFluxType::PartialURangeType JacobianType;
}; // class AdvectionTraitsBase


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueImp,
          class SourceImp,
          class RangeImp>
class AdvectionLocalizableDefault
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp> BaseType;

  static_assert(XT::Functions::is_localizable_flux_function<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_coupling_flux<NumericalCouplingFluxImp>::value,
                "NumericalCouplingFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_boundary_flux<NumericalBoundaryFluxImp>::value,
                "NumericalBoundaryFluxImp has to be derived from LocalNumericalBoundaryFluxInterface!");
  static_assert(is_localizable_boundary_value<BoundaryValueImp>::value,
                "BoundaryValueImp has to be derived from LocalizableBoundaryValueInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceImp>::value,
                "SourceImp has to be derived from LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueImp BoundaryValueType;
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
    this->append(
        local_operator_,
        new XT::Grid::ApplyOn::PartitionSetInnerIntersectionsPrimally<GridLayerType, Dune::Partitions::Interior>());
    this->append(local_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GridLayerType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridLayerType>());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


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
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> Quadrature1dType;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain - 1> IntersectionQuadratureType;

public:
  AdvectionOperatorBase(
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const IntersectionQuadratureType& intersection_quadrature = midpoint_quadrature<DomainFieldType, dimDomain>())
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , intersection_quadrature_(intersection_quadrature)
  {
  }

  template <class SourceType, class RangeType, class... LocalOperatorArgs>
  void apply(SourceType& source,
             RangeType& range,
             const XT::Common::Parameter& param,
             LocalOperatorArgs&&... local_operator_args) const
  {
    std::fill(range.vector().begin(), range.vector().end(), 0.);
    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                SourceType,
                                RangeType>
    localizable_operator(analytical_flux_,
                         boundary_values_,
                         source,
                         range,
                         param,
                         intersection_quadrature_,
                         std::forward<LocalOperatorArgs>(local_operator_args)...);
    localizable_operator.apply(true);
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const IntersectionQuadratureType intersection_quadrature_;
}; // class AdvectionOperatorBase<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_BASE_HH
