// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_KINETIC_HH
#define DUNE_GDT_OPERATORS_FV_KINETIC_HH

#include <dune/gdt/local/fluxes/kinetic.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BasisfunctionImp,
          class GridLayerImp,
          size_t quadratureDim,
          size_t polOrder,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp,
          class Traits>
class AdvectionKineticOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BasisfunctionImp,
          class GridLayerImp,
          size_t quadratureDim,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp>
class AdvectionKineticOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  reconstruction_order,
                                                                  slope_lim,
                                                                  RealizabilityLimiterImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueImp,
                              reconstruction_order,
                              slope_lim,
                              RealizabilityLimiterImp>
      BaseType;

public:
  typedef typename Dune::GDT::
      KineticLocalNumericalCouplingFlux<AnalyticalFluxImp, BasisfunctionImp, GridLayerImp, quadratureDim>
          NumericalCouplingFluxType;
  typedef typename Dune::GDT::KineticLocalNumericalBoundaryFlux<AnalyticalFluxImp,
                                                                BoundaryValueImp,
                                                                BasisfunctionImp,
                                                                GridLayerImp,
                                                                quadratureDim>
      NumericalBoundaryFluxType;
  typedef AdvectionKineticOperator<AnalyticalFluxImp,
                                   BoundaryValueImp,
                                   BasisfunctionImp,
                                   GridLayerImp,
                                   quadratureDim,
                                   reconstruction_order,
                                   slope_lim,
                                   RealizabilityLimiterImp,
                                   AdvectionKineticOperatorTraits>
      derived_type;
}; // class AdvectionKineticOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BasisfunctionImp,
          class GridLayerImp,
          size_t quadratureDim = BasisfunctionImp::dimDomain,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          class RealizabilityLimiterImp = NonLimitingRealizabilityLimiter<typename AnalyticalFluxImp::EntityType>,
          class Traits = internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  BasisfunctionImp,
                                                                  GridLayerImp,
                                                                  quadratureDim,
                                                                  polOrder,
                                                                  slope_lim,
                                                                  RealizabilityLimiterImp>>
class AdvectionKineticOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::OnedQuadratureType;

  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const BasisfunctionImp& basis_functions,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, is_linear)
    , basis_functions_(basis_functions)
  {
  }

  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const BasisfunctionImp& basis_functions,
                           const OnedQuadratureType& quadrature_1d,
                           const std::shared_ptr<RealizabilityLimiterImp>& realizability_limiter = nullptr,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, is_linear, quadrature_1d, realizability_limiter)
    , basis_functions_(basis_functions)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, basis_functions_);
  }

private:
  const BasisfunctionImp& basis_functions_;
}; // class AdvectionKineticOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_KINETIC_HH
