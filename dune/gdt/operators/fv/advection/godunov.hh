// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_GODUNOV_HH
#define DUNE_GDT_OPERATORS_FV_GODUNOV_HH

#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp,
          class Traits>
class AdvectionGodunovOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp>
class AdvectionGodunovOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  BoundaryInfoImp,
                                                                  reconstruction_order,
                                                                  slope_lim,
                                                                  RealizabilityLimiterImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueImp,
                              BoundaryInfoImp,
                              reconstruction_order,
                              slope_lim,
                              RealizabilityLimiterImp>
      BaseType;

public:
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxImp> NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, BoundaryInfoImp>
          NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp,
                                   BoundaryValueImp,
                                   BoundaryInfoImp,
                                   reconstruction_order,
                                   slope_lim,
                                   RealizabilityLimiterImp,
                                   AdvectionGodunovOperatorTraits>
      derived_type;
}; // class AdvectionGodunovOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          class RealizabilityLimiterImp = NonLimitingRealizabilityLimiter<typename AnalyticalFluxImp::EntityType>,
          class Traits = internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  BoundaryInfoImp,
                                                                  polOrder,
                                                                  slope_lim,
                                                                  RealizabilityLimiterImp>>
class AdvectionGodunovOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::OnedQuadratureType;

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const BoundaryInfoType& boundary_info,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, boundary_info, is_linear)
    , is_linear_(is_linear)
  {
  }

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const BoundaryInfoType& boundary_info,
                           const OnedQuadratureType& quadrature_1d,
                           const bool regularize,
                           const std::shared_ptr<RealizabilityLimiterImp>& realizability_limiter = nullptr,
                           const bool is_linear = false)
    : BaseType(
          analytical_flux, boundary_values, boundary_info, is_linear, quadrature_1d, regularize, realizability_limiter)
    , is_linear_(is_linear)
  {
  }

  template <class SourceType, class RangeType>
  void apply(SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, is_linear_);
  }

private:
  const bool is_linear_;
}; // class AdvectionGodunovOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_GODUNOV_HH
