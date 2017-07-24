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

#ifndef DUNE_GDT_OPERATORS_FV_GODUNOV_HH
#define DUNE_GDT_OPERATORS_FV_GODUNOV_HH

#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          class EigenSolverImp,
          class RealizabilityLimiterImp,
          class Traits>
class AdvectionGodunovOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class EigenSolverImp,
          class RealizabilityLimiterImp>
class AdvectionGodunovOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  reconstruction_order,
                                                                  slope_lim,
                                                                  EigenSolverImp,
                                                                  RealizabilityLimiterImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueImp,
                              reconstruction_order,
                              slope_lim,
                              EigenSolverImp,
                              RealizabilityLimiterImp>
      BaseType;

public:
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxImp> NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp>
      NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp,
                                   BoundaryValueImp,
                                   reconstruction_order,
                                   slope_lim,
                                   EigenSolverImp,
                                   RealizabilityLimiterImp,
                                   AdvectionGodunovOperatorTraits>
      derived_type;
}; // class AdvectionGodunovOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::RangeFieldType,
                                                    AnalyticalFluxImp::dimRange,
                                                    AnalyticalFluxImp::dimRangeCols>,
          class RealizabilityLimiterImp = NonLimitingRealizabilityLimiter<typename AnalyticalFluxImp::EntityType>,
          class Traits = internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  polOrder,
                                                                  slope_lim,
                                                                  EigenSolverImp,
                                                                  RealizabilityLimiterImp>>
class AdvectionGodunovOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::Quadrature1dType;
  using typename BaseType::RangeFieldType;

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, is_linear)
    , is_linear_(is_linear)
  {
  }

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const Quadrature1dType& quadrature_1d,
                           const std::shared_ptr<RealizabilityLimiterImp>& realizability_limiter = nullptr,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, is_linear, quadrature_1d, realizability_limiter)
    , is_linear_(is_linear)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, is_linear_);
  }

private:
  const bool is_linear_;
}; // class AdvectionGodunovOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_GODUNOV_HH
