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

#ifndef DUNE_GDT_OPERATORS_FV_MUSTA_HH
#define DUNE_GDT_OPERATORS_FV_MUSTA_HH

#include <dune/gdt/local/fluxes/musta.hh>
#include <dune/gdt/operators/interfaces.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#include "base.hh"

namespace Dune {
namespace GDT {

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          class EigenSolverImp,
          class RealizabilityLimiterImp,
          class Traits>
class AdvectionMustaOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class EigenSolverImp,
          class RealizabilityLimiterImp>
class AdvectionMustaOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                BoundaryValueImp,
                                                                reconstruction_order,
                                                                slope_lim,
                                                                EigenSolverImp,
                                                                RealizabilityLimiterImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueImp,
                              reconstruction_order,
                              slope_lim,
                              EigenSolverImp,
                              RealizabilityLimiterImp>
      BaseType;

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  typedef typename Dune::GDT::MustaLocalNumericalCouplingFlux<AnalyticalFluxType, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::MustaLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                                       BoundaryValueType,
                                                                       LocalizableFunctionType>
      NumericalBoundaryFluxType;
  typedef AdvectionMustaOperator<AnalyticalFluxImp,
                                 BoundaryValueImp,
                                 LocalizableFunctionImp,
                                 reconstruction_order,
                                 slope_lim,
                                 EigenSolverImp,
                                 RealizabilityLimiterImp,
                                 AdvectionMustaOperatorTraits>
      derived_type;
}; // class AdvectionMustaOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::RangeFieldType,
                                                    AnalyticalFluxImp::dimRange,
                                                    AnalyticalFluxImp::dimRangeCols>,
          class RealizabilityLimiterImp = NonLimitingRealizabilityLimiter<typename AnalyticalFluxImp::EntityType>,
          class Traits = internal::AdvectionMustaOperatorTraits<AnalyticalFluxImp,
                                                                BoundaryValueImp,
                                                                LocalizableFunctionImp,
                                                                polOrder,
                                                                slope_lim,
                                                                EigenSolverImp,
                                                                RealizabilityLimiterImp>>
class AdvectionMustaOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::Quadrature1dType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionMustaOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const bool is_linear = false,
                         const size_t num_stages = 2)
    : BaseType(analytical_flux, boundary_values, is_linear)
    , dx_(dx)
    , is_linear_(is_linear)
    , num_stages_(num_stages)
  {
  }

  AdvectionMustaOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const Quadrature1dType& quadrature_1d,
                         const std::shared_ptr<RealizabilityLimiterImp>& realizability_limiter = nullptr,
                         const bool is_linear = false,
                         const size_t num_stages = 2)
    : BaseType(analytical_flux, boundary_values, is_linear, quadrature_1d, realizability_limiter)
    , dx_(dx)
    , is_linear_(is_linear)
    , num_stages_(num_stages)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, num_stages_);
  }

private:
  const LocalizableFunctionType& dx_;
  const bool is_linear_;
  const size_t num_stages_;
}; // class AdvectionMustaOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_MUSTA_HH
