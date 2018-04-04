// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_LAXWENDROFF_HH
#define DUNE_GDT_OPERATORS_FV_LAXWENDROFF_HH

#include <dune/gdt/local/fluxes/laxwendroff.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          class LocalizableFunctionImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp,
          class Traits>
class AdvectionLaxWendroffOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          class LocalizableFunctionImp,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp>
class AdvectionLaxWendroffOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                      BoundaryValueImp,
                                                                      BoundaryInfoImp,
                                                                      reconstruction_order,
                                                                      slope_lim,
                                                                      RealizabilityLimiterImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueImp,
                              BoundaryInfoImp,
                              reconstruction_order,
                              slope_lim,
                              RealizabilityLimiterImp>
      BaseType;

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename Dune::GDT::LaxWendroffLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::LaxWendroffLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                                             BoundaryValueImp,
                                                                             BoundaryInfoImp,
                                                                             LocalizableFunctionType>
      NumericalBoundaryFluxType;
  typedef AdvectionLaxWendroffOperator<AnalyticalFluxImp,
                                       BoundaryValueImp,
                                       BoundaryInfoImp,
                                       LocalizableFunctionImp,
                                       reconstruction_order,
                                       slope_lim,
                                       RealizabilityLimiterImp,
                                       AdvectionLaxWendroffOperatorTraits>
      derived_type;
}; // class AdvectionLaxWendroffOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          class LocalizableFunctionImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          class RealizabilityLimiterImp = NonLimitingRealizabilityLimiter<typename AnalyticalFluxImp::EntityType>,
          class Traits = internal::AdvectionLaxWendroffOperatorTraits<AnalyticalFluxImp,
                                                                      BoundaryValueImp,
                                                                      BoundaryInfoImp,
                                                                      LocalizableFunctionImp,
                                                                      polOrder,
                                                                      slope_lim,
                                                                      RealizabilityLimiterImp>>
class AdvectionLaxWendroffOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::DomainType;
  using typename BaseType::OnedQuadratureType;
  using typename BaseType::RangeFieldType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionLaxWendroffOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const BoundaryInfoType& boundary_info,
                               const LocalizableFunctionType& dx,
                               const bool is_linear = false,
                               const RangeFieldType alpha = BoundaryValueImp::dimDomain)
    : BaseType(analytical_flux, boundary_values, boundary_info, is_linear)
    , dx_(dx)
    , is_linear_(is_linear)
    , alpha_(alpha)
  {
  }

  AdvectionLaxWendroffOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const BoundaryInfoType& boundary_info,
                               const LocalizableFunctionType& dx,
                               const OnedQuadratureType& quadrature_1d,
                               const std::shared_ptr<RealizabilityLimiterImp>& realizability_limiter = nullptr,
                               const bool is_linear = false,
                               const RangeFieldType alpha = BoundaryValueImp::dimDomain)
    : BaseType(analytical_flux, boundary_values, boundary_info, is_linear, quadrature_1d, realizability_limiter)
    , dx_(dx)
    , is_linear_(is_linear)
    , alpha_(alpha)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, is_linear_, alpha_);
  }

private:
  const LocalizableFunctionType& dx_;
  const bool is_linear_;
  const RangeFieldType alpha_;
}; // class AdvectionLaxWendroffOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_LAXWENDROFF_HH
