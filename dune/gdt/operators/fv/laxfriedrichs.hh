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

#ifndef DUNE_GDT_OPERATORS_FV_LAXFRIEDRICHS_HH
#define DUNE_GDT_OPERATORS_FV_LAXFRIEDRICHS_HH

#include <dune/gdt/local/fluxes/laxfriedrichs.hh>
#include <dune/gdt/operators/interfaces.hh>

#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/basisfunctions.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class Traits>
class AdvectionLaxFriedrichsOperator;


namespace internal {


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
  typedef typename Dune::GDT::LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxType, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                                               BoundaryValueType,
                                                                               LocalizableFunctionType>
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


} // namespace internal


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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_LAXFRIEDRICHS_HH
