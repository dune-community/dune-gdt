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

#ifndef DUNE_GDT_OPERATORS_FV_KINETIC_HH
#define DUNE_GDT_OPERATORS_FV_KINETIC_HH

#include <dune/gdt/local/fluxes/kinetic.hh>
#include <dune/gdt/operators/interfaces.hh>

#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/basisfunctions.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t polOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class Traits>
class AdvectionKineticOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp>
class AdvectionKineticOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                  BoundaryValueFunctionImp,
                                                                  reconstructionOrder,
                                                                  slope_lim,
                                                                  realizability_lim,
                                                                  BasisFunctionImp>
{
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
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::BasisFunctionType;
  typedef typename Dune::GDT::KineticLocalNumericalCouplingFlux<AnalyticalFluxType> NumericalCouplingFluxType;
  typedef typename Dune::GDT::KineticLocalNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType>
      NumericalBoundaryFluxType;
  typedef AdvectionKineticOperator<AnalyticalFluxImp,
                                   BoundaryValueFunctionImp,
                                   polOrder,
                                   slope_limiter,
                                   realizability_limiting,
                                   BasisFunctionType,
                                   AdvectionKineticOperatorTraits>
      derived_type;
}; // class AdvectionKineticOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          bool realizability_lim = false,
          class BasisFunctionImp = Hyperbolic::Problems::HatFunctions<typename BoundaryValueFunctionImp::DomainFieldImp,
                                                                      BoundaryValueFunctionImp::dimDomain,
                                                                      typename BoundaryValueFunctionImp::RangeFieldType,
                                                                      BoundaryValueFunctionImp::dimRange,
                                                                      BoundaryValueFunctionImp::dimRangeCols>,
          class Traits = internal::AdvectionKineticOperatorTraits<AnalyticalFluxImp,
                                                                  BoundaryValueFunctionImp,
                                                                  polOrder,
                                                                  slope_lim,
                                                                  realizability_lim,
                                                                  BasisFunctionImp>>
class AdvectionKineticOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;

  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux, const BoundaryValueType& boundary_values)
    : BaseType(analytical_flux, boundary_values)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param);
  }
}; // class AdvectionKineticOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_KINETIC_HH
