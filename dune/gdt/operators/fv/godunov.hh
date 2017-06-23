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
          class EigenSolverImp,
          class Traits>
class AdvectionGodunovOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class EigenSolverImp>
class AdvectionGodunovOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                  BoundaryValueFunctionImp,
                                                                  reconstructionOrder,
                                                                  slope_lim,
                                                                  realizability_lim,
                                                                  BasisFunctionImp,
                                                                  EigenSolverImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp,
                              BoundaryValueFunctionImp,
                              reconstructionOrder,
                              slope_lim,
                              realizability_lim,
                              BasisFunctionImp,
                              EigenSolverImp>
      BaseType;

public:
  using BaseType::polOrder;
  using BaseType::slope_limiter;
  using BaseType::realizability_limiting;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxType> NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType>
      NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp,
                                   BoundaryValueFunctionImp,
                                   polOrder,
                                   slope_limiter,
                                   realizability_limiting,
                                   BasisFunctionImp,
                                   EigenSolverImp,
                                   AdvectionGodunovOperatorTraits>
      derived_type;
}; // class AdvectionGodunovOperatorTraits


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
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::LocalfunctionType>,
          class Traits = internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp,
                                                                  BoundaryValueFunctionImp,
                                                                  polOrder,
                                                                  slope_lim,
                                                                  realizability_lim,
                                                                  BasisFunctionImp,
                                                                  EigenSolverImp>>
class AdvectionGodunovOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values)
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
