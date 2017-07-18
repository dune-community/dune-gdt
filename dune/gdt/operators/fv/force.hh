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

#ifndef DUNE_GDT_OPERATORS_FV_FORCE_HH
#define DUNE_GDT_OPERATORS_FV_FORCE_HH

#include <dune/gdt/local/fluxes/force.hh>
#include <dune/gdt/operators/interfaces.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

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
          class EigenSolverImp,
          class Traits>
class AdvectionForceOperator;


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class EigenSolverImp>
class AdvectionForceOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp,
                                                                BoundaryValueFunctionImp,
                                                                reconstructionOrder,
                                                                slope_lim,
                                                                realizability_lim,
                                                                BasisFunctionImp,
                                                                EigenSolverImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

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
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::BasisFunctionType;
  typedef typename Dune::GDT::ForceLocalNumericalCouplingFlux<AnalyticalFluxType, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::ForceLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType,
                                                                       BoundaryValueType,
                                                                       LocalizableFunctionType>
      NumericalBoundaryFluxType;
  typedef AdvectionForceOperator<AnalyticalFluxImp,
                                 BoundaryValueFunctionImp,
                                 LocalizableFunctionImp,
                                 polOrder,
                                 slope_limiter,
                                 realizability_limiting,
                                 BasisFunctionType,
                                 EigenSolverImp,
                                 AdvectionForceOperatorTraits>
      derived_type;
}; // class AdvectionForceOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t polOrder = 0,
          SlopeLimiters slope_lim = SlopeLimiters::minmod,
          bool realizability_lim = false,
          class BasisFunctionImp =
              Hyperbolic::Problems::HatFunctions<typename BoundaryValueFunctionImp::DomainFieldType,
                                                 BoundaryValueFunctionImp::dimDomain,
                                                 typename BoundaryValueFunctionImp::RangeFieldType,
                                                 BoundaryValueFunctionImp::dimRange,
                                                 BoundaryValueFunctionImp::dimRangeCols>,
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::RangeFieldType,
                                                    AnalyticalFluxImp::dimRange,
                                                    AnalyticalFluxImp::dimRangeCols>,
          class Traits = internal::AdvectionForceOperatorTraits<AnalyticalFluxImp,
                                                                BoundaryValueFunctionImp,
                                                                LocalizableFunctionImp,
                                                                polOrder,
                                                                slope_lim,
                                                                realizability_lim,
                                                                BasisFunctionImp,
                                                                EigenSolverImp>>
class AdvectionForceOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionForceOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const bool is_linear = false)
    : BaseType(analytical_flux, boundary_values, is_linear)
    , dx_(dx)
    , is_linear_(is_linear)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, is_linear_);
  }

private:
  const LocalizableFunctionType& dx_;
  const bool is_linear_;
}; // class AdvectionForceOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_FORCE_HH
