// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH

#include <ostream>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/default.hh>
#include <dune/xt/functions/expression.hh>
#include <dune/xt/functions/checkerboard.hh>

#include <dune/gdt/local/fluxes/interfaces.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {


/** Interface for problem of the form delta_t u + div f(u,x,t) = q(u,x,t) where u: R^d \to R^{r x rC}.
 *  TODO: implement for non-autonomous fluxes.
 *  TODO: implement for rangeDimCols > 1.
 *  TODO: replace TimeDependentExpression by a parametric function interface in dune-xt, once it is available.
 *  TODO: think about SolutionType (remove? use another type?)
 * */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class ProblemInterface
{
  typedef ProblemInterface<E, D, d, R, r, rC> ThisType;

public:
  typedef E EntityType;
  typedef D DomainFieldType;
  static const size_t dimDomain = d;
  typedef R RangeFieldType;
  static const size_t dimRange     = r;
  static const size_t dimRangeCols = rC;

  typedef Dune::GDT::AutonomousAnalyticalFluxInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                                       dimRangeCols>
      FluxType;
  typedef Dune::GDT::RhsEvaluationFluxInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                                dimRangeCols>
      RHSType;
  typedef Dune::XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                                    dimRangeCols>
      InitialValueType;
  typedef typename Dune::XT::Functions::TimeDependentExpressionFunction<EntityType, DomainFieldType, dimDomain,
                                                                   RangeFieldType, dimRange, dimRangeCols, double>
      BoundaryValueType;
  typedef Dune::XT::Common::Configuration ConfigType;
  typedef Dune::XT::Functions::TimeDependentFunctionInterface<
      typename Dune::XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>>
      SolutionType;

  virtual ~ProblemInterface()
  {
  }

  static std::string static_id()
  {
    return "gdt.hyperbolic.problem";
  }

  virtual std::string type() const
  {
    return "gdt.hyperbolic.problem";
  }

  virtual const std::shared_ptr<const FluxType>& flux() const = 0;

  virtual const std::shared_ptr<const RHSType>& rhs() const = 0;

  virtual const std::shared_ptr<const InitialValueType>& initial_values() const = 0;

  virtual const ConfigType grid_config() const = 0;

  virtual const ConfigType boundary_info() const = 0;

  virtual const std::shared_ptr<const BoundaryValueType>& boundary_values() const = 0;

  virtual double CFL() const
  {
    return 0.5;
  }

  virtual double t_end() const
  {
    return 1.0;
  }

  virtual bool has_non_zero_rhs() const
  {
    return false;
  }
}; // ProblemInterface


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH
