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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH

#include <memory>

#include <dune/xt/functions/expression.hh>
#include <dune/xt/functions/checkerboard.hh>

#include <dune/gdt/local/fluxes/analytical.hh>
#include <dune/gdt/local/fluxes/rhs.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** TODO: replace default for initial values, RHS, boundary values (and flux?) by more general types once implemented
 *  TODO: choose correct SolutionType;
 * */
template <class ProblemImp, class E, class D, size_t d, class R, size_t r, size_t rC>
class Default : public ProblemInterface<E, D, d, R, r, rC>
{
  typedef Default<ProblemImp, E, D, d, R, r, rC> ThisType;
  typedef ProblemInterface<E, D, d, R, r, rC> BaseType;

protected:
  // we need an EntityType for the Expression functions that model q(u,x) and f(u). As we do not have a grid for the
  // u-variable, choose an arbitrary EntityType
  typedef typename Dune::template YaspGrid<r>::template Codim<0>::Entity DummyEntityType;
  typedef typename XT::Functions::ExpressionFunction<DummyEntityType, R, r, R, r, d> FluxExpressionFunctionType;
  typedef typename XT::Functions::ExpressionFunction<E, D, d, R, r, rC> InitialValueExpressionFunctionType;
  typedef typename XT::Functions::ExpressionFunction<DummyEntityType, R, r, R, r, rC> RHSExpressionFunctionType;

public:
  typedef ProblemImp ProblemType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxExpressionFunctionType, E, D, d, R, r, rC>
      DefaultFluxType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<InitialValueExpressionFunctionType, E, D, d, R, r, rC>
      DefaultInitialValueType;
  typedef typename BaseType::BoundaryValueType DefaultBoundaryValueType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<RHSExpressionFunctionType, E, D, d, R, r, rC>
      RHSCheckerboardFunctionType;
  typedef typename Dune::GDT::CheckerboardBasedRhsEvaluation<RHSCheckerboardFunctionType, E, D, d, R, r, rC>
      DefaultRHSType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return ProblemType::static_id();
  }

  static std::unique_ptr<ThisType> create(const ConfigType config = ProblemType::default_config())
  {
    const std::shared_ptr<const typename ProblemImp::DefaultFluxType> flux(
        ProblemImp::DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const typename ProblemImp::DefaultRHSType> rhs(
        ProblemImp::DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const typename ProblemImp::DefaultInitialValueType> initial_values(
        ProblemImp::DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const typename ProblemImp::DefaultBoundaryValueType> boundary_values(
        ProblemImp::DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ProblemType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  Default(const std::shared_ptr<const FluxType> flux_ptr,
          const std::shared_ptr<const RHSType> rhs_ptr,
          const std::shared_ptr<const InitialValueType> initial_values_ptr,
          const ConfigType& grid_cfg,
          const ConfigType& boundary_info_cfg,
          const std::shared_ptr<const BoundaryValueType> boundary_vals)
    : flux_(flux_ptr)
    , rhs_(rhs_ptr)
    , initial_values_(initial_values_ptr)
    , grid_config_(grid_cfg)
    , boundary_info_(boundary_info_cfg)
    , boundary_values_(boundary_vals)
  {
  }

  virtual const std::shared_ptr<const FluxType>& flux() const override
  {
    return flux_;
  }

  virtual const std::shared_ptr<const RHSType>& rhs() const override
  {
    return rhs_;
  }

  virtual const std::shared_ptr<const InitialValueType>& initial_values() const override
  {
    return initial_values_;
  }

  virtual const ConfigType grid_config() const override
  {
    return grid_config_;
  }

  virtual const ConfigType boundary_info() const override
  {
    return boundary_info_;
  }

  virtual const std::shared_ptr<const BoundaryValueType>& boundary_values() const override
  {
    return boundary_values_;
  }

  static double CFL()
  {
    return ProblemType::CFL();
  }

  static double t_end()
  {
    return ProblemType::t_end();
  }

  static bool has_non_zero_rhs()
  {
    return ProblemType::has_non_zero_rhs();
  }

private:
  const std::shared_ptr<const FluxType> flux_;
  const std::shared_ptr<const RHSType> rhs_;
  const std::shared_ptr<const InitialValueType> initial_values_;
  const ConfigType grid_config_;
  const ConfigType boundary_info_;
  const std::shared_ptr<const BoundaryValueType> boundary_values_;
}; // class Default<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH
