// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH

#include <memory>

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>

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
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class Default : public ProblemInterface<E, D, d, R, r, rC>
{
  typedef ProblemInterface<E, D, d, R, r, rC> BaseType;
  typedef Default<E, D, d, R, r, rC> ThisType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  // we need an EntityType for the Expression functions that model q(u,x) and f(u). As we do not have a grid for the
  // u-variable,
  // choose an arbitrary EntityType
  typedef typename Dune::template YaspGrid<r>::template Codim<0>::Entity DummyEntityType;
  typedef typename DS::Functions::Expression<DummyEntityType, R, r, R, r, d> FluxExpressionFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxExpressionFunctionType, E, D, d, R, r, rC>
      DefaultFluxType;
  typedef typename DS::Functions::Expression<E, D, d, R, r, rC> InitialValueExpressionFunctionType;
  typedef typename Dune::Stuff::Functions::FunctionCheckerboard<InitialValueExpressionFunctionType, E, D, d, R, r, rC>
      DefaultInitialValueType;
  typedef typename BaseType::BoundaryValueType DefaultBoundaryValueType;
  typedef typename DS::Functions::Expression<DummyEntityType, R, r, R, r, rC> RHSExpressionFunctionType;
  typedef typename DS::Functions::FunctionCheckerboard<RHSExpressionFunctionType, E, D, d, R, r, rC>
      RHSCheckerboardFunctionType;
  typedef
      typename Dune::GDT::CheckerboardBasedRhsEvaluationFluxInterface<RHSCheckerboardFunctionType, E, D, d, R, r, rC>
          DefaultRHSType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".default";
  }

  virtual std::string type() const override
  {
    return BaseType::type() + ".default";
  }

private:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0 0.0 0.0]";
    grid_config["upper_right"]  = "[1.0 1.0 1.0]";
    grid_config["num_elements"] = "[512 60 60]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
  }

public:
  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"]       = DefaultFluxType::static_id();
    flux_config["variable"]   = "u";
    flux_config["expression"] = "[0 0 0]";
    flux_config["order"]      = "0";
    flux_config["gradient"]   = "[0 0 0; 0 0 0; 0 0 0]";
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"]   = "[0.0 0.0 0.0]";
    rhs_config["upper_right"]  = "[1.0 1.0 1.0]";
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["variable"]     = "u";
    rhs_config["values"]       = "[0]";
    rhs_config["name"]         = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["type"]       = InitialValueType::static_id();
    initial_value_config["variable"]   = "x";
    initial_value_config["expression"] = "[0 0 0]";
    initial_value_config["order"]      = "0";
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"]       = BoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = "[0 0 0]";
    boundary_value_config["order"]      = "1";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const ConfigType cfg       = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  Default(const std::shared_ptr<const FluxType> flux_ptr, const std::shared_ptr<const RHSType> rhs_ptr,
          const std::shared_ptr<const InitialValueType> initial_values_ptr =
              std::make_shared<DefaultInitialValueType>("x", "[0 0 0]", 0),
          const ConfigType& grid_cfg          = default_grid_config(),
          const ConfigType& boundary_info_cfg = default_boundary_info_config(),
          const std::shared_ptr<const DefaultBoundaryValueType> boundary_vals =
              std::make_shared<DefaultBoundaryValueType>("x", "[0 0 0]", 0))
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

  virtual double CFL() const override
  {
    return 0.5;
  }

  virtual double t_end() const override
  {
    return 1.0;
  }

  virtual bool is_linear() const override
  {
    return false;
  }

private:
  const std::shared_ptr<const FluxType> flux_;
  const std::shared_ptr<const RHSType> rhs_;
  const std::shared_ptr<const InitialValueType> initial_values_;
  const ConfigType grid_config_;
  const ConfigType boundary_info_;
  const std::shared_ptr<const BoundaryValueType> boundary_values_;
};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_DEFAULT_HH
