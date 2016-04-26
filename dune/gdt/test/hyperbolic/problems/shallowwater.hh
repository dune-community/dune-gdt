// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH

#include <memory>

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class ShallowWater : public Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef ShallowWater<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::DefaultRHSType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".shallowwater";
  }

  std::string type() const override
  {
    return BaseType::type() + ".shallowwater";
  }

  static std::string short_id()
  {
    return "Shallowwater";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0]";
    grid_config["upper_right"]  = "[10.0]";
    grid_config["num_elements"] = "[10]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
  }

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["variable"]   = "u";
    flux_config["expression"] = "[u[1] u[1]*u[1]/u[0]+0.5*u[0]*u[0]]";
    flux_config["order"]      = "2";
    flux_config["gradient.0"] = "[0 1; -1.0*u[1]*u[1]/(u[0]*u[0])+u[0] 2*u[1]/u[0]]";
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"]   = "[0.0]";
    rhs_config["upper_right"]  = "[10.0]";
    rhs_config["num_elements"] = "[1]";
    rhs_config["variable"]     = "u";
    rhs_config["values.0"]     = "[0 0]";
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"]   = "[0.0]";
    initial_value_config["upper_right"]  = "[10.0]";
    initial_value_config["num_elements"] = "[5]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = "[1 0]";
    initial_value_config["values.1"]     = "[1 0]";
    initial_value_config["values.2"]     = "[1+((x[0]-4)^2)*((x[0]-6)^2)*exp(2-((x[0]-4)^2)-((x[0]-6)^2)) 0]";
    initial_value_config["values.3"]     = "[1 0]";
    initial_value_config["values.4"]     = "[1 0]";
    initial_value_config["order"] = "10";
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config    = DefaultBoundaryValueType::default_config();
    boundary_value_config["type"]       = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = "[0 0 0]";
    boundary_value_config["order"] = "0";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
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

  ShallowWater(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> rhs_in,
               const std::shared_ptr<const InitialValueType> initial_values_in, const ConfigType& grid_config_in,
               const ConfigType& boundary_info_in,
               const std::shared_ptr<const DefaultBoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 3;
  }
};


} // namespace Problems

// Test case for shallow water equations, see LeVeque, Finite Volume Methods for Hyperbolic Problems, 2002, Example 13.1
template <class G, class R = double>
class ShallowWaterTestCase
    : public Dune::GDT::Tests::NonStationaryTestCase<G, Problems::ShallowWater<typename G::template Codim<0>::Entity,
                                                                               typename G::ctype, G::dimension, R, 2>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d            = G::dimension;
  static const size_t dimRange     = 2;
  static const size_t dimRangeCols = 1;
  typedef typename Problems::ShallowWater<E, D, d, R, 2> ProblemType;

private:
  typedef typename Dune::GDT::Tests::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  ShallowWaterTestCase(const size_t num_refs = 2, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, Stuff::Grid::Providers::Cube<G>::create(ProblemType::default_grid_config())->grid_ptr(),
               num_refs)
    , problem_(*(ProblemType::create(ProblemType::default_config())))
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual std::bitset<d> periodic_directions() const override final
  {
    std::bitset<d> periodic_dirs;
    periodic_dirs.set();
    return periodic_dirs;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: Shallow Water                                           ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << "||  domain = [0, 10]                                                  ||\n"
        << "||  time = [0, " + DSC::toString(BaseType::t_end())
               + "]                                                   ||\n"
        << "||  flux = [u[1] u[1]*u[1]/u[0]+0.5*u[0]*u[0]]                        ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: solution on finest grid                       ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class ShallowWaterTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
