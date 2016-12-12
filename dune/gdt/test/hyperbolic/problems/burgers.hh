// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class Burgers : public Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef Burgers<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;
  typedef Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  using BaseType::dimDomain;
  using BaseType::dimRange;

  static std::string static_id()
  {
    return BaseType::static_id() + ".burgers";
  }

  std::string type() const override
  {
    return BaseType::type() + ".burgers";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[1.0 1.0 1.0]";
    grid_config["num_elements"] = "[8 8 8]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
  }

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config();
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType flux_config;
    flux_config["variable"] = "u";
    flux_config["expression"] = "[0.5*u[0]*u[0] 0.5*u[0]*u[0] 0.5*u[0]*u[0]]";
    flux_config["order"] = "2";
    if (dimDomain == 1)
      flux_config["gradient"] = "[u[0] 0 0]";
    else {
      flux_config["gradient.0"] = "[u[0] 0 0]";
      flux_config["gradient.1"] = "[u[0] 0 0]";
      flux_config["gradient.2"] = "[u[0] 0 0]";
    }
    config.add(flux_config, "flux", true);
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0 0.0 0.0]";
    initial_value_config["upper_right"] = "[1.0 1.0 1.0]";
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    if (dimDomain == 1)
      initial_value_config["values"] = "sin(pi*x[0])";
    else
      initial_value_config["values.0"] =
          "1.0/40.0*exp(1-(2*pi*x[0]-pi)*(2*pi*x[0]-pi)-(2*pi*x[1]-pi)*(2*pi*x[1]-pi))"; // bump, only in 2D or higher
    initial_value_config["name"] = static_id();
    initial_value_config["order"] = "10";
    config.add(initial_value_config, "initial_values", true);
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Burgers(const std::shared_ptr<const FluxType> flux,
          const std::shared_ptr<const RHSType> rhs,
          const std::shared_ptr<const InitialValueType> initial_values,
          const ConfigType& grid_config,
          const ConfigType& boundary_info,
          const std::shared_ptr<const BoundaryValueType> boundary_values)
    : BaseType(flux, rhs, initial_values, grid_config, boundary_info, boundary_values)
  {
  }

  virtual double CFL() const override
  {
    if (dimDomain == 1)
      return 0.5;
    else
      return 0.1;
  }

  virtual double t_end() const override
  {
    return 1.0;
  }
};


} // namespace Problems


template <class G, class R = double, size_t r = 1, size_t rC = 1>
class BurgersTestCase
    : public Dune::GDT::Test::NonStationaryTestCase<G,
                                                    Problems::Burgers<typename G::template Codim<0>::Entity,
                                                                      typename G::ctype,
                                                                      G::dimension,
                                                                      R,
                                                                      r,
                                                                      rC>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  static const size_t dimRange = r;
  static const size_t dimRangeCols = rC;
  typedef Problems::Burgers<E, D, d, R, r> ProblemType;

private:
  typedef Test::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;

  BurgersTestCase(const size_t num_refs = (d == 1 ? 4 : 1), const double divide_t_end_by = 1.0)
    : BaseType(
          divide_t_end_by, XT::Grid::make_cube_grid<GridType>(ProblemType::default_grid_config()).grid_ptr(), num_refs)
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
    const std::string domainstring = (d == 1)
                                         ? "||  domain = [0, 1]                                                   ||\n"
                                         : "||  domain = [0, 1] x [0, 1]                                          ||\n";
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: Burgers                                                 ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << domainstring
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                   ||\n"
        << "||  flux = 0.5*u[0]^2                                                 ||\n"
        << "||  rhs = 0                                                           ||\n"
        << "||  reference solution: discrete solution on finest grid              ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class BurgersTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_BURGERS_HH
