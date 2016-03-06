// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
#define DUNE_HDD_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH

#include <memory>

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>

#include <dune/gdt/test/nonstationary-eocstudy.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class ShallowWater
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > BaseType;
  typedef ShallowWater< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;

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
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[10.0]";
    grid_config["num_elements"] = "[100]";
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
    flux_config["variable"] = "u";
    flux_config["expression"] = "[u[1] u[1]*u[1]/u[0]+9.81*0.5*u[0]*u[0]]";
    flux_config["order"] = "2";
    flux_config["gradient.0"] = "[0 1; -1.0*u[1]*u[1]/(u[0]*u[0])+9.81*u[0] 2*u[1]/u[0]]";
    config.add(flux_config, "flux");
    ConfigType source_config;
    source_config["lower_left"] = "[0.0]";
    source_config["upper_right"] = "[10.0]";
    source_config["num_elements"] = "[1]";
    source_config["variable"] = "u";
    source_config["values.0"] = "[0 0]";
    source_config["name"] = static_id();
    config.add(source_config, "source");
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[10.0]";
    initial_value_config["num_elements"] = "[5]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = "[1 0]";
    initial_value_config["values.1"] = "[1 0]";
    initial_value_config["values.2"] = "[1+((x[0]-4)^2)*((x[0]-6)^2)*exp(2-((x[0]-4)^2)-((x[0]-6)^2)) 0]";
    initial_value_config["values.3"] = "[1 0]";
    initial_value_config["values.4"] = "[1 0]";
    initial_value_config["order"] = "10";
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config = DefaultBoundaryValueType::default_config();
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
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

  static std::unique_ptr< ThisType > create(const ConfigType cfg = default_config(),
                                            const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr< const DefaultFluxType > flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr< const DefaultRHSType > source(DefaultRHSType::create(config.sub("source")));
    const std::shared_ptr< const DefaultInitialValueType > initial_values(DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr< const DefaultBoundaryValueType > boundary_values(DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique< ThisType >(flux, source, initial_values,
                                                  grid_config, boundary_info, boundary_values);
  } // ... create(...)

  ShallowWater(const std::shared_ptr< const FluxType > flux_in,
               const std::shared_ptr< const RHSType > source_in,
               const std::shared_ptr< const InitialValueType > initial_values_in,
               const ConfigType& grid_config_in,
               const ConfigType& boundary_info_in,
               const std::shared_ptr< const DefaultBoundaryValueType > boundary_values_in)
    : BaseType(flux_in,
               source_in,
               initial_values_in,
               grid_config_in,
               boundary_info_in,
               boundary_values_in)
  {}
};

} // namespace Problems


template< class G, class R = double >
class ShallowWaterTestCase
  : public Dune::GDT::Tests::NonStationaryTestCase< G, Problems::ShallowWater< typename G::template Codim< 0 >::Entity,
                                                                     typename G::ctype, G::dimension, R, 2 > >
{
  typedef typename G::template Codim< 0 >::Entity E;
  typedef typename G::ctype D;
public:
  static const size_t d = G::dimension;
  static const size_t dimRange = 2;
  static const size_t dimRangeCols = 1;
  typedef typename Problems::ShallowWater< E, D, d, R, 2 > ProblemType;
private:
  typedef typename Dune::GDT::Tests::NonStationaryTestCase< G, ProblemType > BaseType;
public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  ShallowWaterTestCase(const size_t num_refs = 2)
    : BaseType(Stuff::Grid::Providers::Cube< G >::create(ProblemType::default_grid_config())->grid_ptr(), num_refs)
    , problem_(*(ProblemType::create(ProblemType::default_config())))
  {}

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual std::bitset< d > periodic_directions() const override final
  {
    std::bitset< d > periodic_dirs;
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
        << "||  flux = [u[1] u[1]*u[1]/u[0]+9.81*0.5*u[0]*u[0]]                   ||\n"
        << "||  source = 0                                                        ||\n"
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

#endif // DUNE_HDD_HYPERBOLIC_PROBLEMS_SHALLOWWATER_HH
