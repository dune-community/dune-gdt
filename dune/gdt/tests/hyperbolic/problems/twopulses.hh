// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/common/string.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class TwoPulses
  : public TwoBeams< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef TwoPulses< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
  typedef TwoBeams< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >  BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultFunctionType;
  using typename BaseType::DefaultSourceType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::SourceType;
  using typename BaseType::FunctionType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".twopulses";
  }

  std::string type() const
  {
    return BaseType::type() + ".twopulses";
  }
private:
  template< size_t N >
  struct CreateSource {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        str += "0";
      }
      str += "]";
      return str;
    }
  };

  // boundary value has to be [50 50 50 ...]*exp(-(t-1)^2/2) at x = 0 and [50 -50 50 -50 ... ]**exp(-(t-1)^2/2) at x = 7
  // simulate with function(x) = 100*(x/7-0.5)*((x/7-0.5)/abs(x/7-0.5))^n*(-1)^n*exp(-(t-1)^2/2)
  template< size_t N >
  struct CreateBoundaryValues {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        str += "100*abs((x[0]/7-0.5))*((x[0]/7-0.5)^" + DSC::toString(cc%2) +")*((-1)^"+ DSC::toString(cc%2) + ")/(abs(x[0]/7-0.5)^" + DSC::toString(cc%2) + ")*exp(-((t-1)^2)/2)";
      }
      str += "]";
      return str;
    }
  };

protected:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[7.0]";
    grid_config["num_elements"] = "[1000]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

public:
  static std::unique_ptr< ThisType > create(const ConfigType cfg = default_config(),
                                            const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr< const DefaultFluxType > flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr< const DefaultSourceType > source(DefaultSourceType::create(config.sub("source")));
    const std::shared_ptr< const DefaultFunctionType > initial_values(DefaultFunctionType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr< const DefaultBoundaryValueType > boundary_values(DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique< ThisType >(flux, source, initial_values,
                                                  grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config();
    config.add(default_grid_config(), "grid", true);
    ConfigType source_config = DefaultSourceType::default_config();
    source_config["lower_left"] = "[0.0]";
    source_config["upper_right"] = "[7.0]";
    source_config["num_elements"] = "[1]";
    source_config["variable"] = "u";
    source_config["values.0"] = CreateSource< dimRange >::value_str();
    source_config["name"] = static_id();
    config.add(source_config, "source", true);
    ConfigType boundary_value_config = DefaultBoundaryValueType::default_config();
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = CreateBoundaryValues< dimRange >::value_str();
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values", true);
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  TwoPulses(const std::shared_ptr< const FluxType > flux,
           const std::shared_ptr< const SourceType > source,
           const std::shared_ptr< const FunctionType > initial_values,
           const ConfigType& grid_config,
           const ConfigType& boundary_info,
           const std::shared_ptr< const BoundaryValueType > boundary_values)
    : BaseType(flux,
               source,
               initial_values,
               grid_config,
               boundary_info,
               boundary_values)
  {}
};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
