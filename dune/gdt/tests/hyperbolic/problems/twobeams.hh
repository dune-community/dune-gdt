// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/affine.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class TwoBeams
  : public Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef TwoBeams< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
  typedef Default< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >  BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  typedef typename Dune::Stuff::Functions::Affine< FluxSourceEntityType,
                                                                 RangeFieldImp,
                                                                 dimRange,
                                                                 RangeFieldImp,
                                                                 dimRange,
                                                                 dimDomain >              DefaultFluxType;
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
    return BaseType::static_id() + ".twobeams";
  }

  std::string type() const
  {
    return BaseType::type() + ".twobeams";
  }
private:
  template< size_t N >
  struct CreateMatrix {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t rr = 0; rr < N; ++rr) {
        if (rr > 0)
          str += "; ";
        for (size_t cc = 0; cc < N; ++cc) {
          if (cc > 0)
            str += " ";
          if (cc == rr - 1)
            str += DSC::toString(double(rr)/(2.0*double(rr)+1.0));
          else if (cc == rr + 1)
            str += DSC::toString((double(rr)+1.0)/(2.0*double(rr)+1.0));
          else
            str += "0";
        }
      }
      str += "]";
      return str;
    }
  };

  template< size_t N >
  struct CreateSource {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        str += "-4*u[" + DSC::toString(cc) + "]";
      }
      str += "]";
      return str;
    }
  };

  template< size_t N >
  struct CreateInitialValues {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        if (cc == 0)
          str += "0.0002";
        else
          str += "0";
      }
      str += "]";
      return str;
    }
  };

  // boundary value has to be [50 50 50 ...] at x = 0 and [50 -50 50 -50 ... ] at x = 1
  // simulate with function(x) = 100*(x-0.5)*((x-0.5)/abs(x-0.5))^n*(-1)^n
  template< size_t N >
  struct CreateBoundaryValues {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        str += "100*abs((x[0]-0.5))*((x[0]-0.5)^" + DSC::toString(cc%2) +")*((-1)^"+ DSC::toString(cc%2) + ")/(abs(x[0]-0.5)^" + DSC::toString(cc%2) + ")";
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
    grid_config["upper_right"] = "[1.0]";
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
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config = DefaultFluxType::default_config();
    flux_config["type"] = DefaultFluxType::static_id();
    flux_config["A"] = CreateMatrix< dimRange >::value_str();
    config.add(flux_config, "flux");
    ConfigType source_config = DefaultSourceType::default_config();
    source_config["lower_left"] = "[0.0]";
    source_config["upper_right"] = "[1.0]";
    source_config["num_elements"] = "[1]";
    source_config["variable"] = "u";
    source_config["values.0"] = CreateSource< dimRange >::value_str();
    source_config["name"] = static_id();
    config.add(source_config, "source");
    ConfigType initial_value_config = DefaultFunctionType::default_config();
    initial_value_config["type"] = DefaultFunctionType::static_id();
    initial_value_config["variable"] = "x";
    initial_value_config["expression"] = CreateInitialValues< dimRange >::value_str();
    initial_value_config["order"] = "0";
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config = BoundaryValueType::default_config();
    boundary_value_config["type"] = BoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = CreateBoundaryValues< dimRange >::value_str();
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  TwoBeams(const std::shared_ptr< const FluxType > flux,
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

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
