// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <memory>
#include <vector>
#include <string>>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/checkerboard.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class SourceBeam
  : public TwoBeams< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef SourceBeam< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
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
    return BaseType::static_id() + ".sourcebeam";
  }

  std::string type() const
  {
    return BaseType::type() + ".sourcebeam";
  }
private:
  // sigma_a = 1 if x <= 2, 0 else
  // T = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  // l-th component of Source is -(sigma_a + T/2*l(l+1))*u[l] + \int_{-1}^1 Q*P_l d\mu), so here
  // Source[l] = -u[l]                                     if x <= 1
  //           = -(1 + l(l+1))u[l] + 2*delta(l)            if 1 < x <= 1.5
  //           = -(1 + l(l+1))u[l]                         if 1.5 < x <= 2
  //           = -5*l(l+1)*u[l]                            else (2 < x <= 3)
  template< size_t N >
  struct CreateSourceValues {
    static void create(ConfigType& source_config)
    {
      for (size_t rr = 0; rr < 6; ++rr) {
        std::string str = "[";
        for (size_t cc = 0; cc < N; ++cc) {
          if (cc > 0)
            str += " ";
          if (rr == 0 || rr == 1) {                     // x <= 1
            str += "-u[" + DSC::toString(cc) + "]";
          } else if (rr == 2) {                         // 1 <= x <= 1.5
            if (cc == 0)
              str += "-(1+" + DSC::toString(cc*(cc+1)) + ")*u[" + DSC::toString(cc) + "]+2";
            else
              str += "-(1+" + DSC::toString(cc*(cc+1)) + ")*u[" + DSC::toString(cc) + "]";
          } else if (rr == 3) {                         // 1.5 <= x <= 2
            str += "-(1+" + DSC::toString(cc*(cc+1)) + ")*u[" + DSC::toString(cc) + "]";
          } else {                                      // 2 <= x <= 3
            str += "-5*" + DSC::toString(cc*(cc+1)) + "*u[" + DSC::toString(cc) + "]";
          }
        }
        str += "]";
        source_config["values." + DSC::toString(rr)] = str;
      }
    }
  };

  // boundary value has to be [0.5 0.5 0.5 ...] at x = 0 and [0.0002 0 0 ... ] at x = 3
  // simulate with function(x) = 0.5 - (0.5-0.0002)*x/3 for first component and 0.5 - 0.5*x/3 for other components
  template< size_t N >
  struct CreateBoundaryValues {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
        if (cc > 0)
          str += " ";
        if (cc == 0)
         str += "0.5-0.4998*x[0]/3.0";
        else
          str += "0.5-0.5*x[0]/3.0";
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
    grid_config["upper_right"] = "[3.0]";
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
    source_config["upper_right"] = "[1.0]";
    source_config["num_elements"] = "[6]";
    source_config["variable"] = "u";
    CreateSourceValues< dimRange >::create(source_config);
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

  SourceBeam(const std::shared_ptr< const FluxType > flux,
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

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
