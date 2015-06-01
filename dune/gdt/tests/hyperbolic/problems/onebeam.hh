// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/grid/provider.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class OneBeam
  : public TwoBeams< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef OneBeam< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
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
    return BaseType::static_id() + ".onebeam";
  }

  std::string type() const
  {
    return BaseType::type() + ".onebeam";
  }
private:
  // precalculated boundary values, for n > 30, the boundary value is less than 1e-28 and thus set to 0;
  static const std::vector< RangeFieldImp > left_boundary_values_;

  static RangeFieldImp get_left_boundary_value(const size_t n) {
    if (n <= 30) {
      return left_boundary_values_[n];
    } else {
      return 0.0;
    }
  }

  // sigma_a = 10 if 0.4 <= x <= 0.7, 0 else
  // T = Q = 0
  // l-th component of Source is -(sigma_a + T/2*l(l+1))*u[l] + \int_{-1}^1 Q*P_l d\mu), so here
  // Source[l] = -10*u[l]       if 0.4 <= x <= 0.7
  //           = 0              else
  template< size_t N >
  struct CreateSourceValues {
    static void create(ConfigType& source_config)
    {
      for (size_t rr = 0; rr < 10; ++rr) {
        std::string str = "[";
        for (size_t cc = 0; cc < N; ++cc) {
          if (cc > 0)
            str += " ";
          if (rr >= 4 && rr <= 7)      // 0.4 <= x <= 0.7
            str += "-10.0*u[" + DSC::toString(cc) + "]";
          else
            str += "0.0";
        }
        str += "]";
        source_config["values." + DSC::toString(rr)] = str;
      }
    }
  };

  // boundary value has to be (l-th component) int_{-1}^1 3*exp(3*m + 3)/(exp(6) - 1) * P_l(m) dm at x = 0 and [0.0002 0 0 ... ] at x = 1
  template< size_t N >
  struct CreateBoundaryValues {
    static std::string value_str()
    {
      std::string str = "[";
      for (size_t cc = 0; cc < N; ++cc) {
          if (cc > 0)
            str += " ";
          if (cc == 0)
            str += DSC::toString(0.0002 - get_left_boundary_value(cc)) + "*x[0]+" + DSC::toString(get_left_boundary_value(cc));
          else
            str += DSC::toString(0.0 - get_left_boundary_value(cc)) + "*x[0]+" + DSC::toString(get_left_boundary_value(cc));
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
    ConfigType config = BaseType::default_config();
    config.add(default_grid_config(), "grid", true);
    ConfigType source_config = DefaultSourceType::default_config();
    source_config["lower_left"] = "[0.0]";
    source_config["upper_right"] = "[1.0]";
    source_config["num_elements"] = "[10]";
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

  OneBeam(const std::shared_ptr< const FluxType > flux,
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
}; // class OneBeam

template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
const std::vector< RangeFieldImp >
OneBeam< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >::left_boundary_values_ = {1.0,
                                                                                                  6.71636489980355855245e-01,
                                                                                                  3.28363510019644144755e-01,
                                                                                                  1.24363973280948905686e-01,
                                                                                                  3.81809056974300592424e-02,
                                                                                                  9.82125618865871928553e-03,
                                                                                                  2.16963300568142518701e-03,
                                                                                                  4.19513164039210756200e-04,
                                                                                                  7.20671854853712433829e-05,
                                                                                                  1.11324462887737092337e-05,
                                                                                                  1.56169232313775636906e-06,
                                                                                                  2.00600026809414359175e-07,
                                                                                                  2.37587842655794863414e-08,
                                                                                                  2.61015792958530426248e-09,
                                                                                                  2.67362899311749478277e-10,
                                                                                                  2.56499029050593908917e-11,
                                                                                                  2.31390262613575206647e-12,
                                                                                                  1.96974017566116948835e-13,
                                                                                                  1.58724211977210871905e-14,
                                                                                                  1.21415612755686868581e-15,
                                                                                                  8.83915394817958527742e-17,
                                                                                                  6.13842130565872485602e-18,
                                                                                                  4.07500767354120497201e-19,
                                                                                                  2.59097953469175454498e-20,
                                                                                                  1.58064025241231349219e-21,
                                                                                                  9.26712241830912425052e-23,
                                                                                                  5.22944129976226306897e-24,
                                                                                                  2.84427887291271140280e-25,
                                                                                                  1.49300327556249155579e-26,
                                                                                                  7.57264934397747012074e-28,
                                                                                                  3.71557124692250023583e-29};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH
