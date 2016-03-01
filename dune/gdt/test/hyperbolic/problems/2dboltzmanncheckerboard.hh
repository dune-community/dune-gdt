// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANNCHECKERBOARD_HH
#define DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANNCHECKERBOARD_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/gdt/playground/functions/entropymomentfunction.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/**
 * Testcase for the Boltzmann equation in two dimensions, see http://www.sciencedirect.com/science/article/pii/S0021999105002275?np=y
 * */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder >
class Boltzmann2DCheckerboard
  : public Boltzmann2DLineSource< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder >
{
  typedef Boltzmann2DCheckerboard< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder > ThisType;
  typedef Boltzmann2DLineSource< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder >  BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::RangeType;
  using typename BaseType::DomainType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".boltzmann2dcheckerboard";
  }

  std::string type() const
  {
    return BaseType::type() + ".boltzmann2dcheckerboard";
  }

  static std::string short_id()
  {
    return "Boltzmann2DCheckerBoard";
  }

protected:
  class GetData
      : public BaseType::GetData
  {
  public:
    using BaseType::GetData::precision; // precision for toString

    // source is q - (\Sigma_s \delta_{l0}\delta{m0} - \Sigma_t) * \psi_l^m
    // see the Checkerboard test case in http://www.sciencedirect.com/science/article/pii/S0021999105002275?np=y
    // q = 0 expect in the center where q = 1. Sigma_s = Sigma_t = 1 in scattering regions, Sigma_s = 0, Sigma_t
    // = 10 in absorbing regions. center is also a scattering region.
    static bool is_absorbing(const size_t row, const size_t col)
    {
      return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4)) || ((row == 3 || row == 5) && (col == 1 || col == 5));
    }

    static void create_source_values(ConfigType& source_config)
    {
      source_config["lower_left"] = "[0.0 0.0]";
      source_config["upper_right"] = "[7.0 7.0]";
      source_config["num_elements"] = "[7 7]";
      MatrixType S;
      RangeFieldImp Sigma_s;
      RangeFieldImp Sigma_t;
      RangeType q;
      for (size_t row = 0; row < 7; ++row) {
        for (size_t col = 0; col < 7; ++col) {
          if (row == 3 && col == 3) { // center
            q *= 0;
            q[0] = 1;
            Sigma_s = 1;
            Sigma_t = 1;
          } else if (is_absorbing(row, col)) { //absorbing regions
            q *= 0;
            Sigma_s = 0;
            Sigma_t = 10;
          } else { // scattering regions (without center)
            q *= 0;
            Sigma_s = 1;
            Sigma_t = 1;
          }
          S *= 0.0;
          S[pos(0, 0)][pos(0, 0)] = Sigma_s - Sigma_t;
          for (size_t l = 1; l <= momentOrder; ++l)
            for (size_t m = 0; m <= l; ++m)
              S[pos(l, m)][pos(l, m)] = -1.0*Sigma_t;
          size_t number = 7*row + col;
          source_config["A." + DSC::toString(number)] = DSC::toString(S, precision);
          source_config["b." + DSC::toString(number)] = DSC::toString(q);
          source_config["sparse." + DSC::toString(number)] = "true";
        }
      }
    } // ... create_source_values(...)

    // initial value is 0
    static std::string create_initial_values()
    {
      std::string str = "[";
      str += "0";
      for (size_t ii = 1; ii < dimRange; ++ii)
        str += " 0";
      str += "]";
      return str;
    } // ... create_initial_values()

  protected:
    using BaseType::GetData::pos;
  }; // class GetData

public:
  // Domain should be [-0.5, 0.5]^2, change it to [0, 1]^2 to make YaspGrid happy
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0 0.0]";
    grid_config["upper_right"] = "[7.0 7.0]";
    grid_config["num_elements"] = "[60 60]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

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

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config = BaseType::default_config().sub("flux");
    config.add(flux_config, "flux");
    ConfigType source_config = DefaultRHSType::default_config();
    GetData::create_source_values(source_config);
    config.add(source_config, "source");
    ConfigType initial_value_config = DefaultInitialValueType::default_config();
    initial_value_config["lower_left"] = "[0.0 0.0]";
    initial_value_config["upper_right"] = "[7.0 7.0]";
    initial_value_config["num_elements"] = "[1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = GetData::create_initial_values();
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config = BaseType::default_config().sub("boundary_values");
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Boltzmann2DCheckerboard(const std::shared_ptr< const FluxType > flux_in,
                        const std::shared_ptr< const RHSType > source_in,
                        const std::shared_ptr< const InitialValueType > initial_values_in,
                        const ConfigType& grid_config_in,
                        const ConfigType& boundary_info_in,
                        const std::shared_ptr< const BoundaryValueType > boundary_values_in)
    : BaseType(flux_in,
               source_in,
               initial_values_in,
               grid_config_in,
               boundary_info_in,
               boundary_values_in)
  {}
}; // ... Boltzmann2DCheckerboard ...


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANNCHECKERBOARD_HH
