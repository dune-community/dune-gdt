// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
#define DUNE_HDD_HYPERBOLIC_PROBLEMS_TWOPULSES_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/common/string.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class TwoPulses : public TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef TwoPulses<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;
  typedef TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".twopulses";
  }

  std::string type() const
  {
    return BaseType::type() + ".twopulses";
  }

protected:
  class GetData : BaseType::GetData
  {
    typedef typename BaseType::GetData GetDataBaseType;

  public:
    using GetDataBaseType::exact_legendre;
    using GetDataBaseType::basefunctions_values_at_minusone;
    using GetDataBaseType::basefunctions_values_at_plusone;

    // q - (sigma_a + T/2*S*M^(-1))*u = Q(x)*base_integrated() - (sigma_a(x)*I_{nxn} + T(x)/2*S*M_inverse)*u = q(x) -
    // A(x)*u
    // here, sigma_a = 0, T = 0 and Q = 0
    // Thus A(x) = 0 and q(x) = 0
    static void create_source_values(ConfigType& source_config)
    {
      source_config["A.0"] = DSC::to_string(MatrixType(0));
      source_config["b.0"] = DSC::to_string(RangeType(0));
    } // ... create_source_values()

    // boundary value of kinetic equation is 100*delta(v-1)**exp(-(t-1)^2/2) at x = 0 and 100*delta(v+1)*exp(-(t-1)^2/2)
    // at x = 7, so k-th component of boundary value has to be 50*\phi_k(1)*exp(-(t-1)^2/2) at x = 0 and
    // 50*\phi_k(-1)*exp(-(t-1)^2/2) at x = 7.
    // Simulate with function(x) = 50*((\phi_k(-1) - \phi_k(1))*x/7 + \phi_k(1))*exp(-(t-1)^2/2),
    // For Legendre polynomials, this is [50 50 50 ...]*exp(-(t-1)^2/2) at x = 0 and
    // [50 -50 50 -50 ... ]*exp(-(t-1)^2/2) at x = 7
    // simulate with function(x) = 50*((-1)^n - 1)*x/7 + 1)*exp(-(t-1)^2/2)
    static std::string create_boundary_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += "50*(" + DSC::to_string(((1.0 - 2.0 * (rr % 2)) - 1.0)) + "*x[0]/7.0+1)*exp((-(t-1)^2)/2)";
        }
        str += "]";
        return str;
      } else {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += "50*(" + DSC::to_string(basefunctions_values_at_minusone()[rr] - basefunctions_values_at_plusone()[rr])
                 + "*x[0]/7.0+" + DSC::to_string(basefunctions_values_at_plusone()[rr]) + ")*exp((-(t-1)^2)/2)";
        }
        str += "]";
        return str;
      }
    } // ... create_boundary_values()
  }; // class GetData

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0]";
    grid_config["upper_right"]  = "[7.0]";
    grid_config["num_elements"] = "[1000]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const DefaultRHSType> source(DefaultRHSType::create(config.sub("source")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(
        flux, source, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static std::unique_ptr<ThisType> create(const std::string basefunctions_file)
  {
    return create(default_config(basefunctions_file), static_id());
  }

  static ConfigType default_config(const std::string basefunctions_file, const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config(basefunctions_file, sub_name);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType source_config      = DefaultRHSType::default_config();
    source_config["lower_left"]   = "[0.0]";
    source_config["upper_right"]  = "[7.0]";
    source_config["num_elements"] = "[1]";
    GetData::create_source_values(source_config);
    config.add(source_config, "source", true);
    ConfigType boundary_value_config    = DefaultBoundaryValueType::default_config();
    boundary_value_config["type"]       = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = GetData::create_boundary_values();
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

  TwoPulses(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> source_in,
            const std::shared_ptr<const InitialValueType> initial_values_in, const ConfigType& grid_config_in,
            const ConfigType& boundary_info_in, const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, source_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }
};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_HDD_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
