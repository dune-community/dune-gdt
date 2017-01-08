// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_RECTANGULARIC_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_RECTANGULARIC_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
class RectangularIC : public TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder>
{
  typedef RectangularIC<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> ThisType;
  typedef TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".rectangularic";
  }

  std::string type() const override
  {
    return BaseType::type() + ".rectangularic";
  }

protected:
  class GetData : BaseType::GetData
  {
    typedef typename BaseType::GetData GetDataBaseType;

  public:
    using GetDataBaseType::exact_legendre;
    using GetDataBaseType::S;
    using GetDataBaseType::M_inverse;
    using GetDataBaseType::base_integrated;

    // initial value of kinetic equation is 10 if 3 <= x <= 4 and 10^(-4) else, thus initial value of the
    // k-th component of the moment vector is 10*base_integrated_k or 10^(-4)*base_integrated_k.
    // For Legendre polynomials, this is 20 or 0.0002 if k == 0 and 0 else
    static void create_initial_values(ConfigType& initial_value_cfg)
    {
      if (exact_legendre()) {
        for (size_t ii = 0; ii < 7; ++ii) {
          std::string str = "[";
          for (size_t rr = 0; rr < dimRange; ++rr) {
            if (rr > 0)
              str += " ";
            if (rr == 0) {
              if (ii == 3)
                str += "20";
              else
                str += "0.0002";
            } else {
              str += "0.0";
            }
          }
          str += "]";
          std::string entry = "values." + Dune::XT::Common::to_string(ii);
          initial_value_cfg[entry] = str;
        }
      } else {
        for (size_t ii = 0; ii < 7; ++ii) {
          std::string str = "[";
          for (size_t rr = 0; rr < dimRange; ++rr) {
            if (rr > 0)
              str += " ";
            if (ii == 3)
              str += Dune::XT::Common::to_string(10 * base_integrated()[rr]);
            else
              str += Dune::XT::Common::to_string(0.0001 * base_integrated()[rr]);
          }
          str += "]";
          std::string entry = "values." + Dune::XT::Common::to_string(ii);
          initial_value_cfg[entry] = str;
        }
      }
    } // ... create_initial_values()

    // q - (sigma_a + T/2*S*M^(-1))*u = Q(x)*base_integrated() - (sigma_a(x)*I_{nxn} + T(x)/2*S*M_inverse)*u = q(x) -
    // A(x)*u
    // here, sigma_a = 0, T = 10^(-2) and Q = 0
    // Thus A(x) = 0.005*S*M_inverse and q(x) = 0
    // For Legendre Polynomials, this gives A[rr][rr] = 0.005*rr*(rr+1), A[rr][cc] = 0 if rr != cc;
    static void create_rhs_values(ConfigType& rhs_config)
    {
      if (exact_legendre()) {
        std::string A_str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            A_str += "; ";
          for (size_t cc = 0; cc < dimRange; ++cc) {
            if (cc > 0)
              A_str += " ";
            if (cc == rr)
              A_str += Dune::XT::Common::to_string(-0.005 * cc * (cc + 1));
            else
              A_str += "0";
          }
        }
        A_str += "]";
        rhs_config["A.0"] = A_str;
        rhs_config["b.0"] = Dune::XT::Common::to_string(RangeType(0));
      } else {
        MatrixType S_M_inverse(S());
        S_M_inverse.rightmultiply(M_inverse());
        S_M_inverse *= -0.005;
        rhs_config["A.0"] = Dune::XT::Common::to_string(S_M_inverse);
        rhs_config["b.0"] = Dune::XT::Common::to_string(RangeType(0));
      }
    } // ... create_rhs_values()

    // boundary value of kinetic equation is 10^(-4) at x = 0 and x = 7,
    // so k-th component of boundary value has to be 10^(-4)*base_integrated_k at x = 0 and x = 7.
    // For Legendre polynomials, this is [0.0002 0 0 ... ] at both sides
    // simulate with constant interpolating function
    static std::string create_boundary_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          if (rr == 0)
            str += "0.0002";
          else
            str += "0";
        }
        str += "]";
        return str;
      } else {
        ;
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += Dune::XT::Common::to_string(0.0001 * base_integrated()[rr]);
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
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[7.0]";
    grid_config["num_elements"] = "[500]";
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
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static std::unique_ptr<ThisType> create(const std::string basefunctions_file)
  {
    return create(default_config(basefunctions_file), static_id());
  } // ... create(...)

  static ConfigType default_config(const std::string basefunctions_file = "", const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config(basefunctions_file, sub_name);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType rhs_config;
    rhs_config["lower_left"] = "[0.0]";
    rhs_config["upper_right"] = "[7.0]";
    rhs_config["num_elements"] = "[1]";
    GetData::create_rhs_values(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs", true);
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[7.0]";
    initial_value_config["num_elements"] = "[7]";
    GetData::create_initial_values(initial_value_config);
    config.add(initial_value_config, "initial_values", true);
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
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

  RectangularIC(const std::shared_ptr<const FluxType> flux_in,
                const std::shared_ptr<const RHSType> rhs_in,
                const std::shared_ptr<const InitialValueType> initial_values_in,
                const ConfigType& grid_config_in,
                const ConfigType& boundary_info_in,
                const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 1.0;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return true;
  }
};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_RECTANGULARIC_HH
