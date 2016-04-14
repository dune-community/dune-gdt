// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/common/string.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class OneBeam : public TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef OneBeam<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;
  typedef TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;

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

  static std::string static_id()
  {
    return BaseType::static_id() + ".onebeam";
  }

  std::string type() const override
  {
    return BaseType::type() + ".onebeam";
  }

protected:
  class GetData : BaseType::GetData
  {
    typedef typename BaseType::GetData GetDataBaseType;

  public:
    using GetDataBaseType::exact_legendre;
    using GetDataBaseType::base_integrated;
    using GetDataBaseType::onebeam_left_boundary_values;

    // q - (sigma_a + T/2*S*M^(-1))*u = Q(x)*base_integrated() - (sigma_a(x)*I_{nxn} + T(x)/2*S*M_inverse)*u = q(x) -
    // A(x)*u
    // here, sigma_a = 10 if 0.4 <= x <= 0.7, 0 else
    // T = Q = 0
    static void create_rhs_values(ConfigType& rhs_config)
    {
      for (size_t ii = 0; ii < 10; ++ii) {
        std::string A_str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            A_str += "; ";
          for (size_t cc = 0; cc < dimRange; ++cc) {
            if (cc > 0)
              A_str += " ";
            if (rr == cc && ii >= 4 && ii <= 6)
              A_str += "-10";
            else
              A_str += "0";
          }
        }
        A_str += "]";
        rhs_config["A." + DSC::to_string(ii)] = A_str;
        rhs_config["b." + DSC::to_string(ii)] = DSC::to_string(RangeType(0));
      }
    } // ... create_rhs_values()

    // boundary value of kinetic equation is 3*exp(3*v + 3)/(exp(6) - 1) at x = 0 and 10^(-4) at x = 1,
    // so k-th component of boundary value has to be (3*exp(3*v + 3)/(exp(6) - 1), \phi_k(v))_v at x = 0 and
    // 10^(-4)*base_integrated_k at x = 1.
    // simulate with linear interpolating function
    static std::string create_boundary_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          if (rr == 0)
            str += DSC::to_string(0.0002 - get_left_boundary_value(rr)) + "*x[0]+"
                   + DSC::to_string(get_left_boundary_value(rr));
          else
            str += DSC::to_string(0.0 - get_left_boundary_value(rr)) + "*x[0]+"
                   + DSC::to_string(get_left_boundary_value(rr));
        }
        str += "]";
        return str;
      } else {
        std::string str = "[";
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (rr > 0)
            str += " ";
          str += DSC::to_string(0.0001 * base_integrated()[rr] - onebeam_left_boundary_values()[rr]) + "*x[0]+"
                 + DSC::to_string(onebeam_left_boundary_values()[rr]);
        }
        str += "]";
        return str;
      }
    } // ... create_boundary_values()

  private:
    // precalculated boundary values if the legendre polynomials are used. For n > 30, the boundary value is less than
    // 1e-28 and thus set to 0
    static const std::vector<RangeFieldImp> left_boundary_values_;

    static RangeFieldImp get_left_boundary_value(const size_t n)
    {
      if (n <= 30) {
        return left_boundary_values_[n];
      } else {
        return 0.0;
      }
    }
  }; // class GetData

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0]";
    grid_config["upper_right"]  = "[1.0]";
    grid_config["num_elements"] = "[25]";
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
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static std::unique_ptr<ThisType> create(const std::string basefunctions_file)
  {
    return create(default_config(basefunctions_file), static_id());
  }

  static ConfigType default_config(const std::string basefunctions_file = "", const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config(basefunctions_file, sub_name);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType rhs_config;
    rhs_config["lower_left"]   = "[0.0]";
    rhs_config["upper_right"]  = "[1.0]";
    rhs_config["num_elements"] = "[10]";
    rhs_config["variable"] = "u";
    GetData::create_rhs_values(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs", true);
    ConfigType boundary_value_config;
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

  OneBeam(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> rhs_in,
          const std::shared_ptr<const InitialValueType> initial_values_in, const ConfigType& grid_config_in,
          const ConfigType& boundary_info_in, const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 4.0;
  }

  virtual bool is_linear() const override
  {
    return true;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return true;
  }
}; // class OneBeam

template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
const std::vector<RangeFieldImp> OneBeam<EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                         rangeDim>::GetData::left_boundary_values_ = {1.0,
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
