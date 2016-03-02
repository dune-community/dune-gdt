// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_HDD_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/checkerboard.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class SourceBeam : public TwoBeams<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef SourceBeam<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;
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

  static std::string static_id()
  {
    return BaseType::static_id() + ".sourcebeam";
  }

  std::string type() const
  {
    return BaseType::type() + ".sourcebeam";
  }

  static std::string short_id()
  {
    return "SourceBeam";
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
    using GetDataBaseType::basefunctions_values_at_plusone;
    using GetDataBaseType::precision;

    // q - (sigma_a + T/2*S*M^(-1))*u = Q(x)*base_integrated() - (sigma_a(x)*I_{nxn} + T(x)/2*S*M_inverse)*u = q(x) -
    // A(x)*u
    // sigma_a = 1 if x <= 2, 0 else
    // T = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
    // Q = 1 if 1 <= x <= 1.5, 0 else
    // Thus A(x) = I_{nxn}                                   if x <= 1
    //           = I_[nxn} + S*M_inverse                     if 1 < x <= 2
    //           = 5*S*M_inverse                             else (2 < x <= 3)
    // and  q(x) = base_integrated                           if 1 < x <= 1.5
    //           = 0                                         else
    // For Legendre Polynomials, l-th component of Source is -(sigma_a + T/2*l(l+1))*u[l] + \int_{-1}^1 Q*P_l d\mu), so
    // here
    // Source[l] = -u[l]                                     if x <= 1
    //           = -(1 + l(l+1))u[l] + 2*delta(l)            if 1 < x <= 1.5
    //           = -(1 + l(l+1))u[l]                         if 1.5 < x <= 2
    //           = -5*l(l+1)*u[l]                            else (2 < x <= 3)
    static void create_source_values(ConfigType& source_config)
    {
      if (exact_legendre()) {
        for (size_t ii = 0; ii < 6; ++ii) {
          std::string A_str = "[";
          std::string q_str = "[";
          for (size_t rr = 0; rr < dimRange; ++rr) {
            if (rr > 0) {
              A_str += "; ";
              q_str += " ";
            }
            if (ii == 2 && rr == 0) // 1 < x <= 1.5
              q_str += "2";
            else
              q_str += "0";
            for (size_t cc = 0; cc < dimRange; ++cc) {
              if (cc > 0)
                A_str += " ";
              if (ii == 0 || ii == 1) { // x <= 1
                if (cc == rr)
                  A_str += "-1";
                else
                  A_str += "0";
              } else if (ii == 2 || ii == 3) { // 1 <= x <= 2
                if (cc == rr)
                  A_str += DSC::to_string(-1.0 - cc * (cc + 1), precision);
                else
                  A_str += "0";
              } else { // 2 <= x <= 3
                if (cc == rr)
                  A_str += DSC::to_string(-5.0 * cc * (cc + 1), precision);
                else
                  A_str += "0";
              }
            }
          }
          A_str += "]";
          q_str += "]";
          source_config["A." + DSC::to_string(ii)] = A_str;
          source_config["b." + DSC::to_string(ii)] = q_str;
          //          std::cout << "A." << DSC::to_string(ii) <<  ": " << A_str << std::endl;
          //          std::cout << "q." << DSC::to_string(ii) <<  ": " << q_str << std::endl;
        }
      } else {
        MatrixType S_M_inverse(S());
        S_M_inverse.rightmultiply(M_inverse());
        for (size_t ii = 0; ii < 6; ++ii) {
          std::string A_str = "[";
          std::string q_str = "[";
          if (ii == 2) // 1 <= x <= 1.5
            q_str = DSC::to_string(base_integrated(), precision);
          for (size_t rr = 0; rr < dimRange; ++rr) {
            if (rr > 0) {
              A_str += "; ";
              if (ii != 2)
                q_str += " ";
            }
            if (ii != 2)
              q_str += "0";
            for (size_t cc = 0; cc < dimRange; ++cc) {
              if (cc > 0)
                A_str += " ";
              if (ii == 0 || ii == 1) { // x <= 1
                if (cc == rr)
                  A_str += "-1";
                else
                  A_str += "0";
              } else if (ii == 2 || ii == 3) { // 1 <= x <= 2
                if (cc == rr)
                  A_str += DSC::to_string(-1.0 - S_M_inverse[rr][cc], precision);
                else
                  A_str += DSC::to_string(-S_M_inverse[rr][cc], precision);
              } else { // 2 <= x <= 3
                A_str += DSC::to_string(-5.0 * S_M_inverse[rr][cc], precision);
              }
            }
          }
          A_str += "]";
          q_str += "]";
          source_config["A." + DSC::to_string(ii)] = A_str;
          source_config["b." + DSC::to_string(ii)] = q_str;
          //          std::cout << "A." << DSC::to_string(ii) <<  ": " << A_str << std::endl;
          //          std::cout << "q." << DSC::to_string(ii) <<  ": " << q_str << std::endl;
        }
      }
    } // ... create_source_values()

    // boundary value of kinetic equation is delta(v-1) at x = 0 and 10^(-4) at x = 3,
    // so k-th component of boundary value has to be 0.5*\phi_k(1) at x = 0 and 0.5*10^(-4)*base_integrated_k at x = 3
    // for Legendre polynomials, this is [0.5 0.5 0.5 ...] at x = 0 and [0.0002 0 0 ... ] at x = 3
    // simulate with linear interpolating function
    static std::string create_boundary_values()
    {
      if (exact_legendre()) {
        std::string str = "[";
        for (size_t cc = 0; cc < dimRange; ++cc) {
          if (cc > 0)
            str += " ";
          if (cc == 0)
            str += "0.5-0.4998*x[0]/3.0";
          else
            str += "0.5-0.5*x[0]/3.0";
        }
        str += "]";
        //        std::cout << str << std::endl;
        return str;
      } else {
        const auto& basefunctions_right = basefunctions_values_at_plusone();
        std::string str = "[";
        for (size_t cc = 0; cc < dimRange; ++cc) {
          if (cc > 0)
            str += " ";
          str += DSC::to_string(0.5 * basefunctions_right[cc], precision) + "+("
                 + DSC::to_string(0.5 * 0.0001 * (base_integrated()[cc]) - 0.5 * basefunctions_right[cc], precision)
                 + ")*x[0]/3.0";
        }
        str += "]";
        //        std::cout << str << std::endl;
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
    grid_config["upper_right"]  = "[3.0]";
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
  } // ... create(...)

  static ConfigType default_config(const std::string basefunctions_file = "", const std::string sub_name = "")
  {
    ConfigType config = BaseType::default_config(basefunctions_file, sub_name);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType source_config      = DefaultRHSType::default_config();
    source_config["lower_left"]   = "[0.0]";
    source_config["upper_right"]  = "[3.0]";
    source_config["num_elements"] = "[6]";
    GetData::create_source_values(source_config);
    source_config["name"] = static_id();
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

  SourceBeam(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> source_in,
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

#endif // DUNE_HDD_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
