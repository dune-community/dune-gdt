// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
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


template <class BasisfunctionImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class RectangularICPn : public KineticTransportEquation<BasisfunctionImp,
                                                        EntityImp,
                                                        DomainFieldImp,
                                                        dimDomain,
                                                        U_,
                                                        RangeFieldImp,
                                                        dimRange>
{
  typedef KineticTransportEquation<BasisfunctionImp, EntityImp, DomainFieldImp, dimDomain, U_, RangeFieldImp, dimRange>
      BaseType;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  RectangularICPn(const BasisfunctionType& basis_functions,
                  const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                  const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "pointsourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[-0.5 -0.5 -0.5]";
    grid_config["upper_right"] = "[0.5 0.5 0.5]";
    grid_config["num_elements"] = "[10 10 10]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{1}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{0.45}),
                                  std::make_pair("num_segments", std::vector<double>{1., 1., 1.})});
  }

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the moments is basis_integrated * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const FieldVector<size_t, 3> num_segments = get_num_segments(parameters());
    static const double sigma = 0.03;
    RangeType basis_integrated = basis_functions_.integrated();
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x) {
          auto ret = basis_integrated;
          ret *= psi_vac_ + 1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm() / (2. * sigma * sigma));
          return ret;
        },
        50);
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::get_num_segments;

  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
}; // class RectangularICPn<...>

template <class GridViewType,
          class BasisfunctionType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class RectangularICMn
    : public RectangularICPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
{
  typedef RectangularICPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
      BaseType;
  typedef RectangularICMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef GDT::EntropyBasedLocalFlux<GridViewType,
                                     BasisfunctionType,
                                     EntityType,
                                     DomainFieldType,
                                     dimDomain,
                                     U_,
                                     RangeFieldType,
                                     dimRange>
      ActualFluxType;
  using typename BaseType::QuadratureType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  RectangularICMn(const BasisfunctionType& basis_functions,
                  const QuadratureType& quadrature,
                  const GridViewType& grid_view,
                  const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                  const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
    , quadrature_(quadrature)
  {
  }

  static std::string static_id()
  {
    return "rectangularicmn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_view_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GridViewType& grid_view_;
  const QuadratureType& quadrature_;
}; // class RectangularICMn<...>

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

  static ConfigType default_config()
  {
    ConfigType config = BaseType::default_config();
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

protected:
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
};

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_RECTANGULARIC_HH
