// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "sourcebeam.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
class PlaneSourcePnLegendre
    : public SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder>
{
  typedef PlaneSourcePnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> ThisType;
  typedef SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
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
    return BaseType::static_id() + ".planesource";
  }

  std::string type() const override
  {
    return BaseType::type() + ".planesource";
  }

  static std::string short_id()
  {
    return "PlaneSource";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[-1.2]";
    grid_config["upper_right"] = "[1.2]";
    grid_config["num_elements"] = "[9600]";
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

  static ConfigType default_config(const RangeFieldImp psi_vac = 5e-9)
  {
    ConfigType config = BaseType::default_config(true, psi_vac);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType initial_value_config;
    create_initial_value_config(initial_value_config, default_grid_config());
    config.add(initial_value_config, "initial_values", true);
    ConfigType rhs_config;
    rhs_config["lower_left"] = default_grid_config()["lower_left"];
    rhs_config["upper_right"] = default_grid_config()["upper_right"];
    rhs_config["num_elements"] = "[1]";
    create_rhs_values_BGK(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs", true);
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = create_realizable_boundary_values(psi_vac);
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values", true);
    return config;
  } // ... default_config(...)

  PlaneSourcePnLegendre(const std::shared_ptr<const FluxType> flux_in,
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
  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  static void create_initial_value_config(ConfigType& initial_value_config,
                                          const ConfigType& grid_config = default_grid_config(),
                                          const RangeFieldImp psi_vac = 5e-9)
  {
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = grid_config["num_elements"];
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = static_id();
    const size_t num_elements = XT::Common::from_string<std::vector<size_t>>(grid_config["num_elements"])[0];
    const RangeFieldImp lower_left = XT::Common::from_string<std::vector<double>>(grid_config["lower_left"])[0];
    const RangeFieldImp upper_right = XT::Common::from_string<std::vector<double>>(grid_config["upper_right"])[0];
    const RangeFieldImp vol_entity = (upper_right - lower_left) / num_elements;
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    for (size_t ii = 0; ii < num_elements; ++ii) {
      std::string str = "[";
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (rr > 0)
          str += " ";
        if (ii == num_elements / 2 || ii == num_elements / 2 - 1) {
          // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
          if (rr == 0)
            str += XT::Common::to_string(2 * (psi_vac + 1. / (2. * vol_entity)), precision);
          else
            str += "0";
        } else {
          if (rr == 0)
            str += XT::Common::to_string(2 * psi_vac, precision);
          else
            str += "0";
        }
      } // rr
      str += "]";
      initial_value_config["values." + XT::Common::to_string(ii)] = str;
    } // ii
  } // ... create_initial_values()

  // n-th component of RHS is -sigma_t u_n + (sigma_s u_0 + 2Q) delta(n)
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // Thus, rhs[n] = (delta(n)-1)u[n]
  static void create_rhs_values_BGK(ConfigType& rhs_config)
  {
    Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
    Dune::FieldVector<RangeFieldImp, dimRange> b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      A[rr][rr] = (rr == 0) - 1;
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
  } // ... create_rhs_values_BGK()

  // Boundary value of kinetic equation is psi_vac at both boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
  // For Legendre-Polynomials, base_integrated_n = 2 delta(n).
  // Modell with constant function.
  static std::string create_realizable_boundary_values(const RangeFieldImp psi_vac = 5e-9)
  {
    std::string str = "[" + XT::Common::to_string(2. * psi_vac, precision);
    for (size_t cc = 1; cc < dimRange; ++cc)
      str += " 0";
    str += "]";
    return str;
  } // ... create_realizable_boundary_values()

}; // class PlaneSourcePnLegendre

/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t num_points>
class PlaneSourcePnHatFunctions
    : public SourceBeamPnHatFunctions<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points>
{
  typedef PlaneSourcePnHatFunctions<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points> ThisType;
  typedef SourceBeamPnHatFunctions<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
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
    return BaseType::static_id() + ".planesource_Pn_hatfunctions";
  }

  std::string type() const override
  {
    return BaseType::type() + ".planesource_Pn_hatfunctions";
  }

  static std::string short_id()
  {
    return "PlaneSourcePnHatFunctions";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[-1.2]";
    grid_config["upper_right"] = "[1.2]";
    grid_config["num_elements"] = "[240]";
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

  static ConfigType default_config(const RangeType& v_points = create_equidistant_points(),
                                   const RangeFieldImp psi_vac = 5e-9)
  {
    ConfigType config = BaseType::default_config(v_points, psi_vac);
    config.add(default_grid_config(), "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType initial_value_config;
    create_initial_value_config(initial_value_config, default_grid_config(), psi_vac, v_points);
    config.add(initial_value_config, "initial_values", true);
    ConfigType rhs_config;
    rhs_config["lower_left"] = default_grid_config()["lower_left"];
    rhs_config["upper_right"] = default_grid_config()["upper_right"];
    rhs_config["num_elements"] = "[1]";
    create_rhs_values(rhs_config, v_points);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs", true);
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = create_boundary_values(v_points, psi_vac);
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values", true);
    return config;
  } // ... default_config(...)

  PlaneSourcePnHatFunctions(const std::shared_ptr<const FluxType> flux_in,
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

  using BaseType::create_equidistant_points;

protected:
  using BaseType::hatfunctions_integrated;
  using BaseType::mass_matrix;
  using BaseType::tridiagonal_matrix_inverse;

  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  static void create_initial_value_config(ConfigType& initial_value_config,
                                          const ConfigType& grid_config = default_grid_config(),
                                          const RangeFieldImp psi_vac = 5e-9,
                                          const RangeType& v_points = create_equidistant_points())
  {
    const auto basis_integrated = hatfunctions_integrated(v_points);
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = grid_config["num_elements"];
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = static_id();
    const size_t num_elements = XT::Common::from_string<std::vector<size_t>>(grid_config["num_elements"])[0];
    const RangeFieldImp lower_left = XT::Common::from_string<std::vector<double>>(grid_config["lower_left"])[0];
    const RangeFieldImp upper_right = XT::Common::from_string<std::vector<double>>(grid_config["upper_right"])[0];
    const RangeFieldImp vol_entity = (upper_right - lower_left) / num_elements;
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    for (size_t ii = 0; ii < num_elements; ++ii) {
      std::string str = "[";
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (rr > 0)
          str += " ";
        // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
        if (ii == num_elements / 2 || ii == num_elements / 2 - 1)
          str += XT::Common::to_string(basis_integrated[rr] * (psi_vac + 1. / (2. * vol_entity)), precision);
        else
          str += XT::Common::to_string(basis_integrated[rr] * psi_vac, precision);
      } // rr
      str += "]";
      initial_value_config["values." + XT::Common::to_string(ii)] = str;
    } // ii
  } // ... create_initial_values()

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  static void create_rhs_values(ConfigType& rhs_config, const RangeType& v_points)
  {
    const RangeType integrated_basis = hatfunctions_integrated(v_points);
    const MatrixType M = mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);

    Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
    Dune::FieldVector<RangeFieldImp, dimRange> b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        A[rr][cc] = 0.5 * integrated_basis[rr] * c[cc] - (rr == cc);
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
  } // ... create_rhs_values()

  // Boundary value of kinetic equation is psi_vac at both boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
  // Modell with constant function.
  static std::string create_boundary_values(const RangeType& v_points, const RangeFieldImp psi_vac = 5e-9)
  {
    const RangeType integrated_basis = hatfunctions_integrated(v_points);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(integrated_basis[rr] * psi_vac, precision);
    }
    str += "]";
    return str;
  } // ... create_boundary_values()

}; // class PlaneSourcePnHatFunctions


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
