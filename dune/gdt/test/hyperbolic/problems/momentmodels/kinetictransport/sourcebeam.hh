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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/pdelab/common/crossproduct.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/math.hh>

#include "kinetictransportequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp,
          class GridLayerImp,
          class EntityImp,
          class DomainFieldImp,
          size_t dimDomain,
          class U_,
          class RangeFieldImp,
          size_t dimRange>
class SourceBeamPn : public KineticTransportEquation<BasisfunctionImp,
                                                     GridLayerImp,
                                                     EntityImp,
                                                     DomainFieldImp,
                                                     dimDomain,
                                                     U_,
                                                     RangeFieldImp,
                                                     dimRange>
{
  typedef KineticTransportEquation<BasisfunctionImp,
                                   GridLayerImp,
                                   EntityImp,
                                   DomainFieldImp,
                                   dimDomain,
                                   U_,
                                   RangeFieldImp,
                                   dimRange>
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
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  SourceBeamPn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const QuadratureType& quadrature = default_quadrature();
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer_, quadrature_, 6, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "sourcebeampn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[3.0]";
    grid_config["num_elements"] = "[300]";
    grid_config["overlap_size"] = "[1]";
    grid_config["num_quad_cells"] = "[50]";
    grid_config["quad_order"] = "50";
    return grid_config;
  }

  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{1, 1, 1, 1, 0, 0}),
                                  std::make_pair("sigma_s", std::vector<double>{0, 0, 2, 2, 10, 10}),
                                  std::make_pair("Q", std::vector<double>{0, 0, 1, 0, 0, 0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{4.0})});
  }

  // Boundary value of kinetic equation is \frac{g}{<g>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, with g(v) = exp(-10^5(v-1)^2), so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // Simulate with linear interpolating function left_boundary_vals * (1-x/3) + base_integrated*psi_vac * x/3
  virtual BoundaryValueType* create_boundary_values() const override
  {
    const RangeType basis_integrated = basis_functions_.integrated();
    const RangeType left_boundary_values =
        helper<BasisfunctionType>::get_left_boundary_values(quadrature_, basis_functions_);
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          RangeType ret = left_boundary_values;
          ret *= 1 - x[0] / 3.;
          RangeType summand2 = basis_integrated;
          summand2 *= psi_vac_ * x[0] / 3.;
          ret += summand2;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  struct helper_base
  {
    // returns the numerator g of the left boundary value (see create_boundary_values)
    static RangeFieldType numerator(const RangeFieldType v)
    {
      return std::exp(-1e5 * (v - 1) * (v - 1));
    }

    // returns the denominator <g> of the left boundary value (see create_boundary_values)
    static RangeFieldType denominator()
    {
      return 1 / 200. * std::sqrt(M_PI / 10) * std::erf(200 * std::sqrt(10));
    }

    // calculates integral from v_l to v_u of numerator g
    static RangeFieldType integral_1(RangeFieldType v_l, RangeFieldType v_u)
    {
      return 1 / 200. * std::sqrt(M_PI / 10)
             * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
    }

    // calculates integral from v_l to v_u of v*g
    static RangeFieldType integral_2(RangeFieldType v_l, RangeFieldType v_u)
    {
      return integral_1(v_l, v_u) - 1. / 2e5 * (numerator(v_u) - numerator(v_l));
    }
  };

  template <class B, class anything = void>
  struct helper : public helper_base
  {
    using helper_base::numerator;
    using helper_base::denominator;

    static RangeType get_left_boundary_values(const QuadratureType& quadrature, const BasisfunctionImp& basis_functions)
    {
      RangeType ret(0);
      for (const auto& quad_point : quadrature) {
        const auto& v = quad_point.position()[0];
        auto summand = basis_functions.evaluate(v);
        summand *= numerator(v) * quad_point.weight();
        ret += summand;
      }
      ret /= denominator();
      return ret;
    }
  };

  template <class anything>
  struct helper<HatFunctions<DomainFieldType, dimDomain, RangeFieldType, dimRange>, anything> : public helper_base
  {
    using helper_base::numerator;
    using helper_base::denominator;
    using helper_base::integral_1;

    static RangeType get_left_boundary_values(const QuadratureType& /*quadrature*/,
                                              const BasisfunctionImp& basis_functions)
    {
      RangeType ret(0);
      for (size_t nn = 0; nn < dimRange; ++nn) {
        const auto& triangulation = basis_functions.triangulation();
        const auto vnm = triangulation[nn - 1];
        const auto vn = triangulation[nn];
        const auto vnp = triangulation[nn + 1];
        if (nn < dimRange - 1)
          ret[nn] += 1. / ((vn - vnp) * denominator())
                     * ((1 - vnp) * integral_1(vn, vnp) - 1. / 2e5 * (numerator(vnp) - numerator(vn)));
        if (nn > 0)
          ret[nn] += 1. / ((vn - vnm) * denominator())
                     * ((1 - vnm) * integral_1(vnm, vn) - 1. / 2e5 * (numerator(vn) - numerator(vnm)));
      }
      return ret;
    }
  };

  template <class anything>
  struct helper<PiecewiseMonomials<DomainFieldType, dimDomain, RangeFieldType, dimRange>, anything> : public helper_base
  {
    using helper_base::denominator;
    using helper_base::integral_1;
    using helper_base::integral_2;

    static RangeType get_left_boundary_values(const QuadratureType& /*quadrature*/,
                                              const BasisfunctionImp& basis_functions)
    {
      const auto& triangulation = basis_functions.triangulation();
      RangeType ret(0);
      for (size_t ii = 0; ii < dimRange / 2; ++ii) {
        ret[2 * ii] = integral_1(triangulation[ii], triangulation[ii + 1]) / denominator();
        ret[2 * ii + 1] = integral_2(triangulation[ii], triangulation[ii + 1]) / denominator();
      }
      return ret;
    }
  };

  using BaseType::basis_functions_;
  using BaseType::quadrature_;
  using BaseType::psi_vac_;
}; // class SourceBeamPn<...>

template <class BasisfunctionType,
          class GridLayerType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class SourceBeamMn : public SourceBeamPn<BasisfunctionType,
                                         GridLayerType,
                                         EntityType,
                                         DomainFieldType,
                                         dimDomain,
                                         U_,
                                         RangeFieldType,
                                         dimRange>
{
  typedef SourceBeamPn<BasisfunctionType,
                       GridLayerType,
                       EntityType,
                       DomainFieldType,
                       dimDomain,
                       U_,
                       RangeFieldType,
                       dimRange>
      BaseType;
  typedef SourceBeamMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef EntropyBasedLocalFlux<BasisfunctionType,
                                GridLayerType,
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
  using BaseType::default_quadrature;

  SourceBeamMn(const BasisfunctionType& basis_functions,
               const GridLayerType& grid_layer,
               const QuadratureType& quadrature = default_quadrature();
               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "sourcebeammn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_layer_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::quadrature_;
}; // class SourceBeamMn<...>


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
