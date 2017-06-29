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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "kineticequation.hh"

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
class OneBeamPn : public KineticFokkerPlanckEquation<BasisfunctionImp,
                                                     EntityImp,
                                                     DomainFieldImp,
                                                     dimDomain,
                                                     U_,
                                                     RangeFieldImp,
                                                     dimRange>
{
  typedef KineticFokkerPlanckEquation<BasisfunctionImp,
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
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;

  using BaseType::default_boundary_cfg;

  OneBeamPn(const BasisfunctionType& basis_functions,
            const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
            const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "onebeampn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-0.5]";
    grid_config["upper_right"] = "[0.5]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    grid_config["num_quad_cells"] = "[20]";
    grid_config["quad_order"] = "50";
    return grid_config;
  }

  // sigma_a = 10 if -0.1 <= x <= 0.2, 0 else
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0, 0, 0, 0, 10, 10, 10, 0, 0, 0}),
                                  std::make_pair("T", std::vector<double>(10, 0)),
                                  std::make_pair("Q", std::vector<double>(10, 0)),
                                  std::make_pair("CFL", std::vector<double>{0.5}),
                                  std::make_pair("t_end", std::vector<double>{4.0}),
                                  std::make_pair("num_segments", std::vector<double>{10})});
  }


  // Boundary value of kinetic equation is 3*exp(3*v + 3)/(exp(6) - 1) at x = -0.5 and psi_vac at x = 0.5,
  // so k-th component of boundary value has to be \int 3*exp(3*v + 3)/(exp(6) - 1) \phi_k(v) dv at x = -0.5 and
  // psi_vac*base_integrated_k at x = 0.5.
  // Model with linear interpolating function.
  virtual BoundaryValueType* create_boundary_values() const override
  {
    RangeType basis_integrated = basis_functions_.integrated();
    RangeType left_boundary_value(0);
    for (const auto& quadpoint : quadrature_) {
      const auto& v = quad_point.position();
      left_boundary_value +=
          basis_functions_.evaluate(v) * 3. * std::exp(3. * v[0] + 3.) / (std::exp(6) - 1) * quadpoint.weight();
    }
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = left_boundary_value;
          ret *= 0.5 - x[0];
          auto right_boundary_value = basis_integrated;
          right_boundary_value *= psi_vac * (0.5 + x[0]);
          ret += right_boundary_value;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  using BaseType::basis_functions_;
  using BaseType::quadrature_;
}; // class OneBeamPn<...>

template <class GridViewType,
          class BasisfunctionType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class OneBeamMn
    : public OneBeamPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
{
  typedef OneBeamPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange> BaseType;
  typedef OneBeamMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using typename BaseType::QuadratureType;
  typedef EntropyBasedLocalFlux<GridViewType,
                                BasisfunctionType,
                                EntityType,
                                DomainFieldType,
                                dimDomain,
                                U_,
                                RangeFieldType,
                                dimRange>
      ActualFluxType;


  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  OneBeamMn(const BasisfunctionType& basis_functions,
            const GridViewType& grid_view,
            const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
            const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {
  }

  static std::string static_id()
  {
    return "onebeammn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_view_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GridViewType& grid_view_;
  using BaseType::quadrature_;
}; // class OneBeamMn<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_ONEBEAM_HH
