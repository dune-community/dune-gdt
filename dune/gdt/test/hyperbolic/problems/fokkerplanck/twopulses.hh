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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH

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
class TwoPulsesPn : public KineticFokkerPlanckEquation<BasisfunctionImp,
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

  TwoPulsesPn(const BasisfunctionType& basis_functions,
              const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
              const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "twopulsespn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0]";
    grid_config["upper_right"] = "[7]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    grid_config["num_quad_cells"] = "[20]";
    grid_config["quad_order"] = "50";
    return grid_config;
  }

  // sigma_a = T = Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("T", std::vector<double>{0}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{7.0}),
                                  std::make_pair("num_segments", std::vector<double>{1})});
  }


  // Boundary value of kinetic equation is 100*delta(v-1)**exp(-(t-1)^2/2) at x = 0 and 100*delta(v+1)*exp(-(t-1)^2/2)
  // at x = 7, so k-th component of boundary value has to be 50*\phi_k(1)*exp(-(t-1)^2/2) at x = 0 and
  // 50*\phi_k(-1)*exp(-(t-1)^2/2) at x = 7.
  // Model with linear interpolating function.
  virtual BoundaryValueType* create_boundary_values() const override
  {
    const auto basis_evaluated_at_one = basis_functions_.evaluate(DomainType(1));
    const auto basis_evaluated_at_minus_one = basis_functions_.evaluate(DomainType(-1));
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter& param) {
          const RangeFieldType t = param.get("t")[0];
          auto ret = basis_evaluated_at_one;
          ret *= 50 * std::exp(-(t - 1) * (t - 1) / 2.) * (1. - x[0] / 7.);
          auto right_boundary_value = basis_evaluated_at_minus_one;
          right_boundary_value *= 50 * std::exp(-(t - 1) * (t - 1) / 2.) * (x[0] / 7.);
          ret += right_boundary_value;
          return ret;
        },
        1);
  } // ... create_boundary_values()

protected:
  using BaseType::basis_functions_;
  using BaseType::quadrature_;
}; // class TwoPulsesPn<...>

template <class GridViewType,
          class BasisfunctionType,
          class EntityType,
          class DomainFieldType,
          size_t dimDomain,
          class U_,
          class RangeFieldType,
          size_t dimRange>
class TwoPulsesMn
    : public TwoPulsesPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange>
{
  typedef TwoPulsesPn<BasisfunctionType, EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange> BaseType;
  typedef TwoPulsesMn ThisType;

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

  TwoPulsesMn(const BasisfunctionType& basis_functions,
              const GridViewType& grid_view,
              const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
              const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {
  }

  static std::string static_id()
  {
    return "twopulsesmn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(grid_view_, quadrature_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GridViewType& grid_view_;
  using BaseType::quadrature_;
}; // class TwoPulsesMn<...>

} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOPULSES_HH
