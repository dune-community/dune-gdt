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
class RectangularICPn : public KineticFokkerPlanckEquation<BasisfunctionImp,
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
    return "rectangularicpn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0]";
    grid_config["upper_right"] = "[7]";
    grid_config["num_elements"] = "[70]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  // sigma_a = 0, T = 0.01 and Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("T", std::vector<double>{0.01}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{1}),
                                  std::make_pair("num_segments", std::vector<double>{1.})});
  }

  // Initial value of the kinetic equation is 10 if 3 <= x <= 4 and psi_vac else, thus initial value of the
  // k-th component of the moment vector is 10*base_integrated_k or psi_vac*base_integrated_k.
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const FieldVector<size_t, 3> num_segments(1);
    RangeType basis_integrated = basis_functions_.integrated();
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x) {
          auto ret = basis_integrated;
          ret *= XT::Common::FloatCmp::ge(x, 3.) && XT::Common::FloatCmp::le(x, 4.) ? 10 : psi_vac_;
          return ret;
        },
        1);
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
                  const GridViewType& grid_view,
                  const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                  const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
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
  using BaseType::quadrature_;
}; // class RectangularICMn<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_RECTANGULARIC_HH
