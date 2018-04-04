// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_PLANESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_PLANESOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace FokkerPlanck {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class PlaneSourcePn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_> BaseType;

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
  using typename BaseType::GridLayerType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  PlaneSourcePn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const QuadratureType& quadrature = default_quadrature(),
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, 1, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "planesourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-1.2]";
    grid_config["upper_right"] = "[1.2]";
    grid_config["num_elements"] = "[240]";
    grid_config["overlap_size"] = "[1]";
    grid_config["num_quad_cells"] = "[20]";
    grid_config["quad_order"] = "20";
    return grid_config;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", std::vector<double>{0}),
                                  std::make_pair("sigma_s", std::vector<double>{1}),
                                  std::make_pair("Q", std::vector<double>{0}),
                                  std::make_pair("CFL", std::vector<double>{0.4}),
                                  std::make_pair("t_end", std::vector<double>{1.0})});
  }

  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  virtual InitialValueType* create_initial_values() const
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const size_t num_elements = XT::Common::from_string<std::vector<size_t>>(grid_cfg_["num_elements"])[0];
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    const RangeFieldType len_domain = upper_right[0] - lower_left[0];
    const RangeFieldType vol_entity = len_domain / num_elements;
    RangeType basis_integrated = basis_functions_.integrated();
    const size_t num_segments = 1;
    const RangeFieldType domain_center = lower_left[0] + len_domain / 2;
    // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x) {
          auto ret = basis_integrated;
          if (XT::Common::FloatCmp::ge(x[0], domain_center - vol_entity)
              && XT::Common::FloatCmp::le(x[0], domain_center + vol_entity))
            ret *= psi_vac_ + 1. / (2. * vol_entity);
          else
            ret *= psi_vac_;
          return ret;
        },
        0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::psi_vac_;
}; // class PlaneSourcePn<...>


} // namespace FokkerPlanck
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_PLANESOURCE_HH
