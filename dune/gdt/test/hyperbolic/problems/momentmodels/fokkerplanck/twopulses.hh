// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_TWOPULSES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_TWOPULSES_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "fokkerplanckequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace FokkerPlanck {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class TwoPulsesPn : public FokkerPlanckEquation<BasisfunctionImp, GridLayerImp, U_>
{
  typedef FokkerPlanckEquation<BasisfunctionImp, GridLayerImp, U_> BaseType;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::QuadratureType;

  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  static std::string static_id()
  {
    return "twopulsespn";
  }

  TwoPulsesPn(const BasisfunctionType& basis_functions,
              const GridLayerType& grid_layer,
              const QuadratureType& quadrature = default_quadrature(),
              const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
              const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, 1, grid_cfg, boundary_cfg)
  {
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
                                  std::make_pair("t_end", std::vector<double>{7.0})});
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


} // namespace FokkerPlanck
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_FOKKERPLANCK_TWOPULSES_HH
