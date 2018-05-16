// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH

#include <dune/xt/common/string.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class PlaneSourcePn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_, 1>
{
  typedef KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_, 1> BaseType;

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

  PlaneSourcePn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const QuadratureType& quadrature = default_quadrature(),
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, {1}, grid_cfg, boundary_cfg)
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
    grid_config["num_quad_cells"] = "[25]";
    grid_config["quad_order"] = "30";
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
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    size_t num_elements = grid_layer_.size(0);
    if (grid_layer_.comm().size() > 1) {
      num_elements -= grid_layer_.overlapSize(0);
      num_elements = grid_layer_.comm().sum(num_elements);
    }
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    const RangeFieldType len_domain = upper_right[0] - lower_left[0];
    const RangeFieldType vol_entity = len_domain / num_elements;
    RangeType basis_integrated = basis_functions_.integrated();
    const RangeFieldType domain_center = lower_left[0] + len_domain / 2;
    // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    initial_vals.emplace_back(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          if (XT::Common::FloatCmp::ge(x[0], domain_center - vol_entity)
              && XT::Common::FloatCmp::le(x[0], domain_center + vol_entity))
            ret *= psi_vac_ + 1. / (2. * vol_entity);
          else
            ret *= psi_vac_;
          return ret;
        },
        0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

protected:
  using BaseType::grid_cfg_;
  using BaseType::grid_layer_;
  using BaseType::basis_functions_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
}; // class PlaneSourcePn<...>


template <class BasisfunctionType, class GridLayerImp, class U_>
class PlaneSourceMn : public PlaneSourcePn<BasisfunctionType, GridLayerImp, U_>
{
  typedef PlaneSourcePn<BasisfunctionType, GridLayerImp, U_> BaseType;
  typedef PlaneSourceMn ThisType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  typedef EntropyBasedLocalFlux<BasisfunctionType, GridLayerImp, U_> ActualFluxType;
  using typename BaseType::QuadratureType;
  using typename BaseType::GridLayerType;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;
  using BaseType::default_quadrature;

  PlaneSourceMn(const BasisfunctionType& basis_functions,
                const GridLayerType& grid_layer,
                const QuadratureType& quadrature = default_quadrature(),
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, quadrature, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "planesourcemn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_, quadrature_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
  using BaseType::quadrature_;
}; // class PlaneSourceMn<...>


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
