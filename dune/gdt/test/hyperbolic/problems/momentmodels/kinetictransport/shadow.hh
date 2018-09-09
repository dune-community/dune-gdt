// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_Shadow_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_Shadow_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/grid/boundaryinfo.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class ShadowPn : public KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>
{
  using BaseType = KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>;

public:
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DirichletBoundaryValueType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;
  using ActualBoundaryValueType =
      MomentModelBoundaryValue<IntersectionType, BasisfunctionType, DirichletBoundaryValueType>;
  static const size_t dimDomain = BaseType::dimDomain;

  using BaseType::default_boundary_cfg;

  ShadowPn(const BasisfunctionType& basis_functions,
           const GridLayerType& grid_layer,
           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, {12, 4, 3}, grid_cfg, boundary_cfg, 1e-8 / (4 * M_PI))
  {
  }

  static std::string static_id()
  {
    return "Shadowpn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0 0 0]";
    grid_config["upper_right"] = "[12 4 3]";
    grid_config["num_elements"] = "[4 4 4]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // sigma_a = 50 on [2,3]x[1,3]x[0,2], sigma_s = 0, Q = 0
  virtual XT::Common::Parameter parameters() const override
  {
    return XT::Common::Parameter({std::make_pair("sigma_a", create_sigma_a()),
                                  std::make_pair("sigma_s", std::vector<double>(12 * 4 * 3, 0.)),
                                  std::make_pair("Q", std::vector<double>(12 * 4 * 3, 0.)),
                                  std::make_pair("CFL", std::vector<double>{0.49 * 1 / std::sqrt(dimDomain)}),
                                  std::make_pair("t_end", std::vector<double>{4})});
  }

  // Boundary value of kinetic equation is \dirac(v - (1, 0, 0)) at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 12, reflecting else
  // so n-th component of \psi_{vac}*base_integrated_n
  // at x = 3.
  virtual BoundaryValueType* create_boundary_values() const override
  {
    auto dirac_integrated = basis_functions_.integrate_dirac_at(DomainType{1, 0, 0});
    auto basis_integrated = basis_functions_.integrated();
    auto dirichlet_boundary_values = std::make_unique<ActualDirichletBoundaryValueType>(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          ret *= psi_vac_;
          if (XT::Common::FloatCmp::eq(x[0], 0.))
            ret += dirac_integrated;
          return ret;
        },
        1);
    auto boundary_info = std::make_unique<XT::Grid::NormalBasedBoundaryInfo<IntersectionType>>(
        1e-10, new XT::Grid::ReflectingBoundary{});
    boundary_info->register_new_normal({-1, 0, 0}, new XT::Grid::DirichletBoundary{});
    boundary_info->register_new_normal({1, 0, 0}, new XT::Grid::DirichletBoundary{});
    return new ActualBoundaryValueType(
        std::move(boundary_info), basis_functions_, std::move(dirichlet_boundary_values));
  } // ... create_boundary_values()

protected:
  static std::vector<double> create_sigma_a()
  {
    std::vector<double> sigma_a(12 * 4 * 3, 0.);
    size_t x = 2;
    for (size_t y = 1; y < 3; ++y)
      for (size_t z = 0; z < 2; ++z)
        sigma_a[x + 12 * y + 12 * 4 * z] = 50.;
    return sigma_a;
  }

  using BaseType::grid_cfg_;
  using BaseType::basis_functions_;
  using BaseType::num_segments_;
  using BaseType::psi_vac_;
  using BaseType::grid_layer_;
}; // class ShadowPn<...>

template <class BasisfunctionType, class GridLayerType, class U_>
class ShadowMn : public ShadowPn<BasisfunctionType, GridLayerType, U_>
{
  using BaseType = ShadowPn<BasisfunctionType, GridLayerType, U_>;
  using ThisType = ShadowMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeType;
  using ActualFluxType = GDT::EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, U_>;

  using BaseType::default_grid_cfg;
  using BaseType::default_boundary_cfg;

  ShadowMn(const BasisfunctionType& basis_functions,
           const GridLayerType& grid_layer,
           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_layer, grid_cfg, boundary_cfg)
  {
  }

  static std::string static_id()
  {
    return "Shadowmn";
  }

  virtual FluxType* create_flux() const
  {
    return new ActualFluxType(basis_functions_, grid_layer_);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_layer_;
}; // class ShadowMn<...>


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_Shadow_HH
