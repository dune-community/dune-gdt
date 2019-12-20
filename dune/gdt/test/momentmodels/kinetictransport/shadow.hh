// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH

#include "base.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class ShadowPn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;

public:
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConstantScalarFunctionType;
  using typename BaseType::DomainType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::GenericScalarFunctionType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::ScalarFunctionType;

  using BaseType::default_boundary_cfg;

  ShadowPn(const MomentBasis& basis_functions,
           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg, 1e-8 / (4 * M_PI))
  {}

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

  RangeFieldType t_end() const override
  {
    return 20;
  }

  // sigma_a = 50 on [2,3]x[1,3]x[0,2], sigma_s = 0, Q = 0
  std::unique_ptr<ScalarFunctionType> sigma_a() const override
  {
    return std::make_unique<GenericScalarFunctionType>(0, [=](const DomainType& x, const XT::Common::Parameter&) {
      return (x[0] < 2 || x[0] > 3 || x[1] < 1 || x[1] > 3 || x[2] < 0 || x[2] > 2) ? 0. : 50.;
    });
  }

  std::unique_ptr<ScalarFunctionType> sigma_s() const override
  {
    return std::make_unique<ConstantScalarFunctionType>(0.);
  }

  std::unique_ptr<ScalarFunctionType> Q() const override
  {
    return std::make_unique<ConstantScalarFunctionType>(0.);
  }

#define USE_DIRAC_BOUNDARY 0
  // Boundary value of kinetic equation is either \dirac(v - (1, 0, 0)) or an isotropic beam with density 2,
  // i.e. 2/(4 pi), at x = 0 and psi_vac else
  std::unique_ptr<BoundaryValueType> boundary_values() const override final
  {
    const auto basis_integrated = basis_functions_.integrated();
    const auto psi_vac = psi_vac_;
    return std::make_unique<GenericFunctionType>(
        0, [basis_integrated, psi_vac](const DomainType& x, const XT::Common::Parameter&) {
          auto ret = basis_integrated;
          ret *= psi_vac;
          // left boundary
          if (XT::Common::FloatCmp::eq(x[0], 0.)) {
#if USE_DIRAC_BOUNDARY
            auto dirac_integrated = basis_functions_.integrate_dirac_at(DomainType{1, 0, 0});
            ret += dirac_integrated;
#else // USE_DIRAC_BOUNDARY
                                                    ret += basis_integrated * 2. / (4 * M_PI);
#endif // USE_DIRAC_BOUNDARY
          }
          return ret;
        });
  } // ... boundary_values()

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_cfg_;
  using BaseType::psi_vac_;
}; // class ShadowPn<...>

template <class GV, class MomentBasis>
class ShadowMn : public ShadowPn<XT::Grid::extract_entity_t<GV>, MomentBasis>
{
  using BaseType = ShadowPn<XT::Grid::extract_entity_t<GV>, MomentBasis>;
  using ThisType = ShadowMn;

public:
  using typename BaseType::FluxType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  ShadowMn(const MomentBasis& basis_functions,
           const GV& grid_view,
           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {}

  static std::string static_id()
  {
    return "Shadowmn";
  }

  std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
}; // class ShadowMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH
