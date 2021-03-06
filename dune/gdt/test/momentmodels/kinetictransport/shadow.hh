// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH

#if HAVE_DUNE_XT_DATA

#  include "base.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class ShadowPn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;

public:
  using typename BaseType::BoundaryDistributionType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConstantScalarFunctionType;
  using typename BaseType::DomainType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::GenericScalarFunctionType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::ScalarFunctionType;

  ShadowPn(const MomentBasis& basis_functions,
           const RangeFieldType psi_vac = 1e-6 / (4 * M_PI),
           const XT::Common::Configuration& grid_cfg = default_grid_cfg())
    : BaseType(basis_functions, psi_vac, grid_cfg)
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

#  define USE_DIRAC_BOUNDARY 0
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
#  if USE_DIRAC_BOUNDARY
            auto dirac_integrated = basis_functions_.integrate_dirac_at(DomainType{1, 0, 0});
            ret += dirac_integrated;
#  else // USE_DIRAC_BOUNDARY
                                                    ret += basis_integrated * 2. / (4 * M_PI);
#  endif // USE_DIRAC_BOUNDARY
          }
          return ret;
        });
  } // ... boundary_values()

  BoundaryDistributionType boundary_distribution() const override final
  {
    return [this](const DomainType& x) -> std::function<RangeFieldType(const DomainType&)> {
      if (XT::Common::FloatCmp::eq(x[0], 0.))
        return [](const DomainType& /*v*/) { return 2. / (4 * M_PI); };
      else
        return [this](const DomainType& /*v*/) { return psi_vac_; };
    };
  }

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
  using typename BaseType::RangeFieldType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_grid_cfg;

  ShadowMn(const MomentBasis& basis_functions,
           const GV& grid_view,
           const RangeFieldType psi_vac = 1e-6 / (4 * M_PI),
           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
           const bool disable_realizability_check = false,
           const RangeFieldType tau = 1e-9)
    : BaseType(basis_functions, psi_vac, grid_cfg)
    , grid_view_(grid_view)
    , disable_realizability_check_(disable_realizability_check)
    , tau_(tau)
  {}

  static std::string static_id()
  {
    return "Shadowmn";
  }

  std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_, tau_, disable_realizability_check_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
  const bool disable_realizability_check_;
  const RangeFieldType tau_;
}; // class ShadowMn<...>


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SHADOW_HH
