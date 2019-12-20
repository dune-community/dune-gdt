// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH

#include "base.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class PointSourcePn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;

public:
  using BaseType::dimDomain;
  using typename BaseType::ConstantScalarFunctionType;
  using typename BaseType::DomainType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::InitialValueType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::ScalarFunctionType;

  using BaseType::default_boundary_cfg;

  PointSourcePn(const MomentBasis& basis_functions,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg, 1e-8 / (4 * M_PI))
  {}

  static std::string static_id()
  {
    return "pointsourcepn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[-1 -1 -1]";
    grid_config["upper_right"] = "[1 1 1]";
    grid_config["num_elements"] = "[10 10 10]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  // Initial value of the kinetic equation is psi_vac + 1/(4 pi^4 sigma^3) * exp(-||x||^2/(pi*sigma^2)).
  std::unique_ptr<InitialValueType> initial_values() const override
  {
    const RangeReturnType basis_integrated = basis_functions_.integrated();
    const auto psi_vac = psi_vac_;
    const auto eval_func = [basis_integrated, psi_vac](const DomainType& x, const XT::Common::Parameter&) {
      static const auto sigma = 0.03;
      static const auto first_factor = 1. / (4 * M_PI * std::pow(M_PI * sigma, 3));
      static const auto second_factor = 1. / (M_PI * std::pow(sigma, 2));
      return basis_integrated * std::max(first_factor * std::exp(-x.two_norm2() * second_factor), psi_vac);
    };
    return std::make_unique<GenericFunctionType>(21, eval_func);
  } // ... initial_values()

  RangeFieldType t_end() const override
  {
    return 0.75;
  }

  // sigma_a = 0, sigma_s = 1, Q = 0
  std::unique_ptr<ScalarFunctionType> sigma_a() const override
  {
    return std::make_unique<ConstantScalarFunctionType>(0.);
  }

  std::unique_ptr<ScalarFunctionType> sigma_s() const override
  {
    return std::make_unique<ConstantScalarFunctionType>(1.);
  }

  std::unique_ptr<ScalarFunctionType> Q() const override
  {
    return std::make_unique<ConstantScalarFunctionType>(0.);
  }

protected:
  using BaseType::basis_functions_;
  using BaseType::grid_cfg_;
  using BaseType::psi_vac_;
}; // class PointSourcePn<...>

template <class GV, class MomentBasis>
class PointSourceMn : public PointSourcePn<XT::Grid::extract_entity_t<GV>, MomentBasis>
{
  using BaseType = PointSourcePn<XT::Grid::extract_entity_t<GV>, MomentBasis>;
  using ThisType = PointSourceMn;

public:
  using typename BaseType::FluxType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  PointSourceMn(const MomentBasis& basis_functions,
                const GV& grid_view,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {}

  static std::string static_id()
  {
    return "pointsourcemn";
  }

  std::unique_ptr<FluxType> flux() const override final
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
}; // class PointSourceMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
