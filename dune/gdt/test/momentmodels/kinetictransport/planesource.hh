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

#include "base.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class PlaneSourcePn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;

public:
  using BaseType::default_boundary_cfg;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::ConstantScalarFunctionType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::InitialValueType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::ScalarFunctionType;

  PlaneSourcePn(const MomentBasis& basis_functions,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
  {}

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
    return grid_config;
  }

  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  std::unique_ptr<InitialValueType> initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    const size_t num_elements = XT::Common::from_string<std::vector<size_t>>(grid_cfg_["num_elements"])[0];
    const RangeFieldType len_domain = upper_right[0] - lower_left[0];
    const RangeFieldType vol_entity = len_domain / num_elements;
    const auto basis_integrated = basis_functions_.integrated();
    const RangeFieldType domain_center = lower_left[0] + len_domain / 2;

    // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
    const auto eval_func = [=](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) {
      ret = basis_integrated;
      if (XT::Common::FloatCmp::ge(x[0], domain_center - vol_entity)
          && XT::Common::FloatCmp::le(x[0], domain_center + vol_entity))
        ret *= psi_vac_ + 1. / (2. * vol_entity);
      else
        ret *= psi_vac_;
    };
    return std::make_unique<GenericFunctionType>(0, eval_func);
  } // ... initial_values()

  RangeFieldType t_end() const override
  {
    return 1.0;
  }

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
}; // class PlaneSourcePn<...>


template <class GV, class MomentBasis>
class PlaneSourceMn : public PlaneSourcePn<XT::Grid::extract_entity_t<GV>, MomentBasis>
{
  using BaseType = PlaneSourcePn<XT::Grid::extract_entity_t<GV>, MomentBasis>;
  using ThisType = PlaneSourceMn;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RangeReturnType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  PlaneSourceMn(const MomentBasis& basis_functions,
                const GV& grid_view,
                const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {}

  static std::string static_id()
  {
    return "planesourcemn";
  }

  std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
}; // class PlaneSourceMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
