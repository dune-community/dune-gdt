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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH

#include "base.hh"

namespace Dune {
namespace GDT {


/**
 * In 2D, this is the Testcase for the Boltzmann equation in two dimensions,
 * see Section 4.2 of Brunner, Holloway, "Two-dimensional time dependent Riemann solvers for neutron transport", Journal
 * of Computational Physics, Volume 210, Issue 1, 2005
 * http://dx.doi.org/10.1016/j.jcp.2005.04.011
 * The 3D version is a straightforward generalization of the setup to 3 dimensions.
 * */
template <class E, class MomentBasisImp>
class CheckerboardPn : public KineticTransportEquationBase<E, MomentBasisImp>
{
  using BaseType = KineticTransportEquationBase<E, MomentBasisImp>;

public:
  using BaseType::dimDomain;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::GenericScalarFunctionType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::ScalarFunctionType;

  using BaseType::default_boundary_cfg;

  CheckerboardPn(const MomentBasis& basis_functions,
                 const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg, 1e-8 / (4 * M_PI))
  {}

  static std::string static_id()
  {
    size_t domainDim = dimDomain; // avoid linker error for dimDomain
    return "checkerboard" + XT::Common::to_string(domainDim) + "d_pn";
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[7.0 7.0 7.0]";
    grid_config["num_elements"] = "[14 14 14]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  virtual RangeFieldType t_end() const override
  {
    return 3.2;
  }

  // Q = 0 except in the center where Q = 1. sigma_s = sigma_t = 1 in scattering regions, sigma_s = 0, sigma_t
  // = 10 in absorbing regions. Center is also a scattering region.
  virtual std::unique_ptr<ScalarFunctionType> sigma_a() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        [](const XT::Common::Parameter&) { return 0; },
        [=](const DomainType& x, const XT::Common::Parameter&) { return is_absorbing(x) ? 10. : 0.; });
  }

  virtual std::unique_ptr<ScalarFunctionType> sigma_s() const override
  {
    return std::make_unique<GenericScalarFunctionType>(
        [](const XT::Common::Parameter&) { return 0; },
        [=](const DomainType& x, const XT::Common::Parameter&) { return is_absorbing(x) ? 0. : 1.; });
  }

  virtual std::unique_ptr<ScalarFunctionType> Q() const override
  {
    return std::make_unique<GenericScalarFunctionType>([](const XT::Common::Parameter&) { return 0; },
                                                       [=](const DomainType& x, const XT::Common::Parameter&) {
                                                         return is_center(x) ? 1. / std::sqrt(4. * M_PI) : 0.;
                                                       });
  }

protected:
  static bool is_absorbing(const FieldVector<RangeFieldType, 2>& x)
  {
    size_t row = XT::Common::numeric_cast<size_t>(std::floor(x[1]));
    if (row == 7)
      row = 6;
    size_t col = XT::Common::numeric_cast<size_t>(std::floor(x[0]));
    if (col == 7)
      col = 6;
    assert(row < 7 && col < 7);
    return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
           || ((row == 3 || row == 5) && (col == 1 || col == 5));
  }

  static bool is_absorbing(const FieldVector<RangeFieldType, 3>& x)
  {
    size_t plane = XT::Common::numeric_cast<size_t>(std::floor(x[2]));
    size_t col = XT::Common::numeric_cast<size_t>(std::floor(x[0]));
    size_t row = XT::Common::numeric_cast<size_t>(std::floor(x[1]));
    assert(plane <= 7 && row <= 7 && col <= 7);
    if (plane == 0 && plane >= 6)
      return false;
    if (row == 0 && row >= 6)
      return false;
    if (col == 0 && col >= 6)
      return false;
    return (plane + row + col) % 2 == 1 && !is_center(x) && !(plane == 3 && row == 5 && col == 3);
  }

  static bool is_center(const FieldVector<RangeFieldType, 2>& x)
  {
    size_t row = XT::Common::numeric_cast<size_t>(std::floor(x[1]));
    size_t col = XT::Common::numeric_cast<size_t>(std::floor(x[0]));
    return row == 3 && col == 3;
  }

  static bool is_center(const FieldVector<RangeFieldType, 3>& x)
  {
    size_t plane = XT::Common::numeric_cast<size_t>(std::floor(x[2]));
    size_t row = XT::Common::numeric_cast<size_t>(std::floor(x[1]));
    size_t col = XT::Common::numeric_cast<size_t>(std::floor(x[0]));
    return plane == 3 && row == 3 && col == 3;
  }
}; // class CheckerboardPn<...>

template <class GV, class MomentBasis>
class CheckerboardMn : public CheckerboardPn<XT::Grid::extract_entity_t<GV>, MomentBasis>
{
  using BaseType = CheckerboardPn<XT::Grid::extract_entity_t<GV>, MomentBasis>;

public:
  using typename BaseType::FluxType;
  using ActualFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  CheckerboardMn(const MomentBasis& basis_functions,
                 const GV& grid_view,
                 const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(basis_functions, grid_cfg, boundary_cfg)
    , grid_view_(grid_view)
  {}

  static std::string static_id()
  {
    return "checkerboardmn";
  }

  virtual std::unique_ptr<FluxType> flux() const override
  {
    return std::make_unique<ActualFluxType>(grid_view_, basis_functions_);
  }

protected:
  using BaseType::basis_functions_;
  const GV& grid_view_;
}; // class CheckerboardMn<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORT_CHECKERBOARD_HH
