// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/functions/flattop.hh>

#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>

#include "inviscid_compressible_flow__euler_1d__explicit.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class GL, class R, size_t m>
using DgOrder1SpaceType = DuneFemDgSpaceWrapper<GL, 1, R, m, 1>;

template <class DF>
using DgOperator = GDT::AdvectionDgOperator<DF, typename DF::SpaceType::GridLayerType>;

using InviscidCompressibleFlowEuler1dExplicitDgOrder1 =
    InviscidCompressibleFlowEuler1dExplicitTest<DgOrder1SpaceType, DgOperator, XT::Grid::Backends::part>;


TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder1, periodic_boundaries)
{
  // use smooth initial values to avoid the Gibbs phenomenon
  ASSERT_NE(initial_values_, nullptr);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> density_minus_one({-0.5}, {0.}, 0.015625, 3.);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> pressure_minus_04({-0.5}, {0.}, 0.015625, 1.2);
  U initial_values(
      [&](const auto& xx, const auto& /*param*/) {
        return euler_tools_.to_conservative(
            XT::Common::hstack(density_minus_one.evaluate(xx) + 1., 0., pressure_minus_04.evaluate(xx) + 0.4));
      },
      1);
  project(initial_values, *initial_values_);

  const double CFL = 0.15;
  this->periodic_boundaries(100 / CFL, CFL, {1e-15, 0.07254});
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder1, impermeable_walls_by_direct_euler_treatment)
{
  // use smooth initial values to avoid the Gibbs phenomenon
  ASSERT_NE(initial_values_, nullptr);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> density_minus_one({-0.5}, {0.}, 0.015625, 3.);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> pressure_minus_04({-0.5}, {0.}, 0.015625, 1.2);
  U initial_values(
      [&](const auto& xx, const auto& /*param*/) {
        return euler_tools_.to_conservative(
            XT::Common::hstack(density_minus_one.evaluate(xx) + 1., 0., pressure_minus_04.evaluate(xx) + 0.4));
      },
      1);
  project(initial_values, *initial_values_);

  const double CFL = 0.15;
  this->impermeable_walls_by_direct_euler_treatment(300 / CFL, CFL, 0.08274);
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder1, impermeable_walls_by_inviscid_mirror_treatment)
{
  // use smooth initial values to avoid the Gibbs phenomenon
  ASSERT_NE(initial_values_, nullptr);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> density_minus_one({-0.5}, {0.}, 0.015625, 3.);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> pressure_minus_04({-0.5}, {0.}, 0.015625, 1.2);
  U initial_values(
      [&](const auto& xx, const auto& /*param*/) {
        return euler_tools_.to_conservative(
            XT::Common::hstack(density_minus_one.evaluate(xx) + 1., 0., pressure_minus_04.evaluate(xx) + 0.4));
      },
      1);
  project(initial_values, *initial_values_);

  const double CFL = 0.15;
  this->impermeable_walls_by_inviscid_mirror_treatment(300 / CFL, CFL, 0.08313);
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder1,
       inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right)
{
  // use smooth initial values to avoid the Gibbs phenomenon
  ASSERT_NE(initial_values_, nullptr);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> density_minus_one({-0.5}, {0.}, 0.015625, 3.);
  XT::Functions::FlatTopFunction<E, D, d, R, 1> pressure_minus_04({-0.5}, {0.}, 0.015625, 1.2);
  U initial_values(
      [&](const auto& xx, const auto& /*param*/) {
        return euler_tools_.to_conservative(
            XT::Common::hstack(density_minus_one.evaluate(xx) + 1., 0., pressure_minus_04.evaluate(xx) + 0.4));
      },
      1);
  project(initial_values, *initial_values_);

  const double CFL = 0.15;
  this->inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right(300 / CFL, CFL, 0.04291);
}
