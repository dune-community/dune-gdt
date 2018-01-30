// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>

#include "inviscid_compressible_flow__euler_1d__explicit.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class GL, class R, size_t m>
using DgOrder0SpaceType = DuneFemDgSpaceWrapper<GL, 0, R, m, 1>;

template <class DF>
using DgOperator = GDT::AdvectionDgOperator<DF, typename DF::SpaceType::GridLayerType>;

using InviscidCompressibleFlowEuler1dExplicitDgOrder0 =
    InviscidCompressibleFlowEuler1dExplicitTest<DgOrder0SpaceType, DgOperator, XT::Grid::Backends::part>;


TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder0, periodic_boundaries)
{
  this->periodic_boundaries();
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder0, impermeable_walls_by_direct_euler_treatment)
{
  this->impermeable_walls_by_direct_euler_treatment();
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder0, impermeable_walls_by_inviscid_mirror_treatment)
{
  this->impermeable_walls_by_inviscid_mirror_treatment();
}
TEST_F(InviscidCompressibleFlowEuler1dExplicitDgOrder0,
       inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right)
{
  this->inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right();
}
