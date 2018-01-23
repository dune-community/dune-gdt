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

#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/operators/advection-fv.hh>

#include "inviscid_compressible_flow__euler_1d_explicit.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class GL, class R, size_t m>
using FvSpaceType = FvSpace<GL, R, m, 1>;

template <class DF>
using FvOperator = GDT::AdvectionFvOperator<DF, typename DF::SpaceType::GridLayerType>;

using InviscidCompressibleFlowEuler1dExplicitFV = InviscidCompressibleFlowEuler1dExplicitTest<FvSpaceType, FvOperator>;


// TEST_F(InviscidCompressibleFlowEuler1dExplicitFV, periodic_boundaries)
//{
//  this->periodic_boundaries();
//}
// TEST_F(InviscidCompressibleFlowEuler1dExplicitFV, impermeable_walls_by_direct_euler_treatment)
//{
//  this->impermeable_walls_by_direct_euler_treatment();
//}
TEST_F(InviscidCompressibleFlowEuler1dExplicitFV, impermeable_walls_by_inviscid_mirror_treatment)
{
  this->impermeable_walls_by_inviscid_mirror_treatment();
}
// TEST_F(InviscidCompressibleFlowEuler1dExplicitFV,
//       inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right)
//{
//  this->inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right();
//}
