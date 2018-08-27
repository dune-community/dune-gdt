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

#include "base.hh"

using namespace Dune;
using namespace Dune::GDT;
using namespace Dune::GDT::Test;


using InviscidCompressibleFlow1dEulerExplicitFvTest =
    InviscidCompressibleFlowEulerExplicitTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest, periodic_boundaries)
{
  // Something like 100 to visualize!
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__periodic_boundaries", 0);
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "vijayasundaram";
  this->run();
}
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest, impermeable_walls_by_direct_euler_treatment)
{
  // Only use one visualization at a time!
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__impermeable_walls_by_direct_euler_treatment", 0);
  this->space_type_ = "fv";
  this->boundary_treatment = "impermeable_walls_by_direct_euler_treatment";
  this->run();
}
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest, impermeable_walls_by_inviscid_mirror_treatment)
{
  // Only use one visualization at a time!
  this->visualization_steps_ =
      DXTC_CONFIG_GET("visualization_steps__impermeable_walls_by_inviscid_mirror_treatment", 0);
  this->space_type_ = "fv";
  this->boundary_treatment = "impermeable_walls_by_inviscid_mirror_treatment";
  this->run();
}
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest,
       inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right)
{
  // Only use one visualization at a time!
  this->visualization_steps_ = DXTC_CONFIG_GET(
      "visualization_steps__inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right", 0);
  this->T_end_ = 2; // We need more time to hit the right wall
  this->space_type_ = "fv";
  this->boundary_treatment = "inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right";
  this->run();
}
