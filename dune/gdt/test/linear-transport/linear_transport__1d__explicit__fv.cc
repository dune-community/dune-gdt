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

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/instationary-eocstudies/hyperbolic-nonconforming.hh>

#include "base.hh"

using namespace Dune;
using namespace Dune::GDT::Test;


using LinearTransport1dExplicitFvTest = LinearTransportExplicitTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_upwind_flux)
{
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__upwind", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_lax_riedrichs_flux)
{
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__lax_friedrichs", 0); // <- Only one of these!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "lax_friedrichs";
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_engquist_osher_flux)
{
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__engquist_osher", 0); // <- Only one of these!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "engquist_osher";
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_vijayasundaram_flux)
{
  this->visualization_steps_ = DXTC_CONFIG_GET("visualization_steps__vijayasundaram", 0); // <- Only one of these!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "vijayasundaram";
  this->run();
}
