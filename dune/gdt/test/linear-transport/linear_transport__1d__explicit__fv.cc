// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

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
  this->visualization_steps_ =
      DXTC_TEST_CONFIG_GET("setup.visualization_steps", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_lax_riedrichs_flux)
{
  this->visualization_steps_ = DXTC_TEST_CONFIG_GET("setup.visualization_steps", 0);
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "lax_friedrichs";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_engquist_osher_flux)
{
  this->visualization_steps_ = DXTC_TEST_CONFIG_GET("setup.visualization_steps", 0);
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "engquist_osher";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_vijayasundaram_flux)
{
  this->visualization_steps_ = DXTC_TEST_CONFIG_GET("setup.visualization_steps", 0);
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "vijayasundaram";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
