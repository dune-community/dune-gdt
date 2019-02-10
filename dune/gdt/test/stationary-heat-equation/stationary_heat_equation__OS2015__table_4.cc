// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>

#include "OS2015.hh"

using namespace Dune;
using namespace Dune::GDT::Test;


using OS015Table4Test = OS2015MultiscaleTest<YASP_2D_EQUIDISTANT_OFFSET>;
TEST_F(OS015Table4Test, columns_1_3_and_5)
{
  this->space_type_ = "dg_p1";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
