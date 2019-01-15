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

#include "ESV2007.hh"

using namespace Dune;
using namespace Dune::GDT::Test;


/**
 * It is quite hard to reproduce the first column of Table 1, since their integration was not too good and ours is by
 * now. So by tweaking the various polynomial degrees and integration orders (need to be low enough), one could get the
 * exact numbers (as we used to some years ago)...
 */
using ESV2007Table1Test = ESV2007DiffusionTest<ALU_2D_SIMPLEX_CONFORMING>;
TEST_F(ESV2007Table1Test, column_1)
{
  this->space_type_ = "dg_p1";
  const auto actual_results = this->run();
  const auto expected_results = DXTC_TEST_CONFIG_SUB("results");
  XT::Test::check_eoc_study_for_success(expected_results, actual_results);
}
