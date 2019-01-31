// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/stokes/stokes-taylorhood.hh>

using namespace Dune;
using namespace Dune::GDT::Test;

using StokesTest = StokesTestcase1<YASP_2D_EQUIDISTANT_OFFSET>;
// using StokesTest = StokesTestcase1<ALU_2D_SIMPLEX_CONFORMING>;
// using StokesTest = StokesTestcase1<UG_2D>;
TEST_F(StokesTest, run)
{
  this->run();
}
