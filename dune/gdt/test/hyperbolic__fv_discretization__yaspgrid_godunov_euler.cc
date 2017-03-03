// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include "hyperbolic/fv-discretization.hh"

using namespace Dune;
using namespace Dune::GDT;

#if HAVE_EIGEN

TYPED_TEST_CASE(hyperbolic_FV_discretization_godunov_euler, YaspGridTestCasesAll);
TYPED_TEST(hyperbolic_FV_discretization_godunov_euler, eoc_study_using_yaspgrid)
{
  this->eoc_study();
}

#else

TEST(DISABLED_hyperbolic_FV_discretization_godunov_euler, eoc_study_using_yaspgrid)
{
}

#endif
