// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include "hyperbolic/fv-discretization.hh"

using namespace Dune;
using namespace Dune::GDT;

#if HAVE_EIGEN

TYPED_TEST_CASE(hyperbolic_FV_discretization_godunovwithreconstruction_euler, YaspGridTestCasesLinear1D);
TYPED_TEST(hyperbolic_FV_discretization_godunovwithreconstruction_euler, eoc_study_using_yaspgrid)
{
#if HAVE_DUNE_XT_GRID
  this->eoc_study();
#else
  EXPECT_DEATH(this->eoc_study(), ".*");
#endif
}

#else

TEST(DISABLED_hyperbolic_FV_discretization_godunovwithreconstruction_euler, eoc_study_using_yaspgrid)
{
}

#endif
