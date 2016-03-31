// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include "hyperbolic/fv-discretization.hh"

using namespace Dune;
using namespace Dune::GDT;

#if HAVE_EIGEN

TYPED_TEST_CASE(hyperbolic_FV_discretization_godunov_adaptiveRK, YaspGridTestCasesPartial);
TYPED_TEST(hyperbolic_FV_discretization_godunov_adaptiveRK, eoc_study_using_yaspgrid)
{
  this->eoc_study();
}

#else

TEST(DISABLED_hyperbolic_FV_discretization_godunov_adaptiveRK, eoc_study_using_yaspgrid)
{
}

#endif
