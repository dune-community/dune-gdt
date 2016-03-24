// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include "linearelliptic/swipdg-discretization.hh"

using namespace Dune;
using namespace Dune::GDT;


#if HAVE_DUNE_PDELAB && HAVE_DUNE_ISTL

TYPED_TEST_CASE(linearelliptic_SWIPDG_discretization, SGridTestCases);
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_pdelab_and_istl_and_sgrid) {
  this->template eoc_study< ChooseSpaceBackend::pdelab, Stuff::LA::ChooseBackend::istl_sparse >();
}

#else

TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_pdelab_and_istl_and_sgrid) {}

#endif
