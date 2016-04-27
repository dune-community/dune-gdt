// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include "linearelliptic/swipdg-discretization.hh"
#include "linearelliptic/swipdg-testcases.hh"

using namespace Dune;
using namespace Dune::GDT;


#if HAVE_DUNE_FEM && HAVE_DUNE_ISTL && HAVE_DUNE_ALUGRID

TYPED_TEST_CASE(linearelliptic_SWIPDG_discretization, AluGridTestCases);
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl_and_alugrid_order_1)
{
  this->template eoc_study<ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse, 1>();
}
TYPED_TEST(linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl_and_alugrid_order_2)
{
  this->template eoc_study<ChooseSpaceBackend::fem, Stuff::LA::ChooseBackend::istl_sparse, 2>();
}

#else

TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl_and_alugrid_order_1)
{
}
TEST(DISABLED_linearelliptic_SWIPDG_discretization, eoc_study_using_fem_and_istl_and_alugrid_order_2)
{
}

#endif
