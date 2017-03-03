// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include "linearelliptic/swipdg-estimators.hh"
#include "linearelliptic/swipdg-estimator-testcases.hh"

using namespace Dune;
using namespace Dune::GDT;


#if 0 && HAVE_DUNE_FEM && HAVE_EIGEN && HAVE_DUNE_ALUGRID

TYPED_TEST_CASE(linearelliptic_SWIPDG_estimators, AluGridTestCases);
TYPED_TEST(linearelliptic_SWIPDG_estimators, eoc_study_using_fem_and_eigen_and_alugrid_order_1)
{
  this->template eoc_study<ChooseSpaceBackend::fem, XT::LA::Backends::eigen_sparse, 1>();
}

#else

TEST(DISABLED_linearelliptic_SWIPDG_estimators, eoc_study_using_fem_and_eigen_and_alugrid_order_1)
{
}
TEST(DISABLED_linearelliptic_SWIPDG_estimators, eoc_study_using_fem_and_eigen_and_alugrid_order_2)
{
}

#endif
