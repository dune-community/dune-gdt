// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include <dune/gdt/test/linearelliptic/swipdg-estimators.hh>
#include <dune/gdt/test/linearelliptic/swipdg-estimator-testcases.hh>

using namespace Dune;
using namespace Dune::GDT;

// clang-format off
{% for TestCase, SpaceBackend, LaBackend, Name in config.permutations %}
// clang-format on

typedef linearelliptic_SWIPDG_estimators<{{TestCase}}>
    linearelliptic_SWIPDG_estimators_{{Name}};

TEST_F(linearelliptic_SWIPDG_estimators_{{Name}}, eoc_study)
{
  this->template eoc_study<Dune::GDT::Backends::{{SpaceBackend}}, Dune::XT::LA::Backends::{{LaBackend}}, 1>();
}

// clang-format off
{% endfor %}
// clang-format on
