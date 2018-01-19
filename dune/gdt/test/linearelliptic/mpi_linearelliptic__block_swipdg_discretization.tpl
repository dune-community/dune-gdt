// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)

#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces.hh>
#include <dune/gdt/test/linearelliptic/eocstudy.hh>
#include <dune/gdt/test/linearelliptic/discretizers/block-ipdg.hh>
#include <dune/gdt/test/linearelliptic/problems.hh>


template <class TestCaseType, Dune::GDT::Backends SpaceBackend, Dune::XT::LA::Backends LaBackend>
struct linearelliptic_block_SWIPDG_discretization : public ::testing::Test
{
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
    TestCaseType test_case;
    test_case.print_header(DXTC_LOG_INFO_0);
    DXTC_LOG_INFO_0 << std::endl;
    typedef LinearElliptic::BlockIpdgDiscretizer<typename TestCaseType::GridType,
                                                 SpaceBackend,
                                                 LaBackend,
                                                 1,
                                                 typename TestCaseType::ProblemType::RangeFieldType,
                                                 1,
                                                 LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor,
                                                 LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DXTC_LOG_INFO_0));
    } catch (Dune::XT::Common::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.swipdg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_block_SWIPDG_discretization


using namespace Dune;
using namespace Dune::GDT;

// clang-format off
{% for TestCase, SpaceBackend, LaBackend, Name in config.permutations %} // clang-format on


typedef linearelliptic_block_SWIPDG_discretization<{{TestCase}},
                                                   Dune::GDT::Backends::{{SpaceBackend}},
                                                   Dune::XT::LA::Backends::{{LaBackend}}>
    linearelliptic_block_SWIPDG_discretization_{{Name}};

TEST_F(linearelliptic_block_SWIPDG_discretization_{{Name}}, eoc_study)
{
  this->eoc_study();
}

// clang-format off
{% endfor %} // clang-format on
