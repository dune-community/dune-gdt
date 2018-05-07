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

#include <dune/xt/common/test/common.hh>
#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/test/linearelliptic/eocstudy.hh>
#include <dune/gdt/test/linearelliptic/discretizers/ipdg.hh>

template <class TestCaseType, Dune::GDT::Backends SpaceBackend, Dune::XT::LA::Backends LaBackend>
struct linearelliptic_SWIPDG_discretization : public ::testing::Test
{
  template <int polOrder>
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
    TestCaseType test_case;
    test_case.print_header(DXTC_LOG_INFO_0);
    DXTC_LOG_INFO_0 << std::endl;
    typedef LinearElliptic::IpdgDiscretizer<typename TestCaseType::GridType,
                                            TestCaseType::layer_type,
                                            SpaceBackend,
                                            LaBackend,
                                            polOrder,
                                            typename TestCaseType::ProblemType::RangeFieldType,
                                            1,
                                            LocalEllipticIpdgIntegrands::Method::swipdg>
        Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DXTC_LOG_INFO_0));
    } catch (Dune::XT::Functions::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.swipdg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_SWIPDG_discretization


using namespace Dune;
using namespace Dune::GDT;

// clang-format off
{% for TestCase, SpaceBackend, LaBackend, Name in config.permutations %}

typedef linearelliptic_SWIPDG_discretization<{{TestCase}}, Dune::GDT::Backends::{{SpaceBackend}},
                                             Dune::XT::LA::Backends::{{LaBackend}}>
    linearelliptic_SWIPDG_discretization{{Name}};

TEST_F(linearelliptic_SWIPDG_discretization{{Name}}, eoc_study_order_1)
{
  this->template eoc_study<1>();
}
TEST_F(linearelliptic_SWIPDG_discretization{{Name}}, eoc_study_order_2)
{
  this->template eoc_study<2>();
}
{% endfor %}
// clang-format on
