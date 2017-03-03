// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include <dune/xt/common/test/common.hh>
#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/test/linearelliptic/eocstudy.hh>
#include <dune/gdt/test/linearelliptic/discretizers/cg.hh>

#include <dune/gdt/test/linearelliptic/problems.hh>
#include <dune/gdt/test/grids.hh>

template <class TestCaseType, Dune::GDT::ChooseSpaceBackend SpaceBackend, Dune::XT::LA::Backends LaBackend>
struct linearelliptic_CG_discretization : public ::testing::Test
{

  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
    TestCaseType test_case;
    test_case.print_header(DXTC_LOG_INFO);
    DXTC_LOG_INFO << std::endl;
    typedef LinearElliptic::CGDiscretizer<typename TestCaseType::GridType,
                                          XT::Grid::Layers::level,
                                          SpaceBackend,
                                          LaBackend,
                                          1,
                                          typename TestCaseType::ProblemType::RangeFieldType,
                                          1>
        Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DXTC_LOG_INFO));
    } catch (Dune::XT::Common::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.cg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_CG_discretization

// clang-format off
{% for TestCase, SpaceBackend, LaBackend, Name in config.permutations %}
// clang-format on

typedef linearelliptic_CG_discretization<{{TestCase}}, Dune::GDT::ChooseSpaceBackend::{{SpaceBackend}},
                                         Dune::XT::LA::Backends::{{LaBackend}}>
    linearelliptic_CG_discretization_{{Name}};

TEST_F(linearelliptic_CG_discretization_{{Name}}, eoc_study)
{
  this->eoc_study();
}
// clang-format off
{% endfor %}
// clang-format on
