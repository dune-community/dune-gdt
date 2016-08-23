// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH

#ifndef THIS_IS_A_BUILDBOT_BUILD
#define THIS_IS_A_BUILDBOT_BUILD 0
#endif

#include <dune/xt/common/test/common.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include "eocstudy.hh"
#include "discretizers/cg.hh"


template <class TestCaseType>
struct linearelliptic_CG_discretization : public ::testing::Test
{
  template <Dune::GDT::ChooseSpaceBackend space_backend, Dune::XT::LA::ChooseBackend la_backend>
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
#if THIS_IS_A_BUILDBOT_BUILD
    TestCaseType test_case(/*num_refs=*/1); // As in: only 1!
#else
    TestCaseType test_case;
#endif
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    typedef LinearElliptic::CGDiscretizer<typename TestCaseType::GridType,
                                          Stuff::Grid::ChooseLayer::level,
                                          space_backend,
                                          la_backend,
                                          1,
                                          typename TestCaseType::ProblemType::RangeFieldType,
                                          1>
        Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
    } catch (Dune::Stuff::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.cg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_CG_discretization

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH
