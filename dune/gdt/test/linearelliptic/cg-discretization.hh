// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH

#ifndef THIS_IS_A_BUILDBOT_BUILD
#define THIS_IS_A_BUILDBOT_BUILD 0
#endif

#include <dune/stuff/test/common.hh>

#include <dune/gdt/spaces/interface.hh>

#include "eocstudy.hh"
#include "discretizers/cg.hh"
#include "problems.hh"


template <class TestCaseType>
struct linearelliptic_CG_discretization : public ::testing::Test
{
  template <Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend>
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
                                          1> Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    Dune::Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
  } // ... eoc_study()
}; // linearelliptic_CG_discretization

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_CG_DISCRETIZATION_HH
