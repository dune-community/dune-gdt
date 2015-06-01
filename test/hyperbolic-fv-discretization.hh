// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef THIS_IS_A_BUILDBOT_BUILD
# define THIS_IS_A_BUILDBOT_BUILD 0
#endif

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/tests/hyperbolic/eocstudy.hh>
#include <dune/gdt/tests/hyperbolic/discretizers/fv.hh>

#include "hyperbolic-testcases.hh"


template< class TestCaseType >
struct hyperbolic_FV_discretization
  : public ::testing::Test
{
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
#if THIS_IS_A_BUILDBOT_BUILD
    TestCaseType test_case(/*num_refs = */ 1);
#else
    TestCaseType test_case(/*num_refs = */ 3);
#endif
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    typedef Hyperbolic::FVDiscretizer< typename TestCaseType::GridType, double, 1, true >  Discretizer;
    Tests::HyperbolicEocStudy< TestCaseType, Discretizer > eoc_study(test_case, {}, "test");
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
  } // ... eoc_study()

}; // hyperbolic_FV_discretization
