// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx>

#ifdef HAVE_FASP
#define HAVE_FASP 0
#endif

#include "elliptic-testcases.hh"
#include "elliptic-swipdg-discretization.hh"


template <class TestCase>
struct EllipticSWIPDGDiscretization : public ::testing::Test
{
  void estimator_study() const
  {
    try {
      const TestCase test_case;
      test_case.print_header(DSC_LOG_INFO);
      DSC_LOG_INFO << std::endl;
      EllipticSWIPDG::EstimatorStudy<TestCase> estimator_study(test_case);
      Dune::Stuff::Test::check_eoc_study_for_success(estimator_study, estimator_study.run(DSC_LOG_INFO));
    } catch (Dune::Stuff::Exceptions::spe10_data_file_missing& ee) {
      std::cerr << ee.what() << std::endl;
    }
  } // ... estimator_study()
}; // struct EllipticSWIPDGDiscretization


TYPED_TEST_CASE(EllipticSWIPDGDiscretization, EllipticEstimatorTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, estimator_study)
{
  this->estimator_study();
}
