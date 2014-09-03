// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
# define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx>

#include "elliptic-cg-discretization.hh"


template< class TestCase >
struct EllipticCGDiscretization
  : public ::testing::Test
{
  void eoc_study() const
  {
    try {
      const TestCase test_case;
      test_case.print_header(DSC_LOG_INFO);
      DSC_LOG_INFO << std::endl;
      EllipticCG::EocStudy< TestCase, 1 > eoc_study(test_case);
      Dune::Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
    } catch (Dune::Stuff::Exceptions::spe10_data_file_missing& ee) {
      DSC_LOG_INFO << ee.what() << std::endl;
    }
  } // ... eoc_study(...)
}; // EllipticCGDiscretization


TYPED_TEST_CASE(EllipticCGDiscretization, EllipticTestCases);
TYPED_TEST(EllipticCGDiscretization, eoc_study) {
  this->eoc_study();
}
