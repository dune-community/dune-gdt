// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx>

#include "elliptic-sipdg-discretization.hh"


template <class TestCase>
struct EllipticSIPDGDiscretization : public ::testing::Test
{
  template <int polOrder>
  void eoc_study() const
  {
    if (std::is_same<TestCase,
                     EllipticTestCase::Spe10Model1<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>::value)
      DSC_LOG_INFO << DSC::colorStringRed("EllipticSIPDGDiscretization does not work for "
                                          "EllipticTestCase::Spe10Model1< Dune::ALUGrid< 2, 2, Dune::simplex, "
                                          "Dune::conforming > >!")
                   << std::endl;
    else {
      const TestCase test_case;
      test_case.print_header(DSC_LOG_INFO);
      DSC_LOG_INFO << std::endl;
      EllipticSIPDG::EocStudy<TestCase, polOrder> eoc_study(test_case);
      Dune::Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
    }
  } // ... eoc_study(...)
}; // EllipticSIPDGDiscretization


TYPED_TEST_CASE(EllipticSIPDGDiscretization, EllipticTestCases);
TYPED_TEST(EllipticSIPDGDiscretization, eoc_study_polorder_1)
{
  this->template eoc_study<1>();
}
TYPED_TEST(EllipticSIPDGDiscretization, eoc_study_polorder_2)
{
  this->template eoc_study<2>();
}

TEST(DISABLED_EllipticSIPDGDiscretization, eoc_study_polorder_1_ESV07_AluConform2d)
{
}
TEST(DISABLED_EllipticSIPDGDiscretization, eoc_study_polorder_2_ESV07_AluConform2d)
{
}
