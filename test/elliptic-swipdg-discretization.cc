// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#include "elliptic-swipdg-discretization.hh"

template <class TestCase>
struct EllipticSWIPDGDiscretization : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticSWIPDG::EocStudy<TestCase, 1> eoc_study_1(test_case);
    auto errors_1 = eoc_study_1.run(test_out);
    for (const auto& norm : eoc_study_1.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(
              errors_1[norm], truncate_vector(eoc_study_1.expected_results(norm), errors_1[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors_1[norm], "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study_1.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW(errors_are_not_as_expected, ss.str());
      }
    test_out << std::endl;
    EllipticSWIPDG::EocStudy<TestCase, 2> eoc_study_2(test_case);
    auto errors_2 = eoc_study_2.run(test_out);
    for (const auto& norm : eoc_study_2.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(
              errors_2[norm], truncate_vector(eoc_study_2.expected_results(norm), errors_2[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors_2[norm], "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study_2.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW(errors_are_not_as_expected, ss.str());
      }
  }
};


TYPED_TEST_CASE(EllipticSWIPDGDiscretization, EllipticTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}


#include <dune/stuff/test/test_main.cxx>
