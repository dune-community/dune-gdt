// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#undef HAVE_FASP

#include "elliptic-testcases.hh"
#include "elliptic-swipdg-discretization.hh"


template <class TestCase>
struct EllipticSWIPDGDiscretization : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    DSC_LOG_INFO << std::endl;
    EllipticSWIPDG::EstimatorStudy<TestCase> estimator_study(test_case);
    auto results = estimator_study.run(test_out);
    std::stringstream ss;
    for (const auto& norm : estimator_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(
              results[norm], truncate_vector(estimator_study.expected_results(norm), results[norm].size()))) {
        Dune::Stuff::Common::print(results[norm], "   errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(estimator_study.expected_results(norm), "   expected results (" + norm + ")", ss);
      }
    const std::string failure = ss.str();
    if (!failure.empty())
      DUNE_THROW(errors_are_not_as_expected, "\n" << failure);
  } // ... produces_correct_results()
}; // struct EllipticSWIPDGDiscretization


TYPED_TEST_CASE(EllipticSWIPDGDiscretization, EllipticEstimatorTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}
