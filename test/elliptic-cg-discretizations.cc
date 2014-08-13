// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#include "elliptic-cg-discretization.hh"


template< class TestCase >
struct EllipticCGDiscretization
  : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticCG::EocStudy< TestCase, 1 > eoc_study(test_case);
    auto errors = eoc_study.run(test_out);
    for (const auto& norm : eoc_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(errors[norm],
                                             truncate_vector(eoc_study.expected_results(norm), errors[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors[norm],                        "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW(errors_are_not_as_expected, ss.str());
      }
  }
}; // EllipticCGDiscretization


TYPED_TEST_CASE(EllipticCGDiscretization, EllipticTestCases);
TYPED_TEST(EllipticCGDiscretization, produces_correct_results) {
  this->produces_correct_results();
}


#include <dune/stuff/test/test_main.hh>
