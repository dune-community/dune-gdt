// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#undef HAVE_FASP

#include <dune/common/exceptions.hh>

#if ENABLE_ALUGRID
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>

#include "elliptic-testcases.hh"
#include "elliptic-swipdg-discretization.hh"

// change this to toggle test_output
std::ostream& test_out = std::cout;
//std::ostream& test_out = DSC_LOG.devnull();

typedef Dune::ALUConformGrid< 2, 2 > AluConform2dGridType;

typedef testing::Types< EllipticTestCase::ESV07< AluConform2dGridType >
//                      , EllipticTestCase::LocalThermalBlock< AluConform2dGridType >
//                      , EllipticTestCase::Spe10Model1< AluConform2dGridType >
                      > AluConform2dTestCases;


template< class TestCase >
struct EllipticSWIPDGDiscretization
  : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticSWIPDG::EstimatorStudy< TestCase > estimator_study(test_case);
    auto results = estimator_study.run(test_out);
    std::stringstream ss;
    for (const auto& norm : estimator_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(results[norm],
                                             truncate_vector(estimator_study.expected_results(norm),
                                                             results[norm].size()))) {
        Dune::Stuff::Common::print(results[norm],                          "   errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(estimator_study.expected_results(norm), "   expected results (" + norm + ")", ss);
      }
    const std::string failure = ss.str();
    if (!failure.empty()) DUNE_THROW_COLORFULLY(errors_are_not_as_expected, "\n" << failure);
  } // ... produces_correct_results()
}; // struct EllipticSWIPDGDiscretization


TYPED_TEST_CASE(EllipticSWIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results) {
  this->produces_correct_results();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}

#else // ENABLE_ALUGRID
#warning "nothing tested in elliptic-estimators.cc because alugrid is missing"
int main(int, char**)
{
  return 0;
}
#endif //ENABLE_ALUGRID
