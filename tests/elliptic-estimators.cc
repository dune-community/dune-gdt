// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#undef HAVE_FASP

#include <dune/common/exceptions.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#else
#error This test requires ALUGrid!
#endif

#include "elliptic-testcases.hh"
#include "elliptic-swipdg-discretization.hh"

// change this to toggle test_output
std::ostream& test_out = std::cout;
// std::ostream& test_out = DSC_LOG.devnull();


class errors_are_not_as_expected : public Dune::Exception
{
};

std::vector<double> truncate_vector(const std::vector<double>& in, const size_t size)
{
  assert(size <= in.size());
  if (size == in.size())
    return in;
  else {
    std::vector<double> ret(size);
    for (size_t ii = 0; ii < size; ++ii)
      ret[ii] = in[ii];
    return ret;
  }
} // ... truncate_vector(...)


typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

typedef testing::Types<EllipticTestCase::ESV07<AluConform2dGridType>> AluConform2dTestCases;


template <class TestCase>
struct EllipticSWIPDGDiscretization : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticSWIPDG::EstimatorStudy<TestCase> estimator_study(test_case);
    auto results = estimator_study.run(test_out);
    for (const auto& norm : estimator_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(
              results[norm], truncate_vector(estimator_study.expected_results(norm), results[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(results[norm], "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(estimator_study.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
      }
  } // ... produces_correct_results()
}; // struct EllipticSWIPDGDiscretization


TYPED_TEST_CASE(EllipticSWIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results)
{
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
