// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#include <utility>

#include <dune/common/exceptions.hh>

#ifndef HAVE_ALUGRID
static_assert(false, "This test requires alugrid!");
#endif
#include <dune/grid/alugrid.hh>

#include "elliptic-testcases.hh"
#include "elliptic-cg-discretization.hh"
#include "elliptic-sipdg-discretization.hh"

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

// change this to toggle output
std::ostream& out = std::cout;
// std::ostream& out = DSC_LOG.devnull();

typedef testing::Types<std::pair<EllipticTestCase::ESV07<AluConform2dGridType>, Int<1>>> AluConform2dTestCases_Order1;

typedef testing::Types<std::pair<EllipticTestCase::ESV07<AluConform2dGridType>, Int<1>>,
                       std::pair<EllipticTestCase::ESV07<AluConform2dGridType>, Int<2>>,
                       std::pair<EllipticTestCase::ESV07<AluConform2dGridType>, Int<3>>>
    AluConform2dTestCases_HigherOrders;

template <class Pair>
struct EllipticCGDiscretization : public ::testing::Test
{
  typedef typename Pair::first_type TestType;
  static const unsigned int polOrder = Pair::second_type::value;

  void check() const
  {
    const TestType test_case;
    test_case.print_header(out);
    out << std::endl;
    EllipticCG::EocStudy<TestType, polOrder> eoc_study(test_case);
    auto errors = eoc_study.run(out);
    for (const auto& norm : eoc_study.provided_norms())
      if (errors[norm] > eoc_study.expected_results(norm))
        DUNE_THROW(errors_are_not_as_expected, "They really ain't (or you have to lower the expectations)!");
  }
};

TYPED_TEST_CASE(EllipticCGDiscretization, AluConform2dTestCases_Order1);
TYPED_TEST(EllipticCGDiscretization, produces_correct_results)
{
  this->check();
}


template <class Pair>
struct EllipticSIPDGDiscretization : public ::testing::Test
{
  typedef typename Pair::first_type TestType;
  static const unsigned int polOrder = Pair::second_type::value;

  void check() const
  {
    const TestType test_case;
    test_case.print_header(out);
    out << std::endl;
    EllipticSIPDG::EocStudy<TestType, polOrder> eoc_study(test_case);
    auto errors = eoc_study.run(out);
    for (const auto& norm : eoc_study.provided_norms())
      if (errors[norm] > eoc_study.expected_results(norm))
        DUNE_THROW(errors_are_not_as_expected, "They really ain't (or you have to lower the expectations)!");
  }
};

TYPED_TEST_CASE(EllipticSIPDGDiscretization, AluConform2dTestCases_HigherOrders);
TYPED_TEST(EllipticSIPDGDiscretization, produces_correct_results)
{
  this->check();
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
