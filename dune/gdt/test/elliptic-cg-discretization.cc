// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#ifndef HAVE_ALUGRID
static_assert(false, "This test requires alugrid!");
#endif

#include <dune/common/exceptions.hh>

#include <dune/grid/alugrid.hh>

#include "elliptic-testcases.hh"
#include "elliptic-cg-discretization.hh"

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

typedef testing::Types<EllipticTestCase::ESV07<AluConform2dGridType>> AluConform2dTestCases;

template <class TestType>
struct EllipticCGDiscretization : public ::testing::Test
{
  void check() const
  {
    // change this to toggle output
    //    std::ostream& out = std::cout;
    std::ostream& out = DSC_LOG.devnull();

    const TestType test_case;
    test_case.print_header(out);
    out << std::endl;
    EllipticCG::EocStudy<TestType, 1> eoc_study_1(test_case);
    auto errors_1 = eoc_study_1.run(out);
    for (const auto& norm : eoc_study_1.provided_norms())
      if (errors_1[norm] > eoc_study_1.expected_results(norm))
        DUNE_THROW(errors_are_not_as_expected, "They really are (or you have to lower the expectations)!");
  }
};

TYPED_TEST_CASE(EllipticCGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticCGDiscretization, ProducesCorrectResults)
{
  this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
