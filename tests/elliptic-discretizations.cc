// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#include <dune/common/exceptions.hh>

#ifndef HAVE_ALUGRID
static_assert(false, "This test requires alugrid!");
#endif
#include <dune/grid/alugrid.hh>

#include <dune/stuff/common/color.hh>

#undef HAVE_FASP

#include "elliptic-testcases.hh"
#include "elliptic-cg-discretization.hh"
#include "elliptic-sipdg-discretization.hh"
#include "elliptic-swipdg-discretization.hh"

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

// change this to toggle output
std::ostream& out = std::cout;
// std::ostream& out = DSC_LOG.devnull();

typedef testing::Types<EllipticTestCase::ESV07<AluConform2dGridType>,
                       EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,
                       EllipticTestCase::ER07<AluConform2dGridType>,
                       EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,
                       EllipticTestCase::Spe10Model1<AluConform2dGridType>> AluConform2dTestCases;

template <class TestCase>
struct EllipticCGDiscretization : public ::testing::Test
{
  void check() const
  {
    const TestCase test_case;
    test_case.print_header(out);
    out << std::endl;
    EllipticCG::EocStudy<TestCase, 1> eoc_study(test_case);
    auto errors = eoc_study.run(out);
    for (const auto& norm : eoc_study.provided_norms())
      if (errors[norm] > eoc_study.expected_results(norm))
        DUNE_THROW(errors_are_not_as_expected, "They really ain't (or you have to lower the expectations)!");
  }
};

TYPED_TEST_CASE(EllipticCGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticCGDiscretization, produces_correct_results)
{
  this->check();
}


template <class TestCase>
struct EllipticSIPDGDiscretization : public ::testing::Test
{
  void check() const
  {
    if (std::is_same<TestCase, EllipticTestCase::Spe10Model1<Dune::ALUConformGrid<2, 2>>>::value) {
      std::cerr << Dune::Stuff::Common::colorStringRed(
                       "EllipticSIPDGDiscretization does not work for "
                       + "EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > >!")
                << std::end;
    } else {
      const TestCase test_case;
      test_case.print_header(out);
      out << std::endl;
      size_t failure = 0;
      EllipticSIPDG::EocStudy<TestCase, 1> eoc_study_1(test_case);
      auto errors_1 = eoc_study_1.run(out);
      out << std::endl;
      EllipticSIPDG::EocStudy<TestCase, 2> eoc_study_2(test_case);
      auto errors_2 = eoc_study_2.run(out);
      for (const auto& norm : eoc_study_1.provided_norms())
        if (errors_1[norm] > eoc_study_1.expected_results(norm))
          ++failure;
      for (const auto& norm : eoc_study_2.provided_norms())
        if (errors_2[norm] > eoc_study_2.expected_results(norm))
          ++failure;
      if (failure)
        DUNE_THROW(errors_are_not_as_expected, "They really ain't (or you have to lower the expectations)!");
    }
  }
};

TYPED_TEST_CASE(EllipticSIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSIPDGDiscretization, produces_correct_results)
{
  this->check();
}


template <class TestCase>
struct EllipticSWIPDGDiscretization : public ::testing::Test
{
  void check() const
  {
    const TestCase test_case;
    test_case.print_header(out);
    out << std::endl;
    size_t failure = 0;
    EllipticSWIPDG::EocStudy<TestCase, 1> eoc_study_1(test_case);
    auto errors_1 = eoc_study_1.run(out);
    out << std::endl;
    EllipticSWIPDG::EocStudy<TestCase, 2> eoc_study_2(test_case);
    auto errors_2 = eoc_study_2.run(out);
    for (const auto& norm : eoc_study_1.provided_norms())
      if (errors_1[norm] > eoc_study_1.expected_results(norm))
        ++failure;
    for (const auto& norm : eoc_study_2.provided_norms())
      if (errors_2[norm] > eoc_study_2.expected_results(norm))
        ++failure;
    if (failure)
      DUNE_THROW(errors_are_not_as_expected, "They really ain't (or you have to lower the expectations)!");
  }
};

TYPED_TEST_CASE(EllipticSWIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results)
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
