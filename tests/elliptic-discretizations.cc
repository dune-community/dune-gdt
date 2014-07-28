// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS

// This one has to come first (includes the config.h)!
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/stuff/test/test_common.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#ifdef HAVE_FASP
#undef HAVE_FASP
#endif

#include <dune/common/exceptions.hh>

#if ENABLE_ALUGRID
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/alugrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>


#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "elliptic-testcases.hh"
#include "elliptic-cg-discretization.hh"
#include "elliptic-sipdg-discretization.hh"
#include "elliptic-swipdg-discretization.hh"

// change this to toggle output
std::ostream& test_out = std::cout;
// std::ostream& test_out = DSC_LOG.devnull();

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

typedef testing::Types<EllipticTestCase::ESV07<AluConform2dGridType>,
                       EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,
                       EllipticTestCase::ER07<AluConform2dGridType>,
                       EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,
                       EllipticTestCase::Spe10Model1<AluConform2dGridType>> AluConform2dTestCases;


template <class TestCase>
struct EllipticCGDiscretization : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticCG::EocStudy<TestCase, 1> eoc_study(test_case);
    auto errors = eoc_study.run(test_out);
    for (const auto& norm : eoc_study.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(errors[norm],
                                             truncate_vector(eoc_study.expected_results(norm), errors[norm].size()))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors[norm], "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
      }
  }
}; // EllipticCGDiscretization

template <class TestCase>
struct EllipticSIPDGDiscretization : public ::testing::Test
{
  void produces_correct_results() const
  {
    if (std::is_same<TestCase, EllipticTestCase::Spe10Model1<Dune::ALUConformGrid<2, 2>>>::value) {
      std::cerr << Dune::Stuff::Common::colorStringRed("EllipticSIPDGDiscretization does not work for "
                                                       "EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > >!")
                << std::endl;
    } else {
      const TestCase test_case;
      test_case.print_header(test_out);
      test_out << std::endl;
      EllipticSIPDG::EocStudy<TestCase, 1> eoc_study_1(test_case);
      auto errors_1 = eoc_study_1.run(test_out);
      for (const auto& norm : eoc_study_1.provided_norms())
        if (!Dune::Stuff::Common::FloatCmp::lt(
                errors_1[norm], truncate_vector(eoc_study_1.expected_results(norm), errors_1[norm].size()))) {
          std::stringstream ss;
          Dune::Stuff::Common::print(errors_1[norm], "errors           (" + norm + ")", ss);
          Dune::Stuff::Common::print(eoc_study_1.expected_results(norm), "   expected results (" + norm + ")", ss);
          DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
        }
      test_out << std::endl;
      EllipticSIPDG::EocStudy<TestCase, 2> eoc_study_2(test_case);
      auto errors_2 = eoc_study_2.run(test_out);
      for (const auto& norm : eoc_study_2.provided_norms())
        if (!Dune::Stuff::Common::FloatCmp::lt(
                errors_2[norm], truncate_vector(eoc_study_2.expected_results(norm), errors_2[norm].size()))) {
          std::stringstream ss;
          Dune::Stuff::Common::print(errors_2[norm], "errors           (" + norm + ")", ss);
          Dune::Stuff::Common::print(eoc_study_2.expected_results(norm), "   expected results (" + norm + ")", ss);
          DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
        }
    }
  }
}; // EllipticSIPDGDiscretization

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
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
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
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
      }
  }
};


TYPED_TEST_CASE(EllipticSIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSIPDGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(EllipticSWIPDGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(EllipticCGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticCGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}

#else // ENABLE_ALUGRID
#warning "nothing tested in elliptic-discretizations.cc because alugrid is missing"
int main(int, char**)
{
  return 0;
}
#endif // ENABLE_ALUGRID
