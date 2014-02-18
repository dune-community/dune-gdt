// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <dune/common/exceptions.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#else
#error This test requires ALUGrid!
#endif

#include <tuple>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/la/solver.hh>

#include "elliptic-testcases.hh"
#include "elliptic-cg-discretization.hh"
//#include "elliptic-sipdg-discretization.hh"
//#include "elliptic-swipdg-discretization.hh"

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;

// change this to toggle output
std::ostream& test_out = std::cout;
// std::ostream& test_out = DSC_LOG.devnull();

typedef testing::
    Types<std::tuple<EllipticTestCase::ESV07<AluConform2dGridType>, Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,
                     Dune::Stuff::LA::EigenDenseVector<double>>,
          std::tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,
                     Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>, Dune::Stuff::LA::EigenDenseVector<double>>,
          std::tuple<EllipticTestCase::ER07<AluConform2dGridType>, Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,
                     Dune::Stuff::LA::EigenDenseVector<double>>,
          std::tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,
                     Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>, Dune::Stuff::LA::EigenDenseVector<double>>,
          std::tuple<EllipticTestCase::Spe10Model1<AluConform2dGridType>,
                     Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>, Dune::Stuff::LA::EigenDenseVector<double>>>
        AluConform2dTestCases;

template <class TestTuple>
struct EllipticCGDiscretization : public ::testing::Test
{
  typedef typename std::tuple_element<0, TestTuple>::type TestCase;
  typedef typename std::tuple_element<1, TestTuple>::type MatrixType;
  typedef typename std::tuple_element<2, TestTuple>::type VectorType;

  typedef EllipticCG::Discretization<typename TestCase::GridPartType, 1, MatrixType, VectorType> DiscretizationType;

  void produces_correct_results() const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const TestCase test_case;
    const auto grid_part = test_case.reference_grid_part();

    DiscretizationType discretization(grid_part,
                                      test_case.boundary_info(),
                                      test_case.diffusion(),
                                      test_case.force(),
                                      test_case.dirichlet(),
                                      test_case.neumann());
    discretization.assemble();
    std::cout << " test case:   " << Stuff::Common::Typename<TestCase>::value() << std::endl;
    std::cout << " matrix type: " << Stuff::Common::Typename<MatrixType>::value() << std::endl;
    std::cout << " vector type: " << Stuff::Common::Typename<VectorType>::value() << std::endl;
    std::cout << " system size: " << discretization.system_matrix().rows() << "x"
              << discretization.system_matrix().cols() << std::endl;

    auto solution_vector = discretization.create_vector();
    auto tmp_vector = discretization.create_vector();
    Stuff::LA::Solver<MatrixType> linear_solver(discretization.system_matrix());

    // print table header
    const size_t min_first_column_width = std::string("solver option").size();
    size_t first_column_width = min_first_column_width;
    for (auto option : linear_solver.options())
      first_column_width = std::max(first_column_width, option.size());
    std::stringstream header;
    std::stringstream delimiter;
    header << " solver option";
    delimiter << "--------------";
    for (size_t ii = 0; ii <= first_column_width - min_first_column_width; ++ii) {
      header << " ";
      delimiter << "-";
    }
    header << "| time (s) | L^oo error (abs|rel) | thrown exception (see dune/stuff/solver.hh) ";
    delimiter << "+----------+----------------------+---------------------------------------------";
    std::cout << Stuff::Common::whitespaceify(header.str(), '=') << std::endl;
    std::cout << header.str() << std::endl;
    std::cout << delimiter.str() << std::endl;

    // loop over all available options
    size_t printed_rows = 0;
    for (auto option : linear_solver.options()) {
      if (option != "qr.sparse" && option.substr(0, 3) != "cg.") {
        const Stuff::Common::ConfigTree config = linear_solver.options(option);
        // print delimiter after every 3rd row
        if (printed_rows == 3) {
          std::cout << delimiter.str() << std::endl;
          printed_rows = 0;
        }
        // print
        std::cout << " " << option;
        for (size_t ii = 0; ii < first_column_width - option.size(); ++ii)
          std::cout << " ";
        // solve the system
        Dune::Timer timer;
        std::stringstream error_msg;
        bool success = true;
        try {
          linear_solver.apply(discretization.rhs_vector(), solution_vector, option);
        } catch (Stuff::Exceptions::linear_solver_failed_bc_matrix_did_not_fulfill_requirements&) {
          error_msg << Stuff::Common::colorStringRed("matrix_did_not_fulfill_requirements");
          if (option.substr(0, 3) == "cg.")
            error_msg << " (OK for this discretization)";
          success = false;
        } catch (Stuff::Exceptions::linear_solver_failed_bc_it_did_not_converge&) {
          error_msg << Stuff::Common::colorStringRed("did_not_converge");
        } catch (Stuff::Exceptions::linear_solver_failed_bc_it_was_not_set_up_correctly&) {
          error_msg << Stuff::Common::colorStringRed("was_not_set_up_correctly");
          success = false;
        } catch (Stuff::Exceptions::linear_solver_failed_bc_the_solution_does_not_solve_the_system&) {
          error_msg << Stuff::Common::colorStringRed("solution_does_not_solve_the_system");
        } catch (Stuff::Exceptions::linear_solver_failed&) {
          error_msg << Stuff::Common::colorStringRed("linear_solver_failed");
          success = false;
        }
        std::cout << " | " << std::setw(8) << std::setprecision(2) << std::fixed << timer.elapsed();
        std::cout << " | ";
        // test solution
        if (success) {
          discretization.system_matrix().mv(solution_vector, tmp_vector);
          tmp_vector -= discretization.rhs_vector();
          double threshhold = config.get<double>("post_check_solves_system");
          if (config.has_key("precision"))
            threshhold = config.get<double>("precision");
          std::stringstream absolute_error;
          absolute_error << std::setw(9) << std::setprecision(2) << std::scientific << tmp_vector.sup_norm();
          if (tmp_vector.sup_norm() < threshhold)
            std::cout << Stuff::Common::colorString(Stuff::Common::toString(absolute_error.str()),
                                                    Stuff::Common::Colors::green);
          else if (tmp_vector.sup_norm() < std::exp(0.5 * std::log(threshhold)))
            std::cout << Stuff::Common::colorString(Stuff::Common::toString(absolute_error.str()),
                                                    Stuff::Common::Colors::brown);
          else
            std::cout << Stuff::Common::colorStringRed(Stuff::Common::toString(absolute_error.str()));
          std::cout << " | " << std::setw(8) << std::setprecision(2) << std::scientific
                    << tmp_vector.sup_norm() / discretization.rhs_vector().sup_norm();
        } else
          std::cout << "                    ";
        std::cout << " | " << error_msg.str() << std::endl;
        ++printed_rows;
      }
    } // loop over all available options
  }
}; // EllipticCGDiscretization


#if 0
template< class TestCase >
struct EllipticSIPDGDiscretization
  : public ::testing::Test
{
  void produces_correct_results() const
  {
    if (std::is_same< TestCase, EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > > >::value) {
      std::cerr
          << Dune::Stuff::Common::colorStringRed("EllipticSIPDGDiscretization does not work for "
                                                 "EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > >!")
          << std::endl;
    } else {
      const TestCase test_case;
      test_case.print_header(test_out);
      test_out << std::endl;
      EllipticSIPDG::EocStudy< TestCase, 1 > eoc_study_1(test_case);
      auto errors_1 = eoc_study_1.run(test_out);
      for (const auto& norm : eoc_study_1.provided_norms()) {
        if (!Dune::Stuff::Common::FloatCmp::lt(errors_1[norm], eoc_study_1.expected_results(norm))) {
          std::stringstream ss;
          Dune::Stuff::Common::print(errors_1[norm],                     "errors           (" + norm + ")", ss);
          Dune::Stuff::Common::print(eoc_study_1.expected_results(norm), "   expected results (" + norm + ")", ss);
          DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
        }
      }
      test_out << std::endl;
      EllipticSIPDG::EocStudy< TestCase, 2 > eoc_study_2(test_case);
      auto errors_2 = eoc_study_2.run(test_out);
      for (const auto& norm : eoc_study_2.provided_norms())
        if (!Dune::Stuff::Common::FloatCmp::lt(errors_2[norm], eoc_study_2.expected_results(norm))) {
          std::stringstream ss;
          Dune::Stuff::Common::print(errors_2[norm],                     "errors           (" + norm + ")", ss);
          Dune::Stuff::Common::print(eoc_study_2.expected_results(norm), "   expected results (" + norm + ")", ss);
          DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
        }
    }
  }
}; // EllipticSIPDGDiscretization

template< class TestCase >
struct EllipticSWIPDGDiscretization
  : public ::testing::Test
{
  void produces_correct_results() const
  {
    const TestCase test_case;
    test_case.print_header(test_out);
    test_out << std::endl;
    EllipticSWIPDG::EocStudy< TestCase, 1 > eoc_study_1(test_case);
    auto errors_1 = eoc_study_1.run(test_out);
    for (const auto& norm : eoc_study_1.provided_norms()) {
      if (!Dune::Stuff::Common::FloatCmp::lt(errors_1[norm], eoc_study_1.expected_results(norm))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors_1[norm],                     "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study_1.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
      }
    }
    test_out << std::endl;
    EllipticSWIPDG::EocStudy< TestCase, 2 > eoc_study_2(test_case);
    auto errors_2 = eoc_study_2.run(test_out);
    for (const auto& norm : eoc_study_2.provided_norms())
      if (!Dune::Stuff::Common::FloatCmp::lt(errors_2[norm], eoc_study_2.expected_results(norm))) {
        std::stringstream ss;
        Dune::Stuff::Common::print(errors_2[norm],                     "errors           (" + norm + ")", ss);
        Dune::Stuff::Common::print(eoc_study_2.expected_results(norm), "   expected results (" + norm + ")", ss);
        DUNE_THROW_COLORFULLY(errors_are_not_as_expected, ss.str());
      }
  }
};
#endif

TYPED_TEST_CASE(EllipticCGDiscretization, AluConform2dTestCases);
TYPED_TEST(EllipticCGDiscretization, produces_correct_results)
{
  this->produces_correct_results();
}

// TYPED_TEST_CASE(EllipticSIPDGDiscretization, AluConform2dTestCases);
// TYPED_TEST(EllipticSIPDGDiscretization, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(EllipticSWIPDGDiscretization, AluConform2dTestCases);
// TYPED_TEST(EllipticSWIPDGDiscretization, produces_correct_results) {
//  this->produces_correct_results();
//}


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
