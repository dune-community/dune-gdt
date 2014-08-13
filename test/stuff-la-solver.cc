// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <dune/common/exceptions.hh>

#include <tuple>

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

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
#include "elliptic-swipdg-discretization.hh"

// +-----------------------------------------------------------------------+
// | Global options. Can be used to disable output or enable slow solvers. |
// +-----------------------------------------------------------------------+

// change this to test all solvers (even really slow ones)
const bool test_all_solvers = false;

// +--------------------------------------------------+
// | 1st we define the test structs that do something |
// +--------------------------------------------------+

struct EllipticDiscretizations
{
  template <class DiscretizationType>
  static void run(const DiscretizationType& discretization, const std::string discretization_id, const bool sec = true)
  {
    try {
      using namespace Dune;
      using namespace Dune::GDT;

      typedef typename DiscretizationType::MatrixType MatrixType;

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
      if (sec) {
        header << "| time (s) ";
        delimiter << "+----------";
      } else {
        header << "| time (ms) ";
        delimiter << "+-----------";
      }
      header << "| L^oo error (abs|rel) | thrown exception (see dune/stuff/solver.hh) ";
      delimiter << "+----------------------+---------------------------------------------";

      Dune::Timer timer;
      discretization.assemble();

      test_out << Stuff::Common::whitespaceify(header.str(), '=') << std::endl;
      test_out << discretization_id << ", polorder " << discretization.polOrder << ", system size "
               << discretization.system_matrix().rows() << "x" << discretization.system_matrix().cols()
               << ", assembly took " << std::setprecision(2) << std::fixed << timer.elapsed() << "s" << std::endl;
      test_out << Stuff::Common::whitespaceify(header.str(), '-') << std::endl;
      test_out << header.str() << std::endl;
      test_out << delimiter.str() << std::endl;

      // loop over all available options
      size_t printed_rows = 0;
      for (auto option : linear_solver.options()) {
        // exclude some solvers that take too long to test
        if (test_all_solvers
            || !(option == "qr.sparse" || option == "bicgstab.identity" || option == "bicgstab.diagonal"
                 || option.substr(0, 3) == "cg.")) {
          const Stuff::Common::Configuration config = linear_solver.options(option);
          // print delimiter after every 3rd row
          if (printed_rows == 3) {
            test_out << delimiter.str() << std::endl;
            printed_rows = 0;
          }
          // print
          test_out << " " << option << std::flush;
          for (size_t ii = 0; ii < first_column_width - option.size(); ++ii)
            test_out << " ";
          // solve the system
          timer.reset();
          std::stringstream error_msg;
          bool success = true;
          try {
            linear_solver.apply(discretization.rhs_vector(), solution_vector, option);
          } catch (Stuff::Exceptions::linear_solver_failed_bc_matrix_did_not_fulfill_requirements&) {
            error_msg << Stuff::Common::colorStringRed("matrix_did_not_fulfill_requirements");
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
          if (sec)
            test_out << " | " << std::setw(8) << std::setprecision(2) << std::fixed << timer.elapsed();
          else
            test_out << " | " << std::setw(9) << std::setprecision(2) << std::fixed << 1000 * timer.elapsed();
          test_out << " | ";
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
              test_out << Stuff::Common::colorString(Stuff::Common::toString(absolute_error.str()),
                                                     Stuff::Common::Colors::green);
            else if (tmp_vector.sup_norm() < std::exp(0.5 * std::log(threshhold)))
              test_out << Stuff::Common::colorString(Stuff::Common::toString(absolute_error.str()),
                                                     Stuff::Common::Colors::brown);
            else
              test_out << Stuff::Common::colorStringRed(Stuff::Common::toString(absolute_error.str()));
            test_out << " | " << std::setw(8) << std::setprecision(2) << std::scientific
                     << tmp_vector.sup_norm() / discretization.rhs_vector().sup_norm();
          } else
            test_out << "                    ";
          test_out << " | " << error_msg.str() << std::endl;
          ++printed_rows;
        }
      } // loop over all available options
    } catch (Dune::IOError&) {
    } // <- the SPE10Model1 function needs a data file on disc
  } // ... run(...)
}; // struct EllipticDiscretizations


template <class TestTuple>
struct SmallEllipticSystems : public ::testing::Test, EllipticDiscretizations
{
  typedef typename std::tuple_element<0, TestTuple>::type TestCase;
  typedef typename std::tuple_element<1, TestTuple>::type MatrixType;
  typedef typename std::tuple_element<2, TestTuple>::type VectorType;

  typedef typename TestCase::GridViewType GridViewType;
  typedef typename TestCase::GridPartType GridPartType;

  void produces_correct_results() const
  {
    using namespace Dune;

    test_out << " test case:   " << Stuff::Common::Typename<TestCase>::value() << std::endl;
    test_out << " matrix type: " << Stuff::Common::Typename<MatrixType>::value() << std::endl;
    test_out << " vector type: " << Stuff::Common::Typename<VectorType>::value() << std::endl;

    const TestCase test_case;
    const auto grid_view = test_case.level_grid_view(1);
    const auto grid_part = test_case.level_grid_part(1);

    run(EllipticCG::Discretization<GridViewType, 1, MatrixType, VectorType>(grid_view,
                                                                            test_case.boundary_info(),
                                                                            test_case.diffusion(),
                                                                            test_case.force(),
                                                                            test_case.dirichlet(),
                                                                            test_case.neumann()),
        "continuous Galerkin",
        false);

    run(EllipticSWIPDG::Discretization<GridPartType, 1, MatrixType, VectorType>(grid_part,
                                                                                test_case.boundary_info(),
                                                                                test_case.diffusion(),
                                                                                test_case.force(),
                                                                                test_case.dirichlet(),
                                                                                test_case.neumann()),
        "SWIP discontinuous Galerkin",
        false);
  } // ... produces_correct_results()
}; // SmallEllipticSystems


template <class TestTuple>
struct DISABLED_LargeEllipticSystems : public ::testing::Test, EllipticDiscretizations
{
  typedef typename std::tuple_element<0, TestTuple>::type TestCase;
  typedef typename std::tuple_element<1, TestTuple>::type MatrixType;
  typedef typename std::tuple_element<2, TestTuple>::type VectorType;

  typedef typename TestCase::GridViewType GridViewType;
  typedef typename TestCase::GridPartType GridPartType;

  void produces_correct_results() const
  {
    using namespace Dune;

    test_out << " test case:   " << Stuff::Common::Typename<TestCase>::value() << std::endl;
    test_out << " matrix type: " << Stuff::Common::Typename<MatrixType>::value() << std::endl;
    test_out << " vector type: " << Stuff::Common::Typename<VectorType>::value() << std::endl;

    const TestCase test_case;
    const auto grid_view = test_case.reference_grid_view();
    const auto grid_part = test_case.reference_grid_part();

    run(EllipticCG::Discretization<GridViewType, 1, MatrixType, VectorType>(grid_view,
                                                                            test_case.boundary_info(),
                                                                            test_case.diffusion(),
                                                                            test_case.force(),
                                                                            test_case.dirichlet(),
                                                                            test_case.neumann()),
        "continuous Galerkin");

    run(EllipticSWIPDG::Discretization<GridPartType, 2, MatrixType, VectorType>(grid_part,
                                                                                test_case.boundary_info(),
                                                                                test_case.diffusion(),
                                                                                test_case.force(),
                                                                                test_case.dirichlet(),
                                                                                test_case.neumann()),
        "SWIP discontinuous Galerkin");
  } // ... produces_correct_results()
}; // LargeEllipticSystems

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

#if HAVE_ALUGRID

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> AluConform2dGridType;

#define ALU_CONFORM_2D_COMMONDENSE_TEST_CASES                                                                                                                                                                                                     \
  /*std::tuple< EllipticTestCase::ESV07< AluConform2dGridType >,*/                                                                                                                                                                                \
  /*Dune::Stuff::LA::CommonDenseMatrix< double >,*/                                                                                                                                                                                               \
  /*Dune::Stuff::LA::CommonDenseVector< double > >*/                                                                                                                                                                                              \
  /*,*/ std::                                                                                                                                                                                                                                     \
      tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                                                                                                                                                            \
            Dune::Stuff::LA::CommonDenseMatrix<double>,                                                                                                                                                                                           \
            Dune::Stuff::LA::                                                                                                                                                                                                                     \
                CommonDenseVector<double>> /*, std::tuple< EllipticTestCase::ER07<                                           \
                                                                AluConform2dGridType >,*/ /*Dune::Stuff::LA::CommonDenseMatrix<                \
                                                                                             double >,*/ /*Dune::Stuff::LA::CommonDenseVector<                       \
                                                                                                                                                                                         double >                                                 \
                                                                                                                                                                                         >*/                                                      \
      ,                                                                                                                                                                                                                                           \
      std::                                                                                                                                                                                                                                       \
          tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                                                                                                                                                       \
                Dune::Stuff::LA::CommonDenseMatrix<double>,                                                                                                                                                                                       \
                Dune::Stuff::LA::                                                                                                                                                                                                                 \
                    CommonDenseVector<double>> /*, std::tuple< EllipticTestCase::Spe10Model1<                                             \
                                                       AluConform2dGridType >,*/ /*Dune::Stuff::LA::CommonDenseMatrix<                         \
                                                                                    double >,*/ /*Dune::Stuff::LA::CommonDenseVector< \
                                                                                                                                                                                                               double > >*/

#define ALU_CONFORM_2D_EIGENDENSE_TEST_CASES                                                                                                                                                                                            \
  /*std::tuple< EllipticTestCase::ESV07< AluConform2dGridType >,*/                                                                                                                                                                      \
  /*Dune::Stuff::LA::EigenDenseMatrix< double >,*/                                                                                                                                                                                      \
  /*Dune::Stuff::LA::EigenDenseVector< double > >*/                                                                                                                                                                                     \
  /*,*/ std::                                                                                                                                                                                                                           \
      tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                                                                                                                                                  \
            Dune::Stuff::LA::EigenDenseMatrix<double>,                                                                                                                                                                                  \
            Dune::Stuff::LA::                                                                                                                                                                                                           \
                EigenDenseVector<double>> /*, std::tuple< EllipticTestCase::ER07<                                          \
                                                               AluConform2dGridType >,*/ /*Dune::Stuff::LA::EigenDenseMatrix<                \
                                                                                            double >,*/ /*Dune::Stuff::LA::EigenDenseVector<                 \
                                                                                                                                                                                      double >                                          \
                                                                                                                                                                                      >*/                                               \
      ,                                                                                                                                                                                                                                 \
      std::                                                                                                                                                                                                                             \
          tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                                                                                                                                             \
                Dune::Stuff::LA::EigenDenseMatrix<double>,                                                                                                                                                                              \
                Dune::Stuff::LA::                                                                                                                                                                                                       \
                    EigenDenseVector<double>> /*, std::tuple< EllipticTestCase::Spe10Model1<                                         \
                                                         AluConform2dGridType >,*/ /*Dune::Stuff::LA::EigenDenseMatrix<                      \
                                                                                      double >,*/ /*Dune::Stuff::LA::EigenDenseVector< \
                                                                                                                                                                                                      double > >*/

#define ALU_CONFORM_2D_EIGENSPARSE_TEST_CASES                                                                          \
  std::tuple<EllipticTestCase::ESV07<AluConform2dGridType>,                                                            \
             Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                       \
             Dune::Stuff::LA::EigenDenseVector<double>>,                                                               \
      std::tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                            \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::ER07<AluConform2dGridType>,                                                         \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                           \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::Spe10Model1<AluConform2dGridType>,                                                  \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>

#define ALU_CONFORM_2D_ISTLSPARSE_TEST_CASES                                                                           \
  std::tuple<EllipticTestCase::ESV07<AluConform2dGridType>,                                                            \
             Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                        \
             Dune::Stuff::LA::IstlDenseVector<double>>,                                                                \
      std::tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                            \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::ER07<AluConform2dGridType>,                                                         \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                           \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::Spe10Model1<AluConform2dGridType>,                                                  \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>

#define ISTL_EIGEN_COMPARISON                                                                                          \
  std::tuple<EllipticTestCase::ESV07<AluConform2dGridType>,                                                            \
             Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                        \
             Dune::Stuff::LA::IstlDenseVector<double>>,                                                                \
      std::tuple<EllipticTestCase::ESV07<AluConform2dGridType>,                                                        \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                            \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::LocalThermalBlock<AluConform2dGridType>,                                            \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::ER07<AluConform2dGridType>,                                                         \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::ER07<AluConform2dGridType>,                                                         \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                           \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::MixedBoundaryTypes<AluConform2dGridType>,                                           \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>,                                                           \
      std::tuple<EllipticTestCase::Spe10Model1<AluConform2dGridType>,                                                  \
                 Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>,                                                    \
                 Dune::Stuff::LA::IstlDenseVector<double>>,                                                            \
      std::tuple<EllipticTestCase::Spe10Model1<AluConform2dGridType>,                                                  \
                 Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,                                                   \
                 Dune::Stuff::LA::EigenDenseVector<double>>

#endif // HAVE_ALUGRID


typedef testing::Types<
#if HAVE_ALUGRID
    ALU_CONFORM_2D_COMMONDENSE_TEST_CASES
#if HAVE_EIGEN
    ,
    ALU_CONFORM_2D_EIGENDENSE_TEST_CASES, ALU_CONFORM_2D_EIGENSPARSE_TEST_CASES
#endif // HAVE_EIGEN
#if HAVE_ISTL
    ,
    ALU_CONFORM_2D_ISTLSPARSE_TEST_CASES
#endif // HAVE_ISTL
#endif // HAVE_ALUGRID
    > Small_TestCases;

typedef testing::Types<
#if HAVE_ALUGRID
    ALU_CONFORM_2D_ISTLSPARSE_TEST_CASES
#if HAVE_EIGEN
    ,
    ALU_CONFORM_2D_EIGENSPARSE_TEST_CASES
#endif // HAVE_EIGEN
#endif // HAVE_ALUGRID
    > Large_TestCases;

// typedef testing::Types< ISTL_EIGEN_COMPARISON
//                      > Large_TestCases;

// +--------------------------------------------------------------------------------------+
// | 3rd we combine the test structs with their appropriate arguments to create the tests |
// +--------------------------------------------------------------------------------------+

TYPED_TEST_CASE(SmallEllipticSystems, Small_TestCases);
TYPED_TEST(SmallEllipticSystems, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(DISABLED_LargeEllipticSystems, Large_TestCases);
TYPED_TEST(DISABLED_LargeEllipticSystems, produces_correct_results)
{
  this->produces_correct_results();
}

// +--------------------------------------------------------------------------------------+
// | 4th we run all the tests                                                             |
// | (run the resulting executable with '--gtest_catch_exceptions=0' to see an exception) |
// +--------------------------------------------------------------------------------------+

#include <dune/stuff/test/test_main.hh>
