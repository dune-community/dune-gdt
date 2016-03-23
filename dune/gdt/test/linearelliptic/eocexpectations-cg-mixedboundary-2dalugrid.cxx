// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "problems/mixedboundary.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                          1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {7.68e-02, 1.92e-02};
      else
        return {8.04e-02, 2.35e-02, 5.95e-03, 1.28e-03};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {3.23e-01, 1.59e-01};
      else
        return {3.35e-01, 1.83e-01, 9.18e-02, 4.14e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex, nonconforming>, double,
                                                                          1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {4.90e-02, 1.09e-02};
      else
        return {5.18e-02, 1.39e-02, 3.42e-03, 7.09e-04};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {2.47e-01, 1.16e-01};
      else
        return {2.55e-01, 1.34e-01, 6.63e-02, 2.98e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
