// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "problems/mixedboundary.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {7.38e-02, 1.84e-02};
      else
        return {7.56e-02, 2.06e-02, 5.56e-03, 1.36e-03};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {3.66e-01, 1.72e-01};
      else
        return {3.76e-01, 1.95e-01, 9.89e-02, 4.67e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune
