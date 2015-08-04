// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include "problems/mixedboundary.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::MixedBoundaryTestCase< YaspGrid< 2, EquidistantOffsetCoordinates<double, 2> >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::MixedBoundaryTestCase< YaspGrid< 2, EquidistantOffsetCoordinates<double, 2> >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {2.86e-03, 6.32e-04};
      else
        return {2.96e-03, 7.44e-04, 1.81e-04, 3.95e-05};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {5.57e-02, 2.48e-02};
      else
        return {5.75e-02, 2.84e-02, 1.39e-02, 6.20e-03};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations< LinearElliptic::MixedBoundaryTestCase< YaspGrid< 2, EquidistantOffsetCoordinates<double, 2> >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;


} // namespace Tests
} // namespace GDT
} // namespace Dune
