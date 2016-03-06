// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/sodshocktube.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::ShockTubeTestCase< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >, double >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::ShockTubeTestCase< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >, double > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        return {6.55e-02, 4.00e-02};
      else
        return {6.76e-02, 4.19e-02, 2.61e-02, 1.61e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations< Hyperbolic::ShockTubeTestCase
                                                 < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >,
                                                   double >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1 >;


} // namespace Tests
} // namespace GDT
} // namespace Dune
