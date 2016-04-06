// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/burgers.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 2 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 2,
                                 Hyperbolic::FluxTimeStepperCombinations::godunov_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 2 >
{
  typedef Hyperbolic::BurgersTestCase< Dune::YaspGrid< 2 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1")
      if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
        return {1.79e-03, 8.44e-04};
      else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0/5.0))
        return {3.60e-04, 1.70e-04};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 2 >,double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          2,
                                          Hyperbolic::FluxTimeStepperCombinations::godunov_euler >;


} // namespace Tests
} // namespace GDT
} // namespace Dune
