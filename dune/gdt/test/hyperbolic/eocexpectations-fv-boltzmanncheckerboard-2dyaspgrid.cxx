// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/2dboltzmann.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Dune::YaspGrid<2>, double, 1>,
                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::Boltzmann2DCheckerboardTestCase<Dune::YaspGrid<2>, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (DSC::FloatCmp::eq(test_case.t_end(), 3.2))
        return {2.27e+00, 1.09e+00};
      else if (DSC::FloatCmp::eq(test_case.t_end(), 3.2 / 5.0))
        return {1.12e-01, 6.11e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Dune::YaspGrid<2>, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 2,
                                         Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
