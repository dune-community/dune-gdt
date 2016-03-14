// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/fokkerplanck/sourcebeam.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double>,
                                Hyperbolic::ChooseDiscretizer::fv, 1, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L1") {
      return {3.25e-01, 1.63e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
