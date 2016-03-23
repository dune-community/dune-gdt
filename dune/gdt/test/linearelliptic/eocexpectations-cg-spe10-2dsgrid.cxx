// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "problems/spe10.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.86e-02, 1.51e-02};
    else if (type == "H1_semi")
      return {3.31e-01, 4.32e-01};
    else if (type == "energy")
      return {9.58e-01, 1.37e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune
