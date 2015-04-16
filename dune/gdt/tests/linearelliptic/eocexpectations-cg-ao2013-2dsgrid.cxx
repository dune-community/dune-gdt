// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "problems/AO2013.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {7.42e-03, 2.09e-03, 5.39e-04, 1.40e-04};
    else if (type == "H1_semi")
      return {4.91e-01, 2.42e-01, 1.23e-01, 6.11e-02};
    else if (type == "energy")
      return {9.00e-02, 5.90e-02, 4.15e-02, 2.81e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
