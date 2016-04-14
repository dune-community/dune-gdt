// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "../problems/mixedboundary.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {9.48e-02, 2.13e-02};
      else
        return {9.77e-02, 2.49e-02, 6.16e-03, 1.38e-03};
    } else if (type == "H1_semi" || type == "energy") {
      if (test_case.num_refinements() == 1)
        return {3.64e-01, 1.67e-01};
      else
        return {3.71e-01, 1.87e-01, 9.31e-02, 4.31e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune
