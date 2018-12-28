// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)

#include <config.h>

#include "cg-mixedboundary-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
    LinearElliptic::ChooseDiscretizer::cg,
    1>::results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg,
                                                    1>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.21e-02, 2.58e-03};
#else
    return {5.00e-02, 1.25e-02, 3.06e-03, 6.57e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.29e-01, 6.07e-02};
#else
    return {2.58e-01, 1.32e-01, 6.63e-02, 3.13e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
