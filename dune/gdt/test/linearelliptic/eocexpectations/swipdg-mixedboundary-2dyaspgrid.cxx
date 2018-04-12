// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#include <config.h>

#include "swipdg-mixedboundary-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.19e-02, 2.64e-03};
#else
    return {7.56e-02, 2.06e-02, 5.56e-03, 1.36e-03};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.35e-01, 6.47e-02};
#else
    return {3.76e-01, 1.95e-01, 9.89e-02, 4.67e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {7.66e-04, 2.12e-04};
#else
    return {8.82e-03, 2.02e-03, 5.24e-04, 1.42e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.54e-02, 1.47e-02};
#else
    return {6.61e-02, 3.27e-02, 1.81e-02, 1.03e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
