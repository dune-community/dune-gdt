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

#include "swipdg-spe10-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.06e+00, 7.58e-01};
#else
    return {2.73e-02, 2.35e-02};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.03e+00, 1.03e+00};
#else
    return {4.64e-01, 6.31e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.39e+00, 2.84e+00};
#else
    return {1.03e+00, 1.49e+00};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.06e+00, 7.47e-01};
#else
    return {1.82e-02, 7.50e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.05e+00, 1.05e+00};
#else
    return {5.53e-01, 4.55e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.29e+00, 1.67e+00};
#else
    return {1.62e+00, 1.33e+00};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
