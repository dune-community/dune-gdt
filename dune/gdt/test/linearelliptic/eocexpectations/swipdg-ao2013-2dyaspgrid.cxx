// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#include <config.h>

#include "swipdg-ao2013-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.39e-01, 6.43e-01};
#else
    return {1.05e-01, 5.00e-02, 1.16e-02, 3.22e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {6.19e-01, 1.65e+00};
#else
    return {6.90e-01, 4.81e-01, 2.28e-01, 1.22e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.68e-01, 3.98e+00};
#else
    return {5.09e-01, 3.44e-01, 2.64e-01, 2.20e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {8.91e-02, 5.90e-01};
#else
    return {8.97e-02, 8.75e-03, 1.89e-03, 6.60e-04};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.37e-01, 1.64e+00};
#else
    return {5.24e-01, 1.47e-01, 7.39e-02, 5.43e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.35e-01, 4.60e+00};
#else
    return {2.70e-01, 2.05e-01, 1.77e-01, 1.46e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
