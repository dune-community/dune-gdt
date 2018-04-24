// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#include <config.h>

#if HAVE_DUNE_ALUGRID

#include "cg-ao2013-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.19e-01, 1.05e-01};
#else
    return {7.92e-02, 4.15e-02, 1.19e-02, 2.72e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {7.02e-01, 5.39e-01};
#else
    return {3.51e-01, 3.02e-01, 1.63e-01, 7.51e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.28e-01, 3.53e-01};
#else
    return {2.72e-01, 1.66e-01, 8.90e-02, 4.27e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.00e-01, 6.51e-02};
#else
    return {2.22e-01, 9.90e-02, 2.96e-02, 6.62e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {6.03e-01, 3.69e-01};
#else
    return {7.00e-01, 4.89e-01, 2.69e-01, 1.25e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.13e-01, 3.21e-01};
#else
    return {4.21e-01, 2.50e-01, 1.34e-01, 6.36e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
