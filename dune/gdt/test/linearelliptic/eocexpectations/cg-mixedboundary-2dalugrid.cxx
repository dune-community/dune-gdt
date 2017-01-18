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

#if HAVE_ALUGRID

#include "cg-mixedboundary-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg,
                              1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.83e-02, 9.80e-03};
#else
    return {8.31e-02, 2.22e-02, 5.52e-03, 1.19e-03};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.00e-01, 1.08e-01};
#else
    return {3.11e-01, 1.64e-01, 8.23e-02, 3.75e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg,
                              1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.82e-02, 8.39e-03};
#else
    return {1.35e-01, 3.99e-02, 1.02e-02, 2.16e-03};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.02e-01, 9.45e-02};
#else
    return {3.75e-01, 2.08e-01, 1.06e-01, 4.87e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_ALUGRID
