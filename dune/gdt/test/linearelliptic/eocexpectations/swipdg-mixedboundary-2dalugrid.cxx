// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#include <config.h>

#if HAVE_ALUGRID

#include "swipdg-mixedboundary-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg,
                              1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.11e-02, 5.11e-03};
#else
    return {4.02e-02, 1.12e-02, 2.83e-03, 6.33e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.68e-01, 8.65e-02};
#else
    return {2.69e-01, 1.39e-01, 6.87e-02, 3.08e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg,
                              2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.53e-03, 3.10e-04};
#else
    return {3.58e-03, 6.25e-04, 1.21e-04, 2.68e-05};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.12e-02, 1.04e-02};
#else
    return {4.81e-02, 1.79e-02, 7.19e-03, 2.85e-03};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, noncoforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg,
                              1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.09e-02, 4.85e-03};
#else
    return {6.84e-02, 2.17e-02, 5.78e-03, 1.27e-03};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.70e-01, 7.64e-02};
#else
    return {3.28e-01, 1.76e-01, 8.72e-02, 3.89e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, noncoforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg,
                              2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.53e-03, 2.68e-04};
#else
    return {9.67e-03, 1.52e-03, 2.75e-04, 5.63e-05};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.14e-02, 1.09e-02};
#else
    return {9.33e-02, 3.19e-02, 1.19e-02, 4.66e-03};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
