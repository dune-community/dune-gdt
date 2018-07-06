// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#include <config.h>

#if HAVE_DUNE_ALUGRID

#include "swipdg-esv2007-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {8.33e-02, 1.82e-02};
#else
    return {1.82e-02, 4.53e-03, 1.12e-03, 2.78e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.04e-01, 1.48e-01};
#else
    return {1.48e-01, 7.28e-02, 3.62e-02, 1.80e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.77e-03, 8.55e-04};
#else
    return {8.55e-04, 1.06e-04, 1.31e-05, 1.63e-06};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.20e-02, 1.41e-02};
#else
    return {1.41e-02, 3.56e-03, 8.91e-04, 2.23e-04};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, nonconforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)

{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {8.33e-02, 2.31e-02};
#else
    return {2.31e-02, 5.96e-03, 1.50e-03, 3.76e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.04e-01, 1.49e-01};
#else
    return {1.49e-01, 7.33e-02, 3.63e-02, 1.81e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, nonconforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.77e-03, 6.78e-04};
#else
    return {6.78e-04, 8.23e-05, 1.02e-05, 1.26e-06};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.20e-02, 1.32e-02};
#else
    return {1.32e-02, 3.30e-03, 8.25e-04, 2.06e-04};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_DUNE_ALUGRID
