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

#include "cg-esv2007-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.57e-01, 3.82e-02};
#else
    return {3.82e-02, 9.64e-03, 2.42e-03, 6.04e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.77e-01, 1.84e-01};
#else
    return {1.84e-01, 9.24e-02, 4.63e-02, 2.31e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.57e-01, 4.22e-02};
#else
    return {4.22e-02, 1.08e-02, 2.70e-03, 6.76e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.77e-01, 1.94e-01};
#else
    return {1.94e-01, 9.79e-02, 4.91e-02, 2.45e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluCube2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluCube2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.53e-02, 8.28e-03};
#else
    return {8.28e-03, 2.04e-03, 5.09e-04, 1.27e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.33e-01, 1.14e-01};
#else
    return {1.14e-01, 5.68e-02, 2.83e-02, 1.42e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
