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

#include "cg-spe10-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {7.72e-01, 8.92e-01};
#else
    return {3.44e-02, 1.01e-02};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {9.27e-01, 9.94e-01};
#else
    return {1.47e-01, 7.78e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.50e+00, 2.60e+00};
#else
    return {1.88e-01, 1.00e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {8.04e-01, 6.93e-01};
#else
    return {5.08e-02, 1.54e-02};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {9.40e-01, 8.83e-01};
#else
    return {2.01e-01, 1.04e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.51e+00, 1.53e+00};
#else
    return {2.30e-01, 1.25e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
} // ... results(...)

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluCube2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluCube2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.04e+00, 7.25e-01};
#else
    return {1.86e-02, 1.51e-02};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.02e+00, 1.01e+00};
#else
    return {3.31e-01, 4.32e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.43e+00, 3.24e+00};
#else
    return {9.58e-01, 1.37e+00};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_DUNE_ALUGRID
