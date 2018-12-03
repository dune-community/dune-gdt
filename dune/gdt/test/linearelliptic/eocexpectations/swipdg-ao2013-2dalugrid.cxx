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

#include "swipdg-ao2013-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1>::results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.05e-01, 1.36e-01};
#else
    return {5.57e-02, 1.99e-02, 5.54e-03, 1.29e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {6.16e-01, 4.65e-01};
#else
    return {4.32e-01, 2.93e-01, 1.50e-01, 6.54e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.58e-01, 8.08e-01};
#else
    return {2.53e-01, 1.46e-01, 7.45e-02, 3.35e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    2>::results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.49e-01, 1.38e-01};
#else
    return {1.18e-02, 2.11e-03, 3.89e-04, 7.76e-05};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.72e-01, 7.52e-01};
#else
    return {1.69e-01, 5.96e-02, 1.94e-02, 6.04e-03};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.47e-01, 1.55e+00};
#else
    return {7.41e-02, 3.36e-02, 1.62e-02, 7.03e-03};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, nonconforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1>::results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.91e-01, 1.20e-01};
#else
    return {1.05e-01, 3.90e-02, 1.27e-02, 3.13e-03};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {5.37e-01, 3.70e-01};
#else
    return {6.52e-01, 4.24e-01, 2.20e-01, 9.59e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.47e-01, 6.11e-01};
#else
    return {3.67e-01, 2.08e-01, 1.06e-01, 4.75e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, nonconforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    2>::results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.23e-01, 1.33e-01};
#else
    return {2.77e-02, 5.53e-03, 8.36e-04, 1.29e-04};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.91e-01, 7.52e-01};
#else
    return {2.75e-01, 1.17e-01, 3.67e-02, 1.10e-02};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.31e-01, 1.38e+00};
#else
    return {1.15e-01, 5.29e-02, 2.34e-02, 1.05e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
