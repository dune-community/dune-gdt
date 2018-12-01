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

#  include "swipdg-er2007-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1>::results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.06e+00, 1.95e-01};
#  else
    return {1.11e-01, 6.09e-02, 1.65e-02};
#  endif
  } else if (type == "H1_semi" || type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {9.37e-01, 4.39e-01};
#  else
    return {3.66e-01, 2.98e-01, 1.46e-01};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    2>::results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.95e-01, 1.51e-01};
#  else
    return {6.35e-02, 6.42e-03, 8.23e-04};
#  endif
  } else if (type == "H1_semi" || type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.04e-01, 3.46e-01};
#  else
    return {2.32e-01, 5.40e-02, 1.41e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, nonconforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1>::results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.06e+00, 2.00e-01};
#  else
    return {2.00e-01, 1.11e-01, 3.29e-02};
#  endif
  } else if (type == "H1_semi" || type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {9.37e-01, 4.41e-01};
#  else
    return {4.41e-01, 4.05e-01, 2.04e-01};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, nonconforming

std::vector<double> LinearEllipticEocExpectations<
    LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    2>::results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>::TestCaseType&,
                const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.95e-01, 1.51e-01};
#  else
    return {1.51e-01, 1.27e-02, 1.50e-03};
#  endif
  } else if (type == "H1_semi" || type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.04e-01, 3.69e-01};
#  else
    return {3.69e-01, 8.28e-02, 2.12e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
