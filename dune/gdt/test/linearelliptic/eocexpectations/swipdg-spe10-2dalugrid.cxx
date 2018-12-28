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

#  include "swipdg-spe10-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.03e+00, 1.58e+00};
#  else
    return {9.48e-03, 2.64e-03};
#  endif
  } else if (type == "H1_semi") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.20e+00, 1.65e+00};
#  else
    return {1.09e-01, 5.36e-02};
#  endif
  } else if (type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {6.88e+00, 8.85e+00};
#  else
    return {1.37e-01, 6.60e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {9.42e-01, 1.44e+00};
#  else
    return {1.05e+00, 1.05e+00};
#  endif
  } else if (type == "H1_semi") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.26e+00, 1.83e+00};
#  else
    return {4.91e-02, 2.26e-02};
#  endif
  } else if (type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.11e+00, 1.07e+01};
#  else
    return {5.60e-02, 2.58e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, nonconforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.05e+00, 1.04e+00};
#  else
    return {9.43e-03, 2.55e-03};
#  endif
  } else if (type == "H1_semi") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.34e+00, 1.15e+00};
#  else
    return {1.43e-01, 7.04e-02};
#  endif
  } else if (type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {6.94e+00, 4.95e+00};
#  else
    return {1.74e-01, 8.53e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, nonconforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#  if DXT_DISABLE_LARGE_TESTS
    return {9.99e-01, 9.74e-01};
#  else
    return {2.13e-03, 6.56e-04};
#  endif
  } else if (type == "H1_semi") {
#  if DXT_DISABLE_LARGE_TESTS
    return {1.27e+00, 1.53e+00};
#  else
    return {6.41e-02, 3.24e-02};
#  endif
  } else if (type == "energy") {
#  if DXT_DISABLE_LARGE_TESTS
    return {7.03e+00, 8.36e+00};
#  else
    return {7.56e-02, 3.81e-02};
#  endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
