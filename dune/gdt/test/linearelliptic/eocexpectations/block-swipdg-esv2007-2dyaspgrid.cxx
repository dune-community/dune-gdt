// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)

#include <config.h>

#if HAVE_DUNE_ALUGRID

#include "block-swipdg-esv2007-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007DdSubdomainsTestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::block_ipdg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007DdSubdomainsTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::block_ipdg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {4.62e-02, 1.13e-02};
#else
    return {1.13e-02, 2.90e-03, 7.41e-04, 1.88e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.56e-01, 1.25e-01};
#else
    return {1.25e-01, 6.25e-02, 3.14e-02, 1.58e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007DdSubdomainsTestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::block_ipdg,
                                                  2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007DdSubdomainsTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::block_ipdg,
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


} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_DUNE_ALUGRID
