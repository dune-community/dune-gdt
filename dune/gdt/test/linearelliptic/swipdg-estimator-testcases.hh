// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/xt/common/test/gtest/gtest.h>

#include "swipdg-estimator-expectations.hh"
#include "problems/AO2013.hh"
#include "problems/ESV2007.hh"
#include "problems/spe10.hh"

#if HAVE_DUNE_ALUGRID


typedef testing::Types<Dune::GDT::LinearElliptic::AO2013TestCase<AluConform2dGridType>,
                       Dune::GDT::LinearElliptic::ESV2007TestCase<AluConform2dGridType>,
                       Dune::GDT::LinearElliptic::Spe10Model1TestCase<AluConform2dGridType>>
    AluGridTestCases;


#endif // HAVE_DUNE_ALUGRID


namespace Dune {
namespace GDT {
namespace Test {

#if HAVE_DUNE_ALUGRID


extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    AO2013TestCase<AluConform2dGridType, double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;
extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;
extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    Spe10Model1TestCase<AluConform2dGridType,
                                                                                        double,
                                                                                        1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;
extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                                    Spe10Model1TestCase<AluSimplex2dGridType,
                                                                                        double,
                                                                                        1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                                1>;


#endif // HAVE_DUNE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
