// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH

#include <dune/grid/sgrid.hh>
#include <dune/alugrid/dgf.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include "swipdg-estimator-expectations.hh"
#include "problems/AO2013.hh"
#include "problems/ESV2007.hh"
#include "problems/spe10.hh"

#if HAVE_DUNE_ALUGRID

typedef testing::Types<Dune::GDT::LinearElliptic::AO2013TestCase<Dune::GDT::Test::Alu2NonConformSimplex>,
                       Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::GDT::Test::Alu2NonConformSimplex>,
                       Dune::GDT::LinearElliptic::
                           Spe10Model1TestCase<Dune::GDT::Test::Alu2NonConformSimplex>> AluGridTestCases;


#endif // HAVE_DUNE_ALUGRID


namespace Dune {
namespace GDT {
namespace Test {

#if HAVE_DUNE_ALUGRID


extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::AO2013TestCase<Alu2NonConformSimplex,
                                                                                               double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<Alu2NonConformSimplex,
                                                                                                double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<Alu2NonConformSimplex,
                                                                                                    double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;


#endif // HAVE_DUNE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
