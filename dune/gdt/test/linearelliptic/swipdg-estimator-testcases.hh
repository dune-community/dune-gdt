// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH

#include <dune/grid/sgrid.hh>
#include <dune/grid/alugrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include "swipdg-estimator-expectations.hh"
#include "problems/AO2013.hh"
#include "problems/ESV2007.hh"
#include "problems/spe10.hh"

#if HAVE_ALUGRID


typedef testing::
    Types<Dune::GDT::LinearElliptic::AO2013TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::Spe10Model1TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>
        AluGridTestCases;


#endif // HAVE_ALUGRID


namespace Dune {
namespace GDT {
namespace Test {

#if HAVE_ALUGRID


extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                                       conforming>,
                                                                                               double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                                        conforming>,
                                                                                                double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2,
                                                                                                            simplex,
                                                                                                            conforming>,
                                                                                                    double, 1>,
                                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>;


#endif // HAVE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_ESTIMATOR_TESTCASES_HH
