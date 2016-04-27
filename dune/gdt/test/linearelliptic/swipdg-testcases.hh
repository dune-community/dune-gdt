// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH

#include <dune/grid/sgrid.hh>
#include <dune/alugrid/dgf.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include "eocexpectations.hh"
#include "problems.hh"


typedef testing::Types<Dune::GDT::LinearElliptic::AO2013TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::ER2007TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::MixedBoundaryTestCase<Dune::SGrid<2, 2>>
#if !THIS_IS_A_BUILDBOT_BUILD
                       ,
                       Dune::GDT::LinearElliptic::Spe10Model1TestCase<Dune::SGrid<2, 2>>
#endif
                       > SGridTestCases;


#if HAVE_DUNE_ALUGRID


typedef testing::
    Types<Dune::GDT::LinearElliptic::AO2013TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::AO2013TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>>,
          Dune::GDT::LinearElliptic::ER2007TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::ER2007TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>>,
          Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>>,
          Dune::GDT::LinearElliptic::MixedBoundaryTestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::MixedBoundaryTestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>>,
          Dune::GDT::LinearElliptic::Spe10Model1TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>,
          Dune::GDT::LinearElliptic::Spe10Model1TestCase<Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>>>
        AluGridTestCases;


#endif // HAVE_DUNE_ALUGRID


namespace Dune {
namespace GDT {
namespace Test {


// SGrid< 2, 2 >, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

// SGrid< 2, 2 >, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;


#if HAVE_DUNE_ALUGRID


// ALUGrid< 2, 2, simplex, conforming >, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  conforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                conforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;


// ALUGrid< 2, 2, simplex, nonconforming >, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                            nonconforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  nonconforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                nonconforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;


// ALUGrid< 2, 2, simplex, conforming >, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  conforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                conforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;


// ALUGrid< 2, 2, simplex, nonconforming >, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                            nonconforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  nonconforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                nonconforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;


#endif // HAVE_DUNE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH
