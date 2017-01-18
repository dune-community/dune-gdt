// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/alugrid.hh>

#include <dune/xt/common/test/gtest/gtest.h>

#include "eocexpectations.hh"
#include "problems.hh"


typedef testing::Types<Dune::GDT::LinearElliptic::AO2013TestCase<Yasp2Grid>,
                       Dune::GDT::LinearElliptic::ER2007TestCase<Yasp2Grid>,
                       Dune::GDT::LinearElliptic::ESV2007TestCase<Yasp2Grid>,
                       Dune::GDT::LinearElliptic::MixedBoundaryTestCase<Yasp2Grid>
#if !DXT_DISABLE_LARGE_TESTS
                       ,
                       Dune::GDT::LinearElliptic::Spe10Model1TestCase<Yasp2Grid>
#endif
                       >
    YaspGridTestCases;


#if HAVE_ALUGRID


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


#endif // HAVE_ALUGRID


namespace Dune {
namespace GDT {
namespace Test {


// YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      MixedBoundaryTestCase<Dune::YaspGrid<2,
                                                                           Dune::EquidistantOffsetCoordinates<double,
                                                                                                              2>>,
                                                            double,
                                                            1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                  1>;

extern template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      Spe10Model1TestCase<Dune::YaspGrid<2,
                                                                         Dune::EquidistantOffsetCoordinates<double, 2>>,
                                                          double,
                                                          1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                  1>;

// YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      MixedBoundaryTestCase<Dune::YaspGrid<2,
                                                                           Dune::EquidistantOffsetCoordinates<double,
                                                                                                              2>>,
                                                            double,
                                                            1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                  2>;

extern template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      Spe10Model1TestCase<Dune::YaspGrid<2,
                                                                         Dune::EquidistantOffsetCoordinates<double, 2>>,
                                                          double,
                                                          1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg,
                                  2>;


#if HAVE_ALUGRID


// ALUGrid< 2, 2, simplex, conforming >, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        MixedBoundaryTestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                              double,
                                                                              1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                            double,
                                                                            1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;


// ALUGrid< 2, 2, simplex, nonconforming >, polorder 1

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                       double,
                                                                       1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                       double,
                                                                       1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                        double,
                                                                        1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        MixedBoundaryTestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                              double,
                                                                              1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double,
                                                                            1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    1>;


// ALUGrid< 2, 2, simplex, conforming >, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        MixedBoundaryTestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                              double,
                                                                              1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                            double,
                                                                            1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;


// ALUGrid< 2, 2, simplex, nonconforming >, polorder 2

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                       double,
                                                                       1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                       double,
                                                                       1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                        double,
                                                                        1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        MixedBoundaryTestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                              double,
                                                                              1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::
                                                        Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double,
                                                                            1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                                    2>;


#endif // HAVE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_SWIPDG_TESTCASES_HH
