// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/grid/sgrid.hh>
#include <dune/grid/alugrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include "eocexpectations.hh"
#include "problems/AO2013.hh"
#include "problems/ER2007.hh"
#include "problems/ESV2007.hh"
#include "problems/mixedboundary.hh"
#include "problems/spe10.hh"


typedef testing::Types<Dune::GDT::LinearElliptic::AO2013TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::ER2007TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::ESV2007TestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::MixedBoundaryTestCase<Dune::SGrid<2, 2>>,
                       Dune::GDT::LinearElliptic::Spe10Model1TestCase<Dune::SGrid<2, 2>>> SGridTestCases;


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


// CG, polorder 1
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<SGrid<2, 2>, double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

// SWIPDG, polorder 1
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

// SWIPDG, polorder 2
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


#if HAVE_ALUGRID


// CG, polorder 1
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                            nonconforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  conforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  nonconforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                conforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                nonconforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::cg, 1>;

// SWIPDG, polorder 1
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                            nonconforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  conforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  nonconforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                conforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;
extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                nonconforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>;

// SWIPDG, polorder 2
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex,
                                                                                           nonconforming>,
                                                                                   double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;
extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                            nonconforming>,
                                                                                    double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  conforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;
extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<ALUGrid<2, 2, simplex,
                                                                                                  nonconforming>,
                                                                                          double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                conforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;
extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                nonconforming>,
                                                                                        double, 1>,
                                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>;


#endif // HAVE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune
