// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/grid/sgrid.hh>
#include <dune/grid/alugrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/tests/linearelliptic/eocexpectations.hh>
#include <dune/gdt/tests/linearelliptic/problems/AO2013.hh>
#include <dune/gdt/tests/linearelliptic/problems/ER2007.hh>
#include <dune/gdt/tests/linearelliptic/problems/ESV2007.hh>
#include <dune/gdt/tests/linearelliptic/problems/mixedboundary.hh>
#include <dune/gdt/tests/linearelliptic/problems/spe10.hh>


typedef testing::Types< Dune::GDT::LinearElliptic::AO2013TestCase< Dune::YaspGrid< 2 > >
                      , Dune::GDT::LinearElliptic::ER2007TestCase< Dune::YaspGrid< 2 > >
                      , Dune::GDT::LinearElliptic::ESV2007TestCase< Dune::YaspGrid< 2 > >
                      , Dune::GDT::LinearElliptic::MixedBoundaryTestCase< Dune::YaspGrid< 2 > >
                      , Dune::GDT::LinearElliptic::Spe10Model1TestCase< Dune::YaspGrid< 2 > >
                      > SGridTestCases;


#if HAVE_ALUGRID


typedef testing::Types< Dune::GDT::LinearElliptic::AO2013TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >
                      , Dune::GDT::LinearElliptic::AO2013TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > >
                      , Dune::GDT::LinearElliptic::ER2007TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >
                      , Dune::GDT::LinearElliptic::ER2007TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > >
                      , Dune::GDT::LinearElliptic::ESV2007TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >
                      , Dune::GDT::LinearElliptic::ESV2007TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > >
                      , Dune::GDT::LinearElliptic::MixedBoundaryTestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >
                      , Dune::GDT::LinearElliptic::MixedBoundaryTestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > >
                      , Dune::GDT::LinearElliptic::Spe10Model1TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >
                      , Dune::GDT::LinearElliptic::Spe10Model1TestCase< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > >
                      > AluGridTestCases;


#endif // HAVE_ALUGRID


namespace Dune {
namespace GDT {
namespace Tests {


extern template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< YaspGrid< 2 >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::ER2007TestCase< YaspGrid< 2 >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< YaspGrid< 2 >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::MixedBoundaryTestCase< YaspGrid< 2 >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< YaspGrid< 2 >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;



#if HAVE_ALUGRID


extern template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;
extern template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::ER2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;
extern template class LinearEllipticEocExpectations< LinearElliptic::ER2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;
extern template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::MixedBoundaryTestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;
extern template class LinearEllipticEocExpectations< LinearElliptic::MixedBoundaryTestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;

extern template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;
extern template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                                     LinearElliptic::ChooseDiscretizer::cg,
                                                     1 >;



#endif // HAVE_ALUGRID

} // namespace Tests
} // namespace GDT
} // namespace Dune
