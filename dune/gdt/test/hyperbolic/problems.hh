// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/test/hyperbolic/eocexpectations.hh>
#include <dune/gdt/test/hyperbolic/problems/burgers.hh>
#include <dune/gdt/test/hyperbolic/problems/transport.hh>
#include <dune/gdt/test/hyperbolic/problems/shallowwater.hh>
#include <dune/gdt/test/hyperbolic/problems/sodshocktube.hh>


typedef testing::Types<    Dune::GDT::Hyperbolic::BurgersTestCase
                                                < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > > >
                         , Dune::GDT::Hyperbolic::BurgersTestCase
                                                < Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > > >
                         , Dune::GDT::Hyperbolic::TransportTestCase
                                                < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >, double, 1, 1 >
                         , Dune::GDT::Hyperbolic::TransportTestCase
                                                < Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >, double, 1, 1 >
                         , Dune::GDT::Hyperbolic::ShallowWaterTestCase
                                                < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > > >
                         , Dune::GDT::Hyperbolic::ShockTubeTestCase
                                                < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > > >
                      > YaspGridTestCases;


namespace Dune {
namespace GDT {
namespace Tests {


extern template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase
                                                 < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >,
                                                   double,
                                                   1 >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 1 >;

extern template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase
                                                 < Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >,
                                                   double,
                                                   1 >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 2 >;

extern template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase
                                                 < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >,
                                                   double,
                                                   1,
                                                   1 >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 1 >;

extern template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase
                                                 < Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >,
                                                   double,
                                                   1,
                                                   1 >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 2 >;

extern template class HyperbolicEocExpectations< Hyperbolic::ShallowWaterTestCase
                                                 < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >,
                                                   double >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 1 >;

extern template class HyperbolicEocExpectations< Hyperbolic::ShockTubeTestCase
                                                 < Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >,
                                                   double >,
                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                 1 >;

} // namespace Tests
} // namespace GDT
} // namespace Dune
