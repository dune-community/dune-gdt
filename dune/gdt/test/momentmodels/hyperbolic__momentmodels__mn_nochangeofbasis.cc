// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#define ENTROPY_FLUX_UNSPECIALIZED_USE_ADAPTIVE_CHANGE_OF_BASIS 0

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
#include <dune/gdt/test/momentmodels/mn-discretization.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

using YaspGridTestCasesNoBasisChange = testing::Types<
#if HAVE_CLP
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, true>,
#endif
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true>
#if !DXT_DISABLE_LARGE_TESTS
    ,
#  if HAVE_CLP
    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, false>,
#  endif
    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false>
#endif
    >;

TYPED_TEST_SUITE(HyperbolicMnTest, YaspGridTestCasesNoBasisChange);
TYPED_TEST(HyperbolicMnTest, check)
{
  this->run();
}
