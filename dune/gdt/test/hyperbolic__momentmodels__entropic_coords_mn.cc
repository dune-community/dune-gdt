// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#define USE_LP_POSITIVITY_LIMITER 1

#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
#include <dune/gdt/test/entropic-coords-mn-discretization.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

using YaspGridTestCasesAll = testing::Types<
#if HAVE_CLP
//    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
//    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 20>, true>
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 20>, true>
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 21>, true>
#endif
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false, true>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true, true>
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 20, 1, 1>, false>,
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 100, 1, 1>, true>
//    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>
//    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>,
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>
//    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>
#if !DXT_DISABLE_LARGE_TESTS
//    ,
#  if HAVE_CLP
//    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>,
//    false>, Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2,
//    3>, true> Dune::GDT::CheckerboardMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>,
//    true> Dune::GDT::CheckerboardMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>,
//    false>, Dune::GDT::ShadowMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>,
//    false>, Dune::GDT::ShadowMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 1, 1, 3>, true>
#  endif
//    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false>
// Our shifted qr eigensolver fails for this problem, needs better shifting strategy
#  if HAVE_MKL || HAVE_LAPACKE || HAVE_EIGEN
//    ,
//    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, true>
#  endif
#  if HAVE_QHULL
//    ,
//    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, false>,
//    Dune::GDT::PointSourceMnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, true>
#  endif
#endif
    >;

TYPED_TEST_CASE(HyperbolicEntropicCoordsMnTest, YaspGridTestCasesAll);
TYPED_TEST(HyperbolicEntropicCoordsMnTest, check)
{
  this->run();
}
