// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_MOMENTMODELS_PN_HH
#define DUNE_GDT_TEST_HYPERBOLIC_MOMENTMODELS_PN_HH

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

using YaspGridTestCasesWithoutReconstruction = testing::Types<
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false>,
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, false>,
    Dune::GDT::CheckerboardPnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, false>,
    Dune::GDT::ShadowPnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, false>
#if !DXT_DISABLE_LARGE_TESTS
    ,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 1, 1, 3>, false>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, false>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 1, 1, 3>, false>
#endif
    >;

using YaspGridTestCasesWithReconstruction = testing::Types<
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, true>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, true>,
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true>,
    Dune::GDT::SourceBeamPnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>,
    Dune::GDT::PlaneSourcePnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>
#if !DXT_DISABLE_LARGE_TESTS
    ,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, true>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, true>,
    Dune::GDT::PointSourcePnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, true>
#endif
    >;

#endif // DUNE_GDT_TEST_HYPERBOLIC_MOMENTMODELS_PN_HH
