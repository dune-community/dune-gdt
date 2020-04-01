// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#define USE_LP_POSITIVITY_LIMITER 0

#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
#include <dune/gdt/test/momentmodels/mn-discretization.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

using YaspGridTestCasesOrd1 = testing::Types<
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 1>, false, false>,
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 10>, false, false>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 20>, false, false>
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 30>, false, false>,
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 40>, false, false>,
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 50>, false, false>,
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 2, 1, 1>, false,
    // false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 10, 1, 1>,
    // false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 20,
    // 1, 1>, false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1,
    // double, 30, 1, 1>, false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1,
    // Dune::GDT::HatFunctionMomentBasis<double, 1, double, 40, 1, 1>, false, false>,
    // Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 50, 1, 1>, false,
    // false>, Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 50, 1, 1>, false,
    // false> Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 10, 1, 1>, false,
    // false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 20, 1, 1>,
    // false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 30, 1,
    // 1>, false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 40,
    // 1, 1>, false, false>, Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double,
    // 50, 1, 1>, false, false>
    >;

TYPED_TEST_CASE(HyperbolicMnTest, YaspGridTestCasesOrd1);
TYPED_TEST(HyperbolicMnTest, check)
{
  this->run();
}
