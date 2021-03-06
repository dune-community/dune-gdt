// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#define ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING 1

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#if HAVE_DUNE_XT_DATA

#  include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
#  include <dune/gdt/test/momentmodels/entropic-coords-mn-discretization.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

using YaspGridTestCasesAll = testing::Types<
    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false, true>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false, true>,
    Dune::GDT::
        CheckerboardMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false, true>,
    Dune::GDT::ShadowMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false, true>>;

TYPED_TEST_SUITE(HyperbolicEntropicCoordsMnTest, YaspGridTestCasesAll);
TYPED_TEST(HyperbolicEntropicCoordsMnTest, check)
{
  this->run();
}

#else // HAVE_DUNE_XT_DATA

GTEST_TEST(HyperbolicEntropicCoordsMnTest, YaspGridTestCasesAll)
{
  std::cerr << "Test disabled, missing dune-xt-data!" << std::endl;
}

#endif // HAVE_DUNE_XT_DATA
