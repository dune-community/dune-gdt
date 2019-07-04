// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
#include <dune/gdt/test/mn-discretization.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using Yasp3 = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;

constexpr Dune::GDT::EntropyType entropy = Dune::GDT::EntropyType::BoseEinstein;

using YaspGridTestCasesAll = testing::Types<
#if HAVE_CLP
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7, 1, entropy>, false>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7, 1, entropy>, true>,
    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7, 1, entropy>, false>,
    Dune::GDT::PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7, 1, entropy>, true>,
#endif
    Dune::GDT::
        SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1, 1, entropy>, false>,
    Dune::GDT::SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1, 1, entropy>, true>,
    Dune::GDT::
        PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3, entropy>, false>,
    Dune::GDT::
        PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3, entropy>, true>>;

TYPED_TEST_CASE(HyperbolicMnTest, YaspGridTestCasesAll);
TYPED_TEST(HyperbolicMnTest, check)
{
  this->run();
}
