// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/testcases.hh>
#include <dune/gdt/test/hyperbolic/mn-discretization.hh>

using Yasp1 = Yasp1Grid;
using Yasp2 = Yasp2Grid;
using Yasp3 = Yasp3Grid;

using YaspGridTestCasesAll = testing::
    Types<Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, true>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::LegendreMomentBasis<double, double, 7>, true>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, true>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              SourceBeamMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PlaneSourceMnTestCase<Yasp1, Dune::GDT::PartialMomentBasis<double, 1, double, 8, 1, 1>, true>
#if HAVE_CLP
          ,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, true>,
#endif
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, true>
#if HAVE_QHULL
          ,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, false>,
          Dune::GDT::Hyperbolic::Problems::KineticTransport::
              PointSourceMnTestCase<Yasp3, Dune::GDT::PartialMomentBasis<double, 3, double, 0, 1, 3>, true>
#endif
          >;

TYPED_TEST_CASE(HyperbolicMnTest, YaspGridTestCasesAll);
TYPED_TEST(HyperbolicMnTest, check)
{
  this->run();
}
