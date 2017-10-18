// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#include <dune/gdt/test/hyperbolic/eocexpectations.hh>
#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/hyperbolic/problems/burgers.hh>
#include <dune/gdt/test/hyperbolic/problems/transport.hh>
#include <dune/gdt/test/hyperbolic/problems/shallowwater.hh>
#include <dune/gdt/test/hyperbolic/problems/sodshocktube.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/checkerboard.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/fokkerplanck/sourcebeam.hh>

using Yasp1 = Yasp1Grid;
using Yasp2 = Yasp2Grid;
using Yasp3 = Yasp3Grid;

typedef testing::Types<Dune::GDT::Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                       Dune::GDT::Hyperbolic::CheckerboardTestCase<Yasp2, double, 1>,
                       Dune::GDT::Hyperbolic::BurgersTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::BurgersTestCase<Yasp2>,
                       Dune::GDT::Hyperbolic::ShallowWaterTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::ShockTubeTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp2, double, 1>>
    YaspGridTestCasesAll;

typedef testing::Types<Dune::GDT::Hyperbolic::BurgersTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::ShockTubeTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp2, double, 1>>
    YaspGridTestCasesPartial;

typedef testing::Types<Dune::GDT::Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1>>
    YaspGridTestCasesLinear1D;


namespace Dune {
namespace GDT {
namespace Test {


extern template class HyperbolicEocExpectations<Hyperbolic::CheckerboardTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                2,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                2,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShallowWaterTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler,
                                                0>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler,
                                                1>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler,
                                                0>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                1,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler,
                                                1>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                2,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                2,
                                                NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv,
                                                2,
                                                NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler,
                                                TimeStepperMethods::explicit_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
