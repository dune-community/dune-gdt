// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/gdt/test/hyperbolic/eocexpectations.hh>
#include <dune/gdt/test/hyperbolic/problems/2dboltzmann.hh>
#include <dune/gdt/test/hyperbolic/problems/burgers.hh>
#include <dune/gdt/test/hyperbolic/problems/transport.hh>
#include <dune/gdt/test/hyperbolic/problems/shallowwater.hh>
#include <dune/gdt/test/hyperbolic/problems/sodshocktube.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
using Yasp2 = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;

typedef testing::Types<Dune::GDT::Hyperbolic::Boltzmann2DCheckerboardTestCase<Yasp2, double, 1>,
                       Dune::GDT::Hyperbolic::BurgersTestCase<Yasp1>, Dune::GDT::Hyperbolic::BurgersTestCase<Yasp2>,
                       Dune::GDT::Hyperbolic::ShallowWaterTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::ShockTubeTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::SourceBeamTestCase<Yasp1, double, 5>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp2, double, 1, 1>>
    YaspGridTestCasesAll;

typedef testing::Types<Dune::GDT::Hyperbolic::BurgersTestCase<Yasp1>, Dune::GDT::Hyperbolic::ShockTubeTestCase<Yasp1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp2, double, 1, 1>>
    YaspGridTestCasesPartial;

typedef testing::Types<Dune::GDT::Hyperbolic::SourceBeamTestCase<Yasp1, double, 5>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Yasp1, double, 1, 1>>
    YaspGridTestCasesLinear1D;


namespace Dune {
namespace GDT {
namespace Tests {

extern template class HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShallowWaterTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                NumericalFluxes::godunov_with_reconstruction,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                NumericalFluxes::godunov_with_reconstruction,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2, NumericalFluxes::godunov,
                                                TimeStepperMethods::explicit_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2, NumericalFluxes::godunov,
                                                TimeStepperMethods::dormand_prince>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2, NumericalFluxes::laxfriedrichs,
                                                TimeStepperMethods::explicit_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
