// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/test/hyperbolic/eocexpectations.hh>
#include <dune/gdt/test/hyperbolic/problems/2dboltzmann.hh>
#include <dune/gdt/test/hyperbolic/problems/burgers.hh>
#include <dune/gdt/test/hyperbolic/problems/transport.hh>
#include <dune/gdt/test/hyperbolic/problems/shallowwater.hh>
#include <dune/gdt/test/hyperbolic/problems/sodshocktube.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>


typedef testing::Types<Dune::GDT::Hyperbolic::Boltzmann2DCheckerboardTestCase<Dune::YaspGrid<2>, double, 1>,
                       Dune::GDT::Hyperbolic::BurgersTestCase<Dune::YaspGrid<1>>,
                       Dune::GDT::Hyperbolic::BurgersTestCase<Dune::YaspGrid<2>>,
                       Dune::GDT::Hyperbolic::ShallowWaterTestCase<Dune::YaspGrid<1>>,
                       Dune::GDT::Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>>,
                       Dune::GDT::Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double, 5>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1, 1>> YaspGridTestCasesAll;

typedef testing::Types<Dune::GDT::Hyperbolic::BurgersTestCase<Dune::YaspGrid<1>>,
                       Dune::GDT::Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1, 1>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1, 1>>
    YaspGridTestCasesPartial;

typedef testing::Types<Dune::GDT::Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double, 5>,
                       Dune::GDT::Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1, 1>>
    YaspGridTestCasesLinear1D;


namespace Dune {
namespace GDT {
namespace Tests {

extern template class HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Dune::YaspGrid<2>, double,
                                                                                            1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Dune::YaspGrid<2>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShallowWaterTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

extern template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Dune::YaspGrid<1>, double>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::
                                                    godunovwithreconstruction_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<1>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                                Hyperbolic::FluxTimeStepperCombinations::
                                                    godunovwithreconstruction_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

extern template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
