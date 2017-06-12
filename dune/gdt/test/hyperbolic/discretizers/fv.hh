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

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/spaces/fv/product.hh>

#include <dune/gdt/test/hyperbolic/problems/interface.hh>
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {


template <class TestCaseType,
          class GridType,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          NumericalFluxes numerical_flux = NumericalFluxes::godunov,
          TimeStepperMethods time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp,
          TimeStepperMethods rhs_time_stepper_method = TimeStepperMethods::matrix_exponential>
class FvDiscretizer
{
public:
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::fv;
  static const constexpr NumericalFluxes numerical_flux_type = numerical_flux;
  static const constexpr TimeStepperMethods time_stepper_type = time_stepper_method;
  static const constexpr TimeStepperMethods rhs_time_stepper_type = rhs_time_stepper_method;

  typedef typename XT::Grid::PeriodicGridView<typename XT::Grid::GridProvider<GridType>::LevelGridViewType>
      GridLayerImp;
  typedef Dune::
      GridView<XT::Grid::internal::PeriodicGridViewTraits<typename XT::Grid::GridProvider<GridType>::LevelGridViewType,
                                                          false>>
          GridLayerType;
  typedef FvProductSpace<GridLayerType, RangeFieldType, dimRange, dimRangeCols> FVSpaceType;
  typedef HyperbolicFVDefaultDiscretization<TestCaseType,
                                            FVSpaceType,
                                            numerical_flux,
                                            time_stepper_method,
                                            rhs_time_stepper_method>
      DiscretizationType;

  typedef Hyperbolic::ProblemInterface<typename GridType::template Codim<0>::Entity,
                                       typename GridType::ctype,
                                       GridType::dimension,
                                       typename DiscretizationType::DiscreteFunctionType,
                                       RangeFieldType,
                                       dimRange>
      ProblemType;

  static std::string static_id()
  { // int() needed, otherwise we get a linker error
    return std::string("gdt.hyperbolic.discretization.fv.dim") + Dune::XT::Common::to_string(int(GridType::dimension));
  }

  static DiscretizationType
  discretize(XT::Grid::GridProvider<GridType>& grid_provider,
             const TestCaseType& test_case,
             const int level = 0,
             const std::bitset<GridType::dimension> periodic_directions = std::bitset<GridType::dimension>())
  {
    auto logger = XT::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    GridLayerImp imp(grid_provider.level_view(level), periodic_directions);
    auto space = std::make_shared<const FVSpaceType>(GridLayerType(imp));
    logger.debug() << "grid has " << space->grid_layer().indexSet().size(0) << " elements" << std::endl;
    return DiscretizationType(test_case, space);
  } // ... discretize(...)
}; // class HyperbolicFvDiscretizer


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
