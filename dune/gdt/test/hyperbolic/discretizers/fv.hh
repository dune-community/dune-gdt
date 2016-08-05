// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/periodicview.hh>

#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/spaces/fv/product.hh>

#include <dune/gdt/test/hyperbolic/problems/interface.hh>
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {


template <class TestCaseType, class GridType, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1,
          NumericalFluxes numerical_flux                                                                 = NumericalFluxes::godunov,
          TimeStepperMethods time_stepper_method                                                         = TimeStepperMethods::explicit_euler>
class FvDiscretizer
{
public:
  typedef Hyperbolic::ProblemInterface<typename GridType::template Codim<0>::Entity, typename GridType::ctype,
                                       GridType::dimension, RangeFieldType, dimRange, dimRangeCols> ProblemType;
      ProblemType;
  static const constexpr ChooseDiscretizer type               = ChooseDiscretizer::fv;
  static const constexpr NumericalFluxes numerical_flux_type  = numerical_flux;
  static const constexpr TimeStepperMethods time_stepper_type = time_stepper_method;

  typedef typename DSG::PeriodicGridView<typename Stuff::Grid::ProviderInterface<GridType>::LevelGridViewType>
      GridViewType;
  typedef FvProductSpace<GridViewType, RangeFieldType, dimRange, dimRangeCols> FVSpaceType;
  typedef HyperbolicFVDefaultDiscretization<TestCaseType, FVSpaceType, numerical_flux, time_stepper_method,
                                            time_stepper_method> DiscretizationType;
      DiscretizationType;

  static std::string static_id()
  { // int() needed, otherwise we get a linker error
    return std::string("gdt.hyperbolic.discretization.fv.dim") + DSC::to_string(int(GridType::dimension));
  }

  static DiscretizationType
  discretize(Stuff::Grid::ProviderInterface<GridType>& grid_provider, const TestCaseType& test_case,
             const int level                                            = 0,
             const std::bitset<GridType::dimension> periodic_directions = std::bitset<GridType::dimension>())
  {
    auto logger = Stuff::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space =
        std::make_shared<const FVSpaceType>(GridViewType(grid_provider.level_view(level), periodic_directions));
    logger.debug() << "grid has " << space->grid_view().indexSet().size(0) << " elements" << std::endl;
    return DiscretizationType(test_case, space);
  } // ... discretize(...)
}; // class HyperbolicFvDiscretizer


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
