// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/periodicview.hh>

#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/spaces/fv/defaultproduct.hh>

#include <dune/gdt/test/hyperbolic/problems/interface.hh>
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {

template <class TestCaseType, class GridType, class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1,
          bool use_lax_friedrichs_flux = false, bool use_adaptive_timestepper = false,
          bool use_linear_reconstruction = false>
class FVDiscretizer
{
public:
  typedef ProblemInterface<typename GridType::template Codim<0>::Entity, typename GridType::ctype, GridType::dimension,
                           RangeFieldType, dimRange, dimRangeCols> ProblemType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::fv;
  static const constexpr FluxTimeStepperCombinations flux_and_timestepper_type =
      use_lax_friedrichs_flux
          ? FluxTimeStepperCombinations::laxfriedrichs_euler
          : (use_linear_reconstruction ? FluxTimeStepperCombinations::godunovwithreconstruction_euler
                                       : (use_adaptive_timestepper ? FluxTimeStepperCombinations::godunov_adaptiveRK
                                                                   : FluxTimeStepperCombinations::godunov_euler));
  typedef
      typename DSG::PeriodicGridView<typename Stuff::Grid::ProviderInterface<GridType>::LevelGridViewType> GridViewType;
  typedef typename Spaces::FV::DefaultProduct<GridViewType, RangeFieldType, dimRange, dimRangeCols> FVSpaceType;
  typedef Discretizations::NonStationaryDefault<TestCaseType, FVSpaceType, use_lax_friedrichs_flux,
                                                use_adaptive_timestepper, use_linear_reconstruction> DiscretizationType;

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
}; // class FVDiscretizer


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
