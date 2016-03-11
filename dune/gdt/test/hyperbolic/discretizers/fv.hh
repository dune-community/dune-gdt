// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/periodicview.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/playground/operators/elliptic-cg.hh>
#include <dune/gdt/spaces/fv/defaultproduct.hh>
#include <dune/gdt/spaces/constraints.hh>

#include <dune/gdt/test/hyperbolic/problems/interface.hh>
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {

template< class GridType,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1 >
class FVDiscretizer
{
public:
  typedef ProblemInterface< typename GridType::template Codim< 0 >::Entity,
                            typename GridType::ctype,
                            GridType::dimension,
                            RangeFieldType,
                            dimRange,
                            dimRangeCols> ProblemType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::fv;
  typedef typename DSG::PeriodicGridView
      < typename Stuff::Grid::ProviderInterface< GridType >::LevelGridViewType >             GridViewType;
  typedef typename Spaces::FV::DefaultProduct< GridViewType, RangeFieldType, dimRange, dimRangeCols > FVSpaceType;
  typedef Discretizations::NonStationaryDefault< ProblemType, FVSpaceType >             DiscretizationType;

  static std::string static_id()
  {                                                                                // int() needed, otherwise we get a linker error
    return std::string("gdt.hyperbolic.discretization.fv.dim") + DSC::to_string(int(GridType::dimension));
  }

  static DiscretizationType discretize(Stuff::Grid::ProviderInterface< GridType >& grid_provider,
                                       const ProblemType& problem,
                                       const int level = 0,
                                       const std::bitset< GridType::dimension > periodic_directions = std::bitset< GridType::dimension >())
  {
    auto logger = Stuff::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space = std::make_shared< const FVSpaceType >(GridViewType(grid_provider.level_view(level), periodic_directions));
    logger.debug() << "grid has " << space->grid_view().indexSet().size(0) << " elements" << std::endl;
    return std::move(DiscretizationType(problem, space));
  } // ... discretize(...)
}; // class FVDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
