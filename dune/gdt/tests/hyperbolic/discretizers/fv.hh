// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH

#include <config.h>

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
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/operators/elliptic-cg.hh>
#include <dune/gdt/spaces/fv/defaultproduct.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "../problems/interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {

template< class GridType, bool periodic = false >
struct ChooseGridView {
  typedef typename Stuff::Grid::ProviderInterface< GridType >::LevelGridViewType type;
};

template< class GridType >
struct ChooseGridView< GridType, true > {
  typedef typename DSG::PeriodicGridView
      < typename Stuff::Grid::ProviderInterface< GridType >::LevelGridViewType > type;
};

/**
 * \brief Discretizes a linear elliptic PDE using a continuous Galerkin Finite Element method.
 * \tparam GG GridType
 * \tparam ll layer
 * \tparam ss space_backend
 * \tparam la la_backend
 * \tparam pp polOrder
 * \tparam RR RangeFieldType
 * \tparam rr dimRange
 */
template< class GridType,
          class RangeFieldType = double,
          size_t dimRange = 1,
          bool use_periodic_grid_view = false >
class FVDiscretizer
{
public:
  typedef ProblemInterface< typename GridType::template Codim< 0 >::Entity,
                            typename GridType::ctype,
                            GridType::dimension,
                            RangeFieldType,
                            dimRange > ProblemType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::fv;
  typedef typename ChooseGridView< GridType, use_periodic_grid_view >::type             GridViewType;
  typedef typename Spaces::FV::DefaultProduct< GridViewType, RangeFieldType, dimRange > FVSpaceType;
  typedef Discretizations::NonStationaryDefault< ProblemType, FVSpaceType >             DiscretizationType;

  static std::string static_id()
  {                                                                                // int() needed, otherwise we get a linker error
    return std::string("gdt.hyperbolic.discretization.fv.dim") + DSC::toString(int(GridType::dimension));
  }

  static DiscretizationType discretize(Stuff::Grid::ProviderInterface< GridType >& grid_provider,
                                       const ProblemType& problem,
                                       const int level = 0)
  {
    auto logger = Stuff::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    FVSpaceType space(GridViewType(grid_provider.level_view(level)));
    logger.debug() << "grid has " << space.grid_view().indexSet().size(0) << " elements" << std::endl;
    return DiscretizationType(problem, space);
  } // ... discretize(...)
}; // class CGDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_FV_HH
