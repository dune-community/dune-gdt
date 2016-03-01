// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_HH

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/layers.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/cg/interface.hh>

#include "projections/dirichlet.hh"
#include "projections/l2.hh"
#include "projections/lagrange.hh"

namespace Dune {
namespace GDT {


template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && is_cg_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project(const GridViewType& grid_view,
        const SourceType& source,
        DiscreteFunction< SpaceType, VectorType >& range,
        const size_t /*over_integrate*/ = 0)
{
  project_lagrange(grid_view, source, range);
}

template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && !is_cg_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project(const GridViewType& grid_view,
        const SourceType& source,
        DiscreteFunction< SpaceType, VectorType >& range,
        const size_t over_integrate = 0)
{
  project_l2(grid_view, source, range, over_integrate);
}


template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && is_cg_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project(const SourceType& source, DiscreteFunction< SpaceType, VectorType >& range, const size_t /*over_integrate*/ = 0)
{
  project_lagrange(source, range);
}

template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && !is_cg_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project(const SourceType& source, DiscreteFunction< SpaceType, VectorType >& range, const size_t over_integrate = 0)
{
  project_l2(source, range, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_HH
