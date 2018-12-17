// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_HH

#include <dune/xt/grid/layers.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/type_traits.hh>

#include "projections/dirichlet.hh"
#include "projections/l2.hh"
#include "projections/lagrange.hh"

namespace Dune {
namespace GDT {


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value
                            && is_cg_space<SpaceType>::value && XT::LA::is_vector<VectorType>::value,
                        void>::type
project(const GridLayerType& grid_layer,
        const SourceType& source,
        DiscreteFunction<SpaceType, VectorType>& range,
        const size_t /*over_integrate*/ = 0)
{
  project_lagrange(grid_layer, source, range);
}

template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value
                            && !is_cg_space<SpaceType>::value && XT::LA::is_vector<VectorType>::value,
                        void>::type
project(const GridLayerType& grid_layer,
        const SourceType& source,
        DiscreteFunction<SpaceType, VectorType>& range,
        const size_t over_integrate = 0)
{
  project_l2(grid_layer, source, range, over_integrate);
}


template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_cg_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project(const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range, const size_t /*over_integrate*/ = 0)
{
  project_lagrange(source, range);
}

template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && !is_cg_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project(const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range, const size_t over_integrate = 0)
{
  project_l2(source, range, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_HH
