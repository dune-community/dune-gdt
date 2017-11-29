// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_PROLONGATIONS_HH
#define DUNE_GDT_PROLONGATIONS_HH

#include <dune/xt/grid/layers.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/cg/interface.hh>

#include "prolongations/l2.hh"
#include "prolongations/lagrange.hh"

namespace Dune {
namespace GDT {


template <class GridLayerType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value && is_cg_space<RangeSpaceType>::value, void>::type
prolong(const GridLayerType& grid_layer,
        const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
        const size_t /*over_integrate*/ = 0)
{
  prolong_lagrange(grid_layer, source, range);
}

template <class GridLayerType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value && !is_cg_space<RangeSpaceType>::value, void>::type
prolong(const GridLayerType& grid_layer,
        const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
        const size_t over_integrate = 0)
{
  prolong_l2(grid_layer, source, range, over_integrate);
}


template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
typename std::enable_if<is_cg_space<RangeSpaceType>::value, void>::type
prolong(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
        const size_t /*over_integrate*/ = 0)
{
  prolong_lagrange(source, range);
}

template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
typename std::enable_if<!is_cg_space<RangeSpaceType>::value, void>::type
prolong(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
        const size_t over_integrate = 0)
{
  prolong_l2(source, range, over_integrate);
}


template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
void prolong(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
             DiscreteFunction<BlockSpace<RangeSpaceType>, RangeVectorType>& range,
             const size_t over_integrate = 0)
{
  const auto& block_space = range.space();
  for (size_t subdomain = 0; subdomain < block_space.dd_grid().size(); ++subdomain) {
    const auto& local_space = block_space.local_space(subdomain);
    auto subdomain_range = make_discrete_function<RangeVectorType>(local_space);
    prolong(source, subdomain_range, over_integrate);
    for (size_t ii = 0; ii < local_space.mapper().size(); ++ii)
      range.vector()[block_space.mapper().mapToGlobal(subdomain, ii)] = subdomain_range.vector()[ii];
  }
} // ... prolong(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_HH
