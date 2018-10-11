// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk       (2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_FV_FROM_GRIDFUNCTOR_HH
#define DUNE_GDT_DISCRETEFUNCTION_FV_FROM_GRIDFUNCTOR_HH

#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/fv.hh>

namespace Dune {
namespace GDT {

template <class VectorType, class GridLayerType, class RangeFieldType>
std::pair<FvSpace<GridLayerType, RangeFieldType, 1>, VectorType>
evaluate_into_fv_vector(const GridLayerType& grid_layer,
                        XT::Grid::Functor::Codim0Return<GridLayerType, RangeFieldType>& functor,
                        std::string name,
                        bool visualize = false)
{
  using Space = FvSpace<GridLayerType, RangeFieldType, 1>;
  const Space space{grid_layer};
  auto function = make_discrete_function<VectorType>(space, name);
  XT::Grid::Walker<GridLayerType> grid_walker(grid_layer);
  auto wrap = [&](const auto& element) {
    auto el_func = function.local_discrete_function(element);
    el_func->vector().set(0, functor.compute_locally(element));
  };
  grid_walker.append(wrap);
  grid_walker.walk();
  if (visualize) {
    function.visualize(name, false);
  }
  return std::make_pair(space, function.vector());
}
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_FV_FROM_GRIDFUNCTOR_HH
