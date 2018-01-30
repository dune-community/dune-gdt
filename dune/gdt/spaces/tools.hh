// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014, 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACE_TOOLS_HH
#define DUNE_GDT_SPACE_TOOLS_HH

#warning This header is deprecated, all information is in the space now (04.04.2017)!

#include <memory>
#include <type_traits>

#include "interface.hh"

#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace GDT {
namespace SpaceTools {


template <class GridType, bool view = true>
struct DUNE_DEPRECATED_MSG("Do not use this any more, all information is in the space (04.04.2017)!") LeafGridPartView
{
  typedef typename GridType::LeafGridView Type;

  static Type create(GridType& grid)
  {
    return Type(grid.leafGridView());
  }
}; // struct LeafGridPartView< ..., true >


template <class GridType, bool view = true>
struct DUNE_DEPRECATED_MSG("Do not use this any more, all information is in the space (04.04.2017)!") LevelGridPartView
{
  typedef typename GridType::LevelGridView Type;

  static Type create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return Type(grid.levelGridView(level));
  }
}; // struct LevelGridPartView< ..., true >


template <class SpaceType>
class DUNE_DEPRECATED_MSG("Do not use this any more, all information is in the space (04.04.2017)!") GridPartView
{
  static_assert(std::is_base_of<SpaceInterface<typename SpaceType::Traits,
                                               SpaceType::dimDomain,
                                               SpaceType::dimRange,
                                               SpaceType::dimRangeCols>,
                                SpaceType>::value,
                "SpaceType has to be derived from SpaceInterface!");
  static const bool needs_grid_view = (SpaceType::layer_backend == XT::Grid::Backends::view);

public:
  using GridType = XT::Grid::extract_grid_t<typename SpaceType::GridLayerType>;
  typedef typename LeafGridPartView<GridType, needs_grid_view>::Type LeafGridLayerType;
  typedef typename LevelGridPartView<GridType, needs_grid_view>::Type LevelGridLayerType;

  static LeafGridLayerType create_leaf(GridType& grid)
  {
    return LeafGridPartView<GridType, needs_grid_view>::create(grid);
  }

  static LevelGridLayerType create_level(GridType& grid, const int level)
  {
    return LevelGridPartView<GridType, needs_grid_view>::create(grid, level);
  }
}; // struct GridPartView


} // namespace SpaceTools
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_TOOLS_HH
