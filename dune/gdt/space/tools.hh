// This file is view of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_TOOLS_HH
#define DUNE_GDT_SPACE_TOOLS_HH

#include <memory>
#include <type_traits>

#include "interface.hh"

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/levelgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#endif // HAVE_DUNE_FEM

namespace Dune {
namespace GDT {
namespace SpaceTools {


template <class GridType, bool view = true>
struct LeafGridPartView
{
  typedef typename GridType::LeafGridView Type;

  static std::shared_ptr<Type> create(GridType& grid)
  {
    return std::make_shared<Type>(grid.leafGridView());
  }
}; // struct LeafGridPartView< ..., true >


template <class GridType, bool view = true>
struct LevelGridPartView
{
  typedef typename GridType::LevelGridView Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid.levelGridView(level));
  }
}; // struct LevelGridPartView< ..., true >


#if HAVE_DUNE_FEM


template <class GridType>
struct LeafGridPartView<GridType, false>
{
  typedef Dune::Fem::LeafGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid)
  {
    return std::make_shared<Type>(grid);
  }
}; // struct LeafGridPartView< ..., false >


template <class GridType>
struct LevelGridPartView<GridType, false>
{
  typedef Dune::Fem::LevelGridPart<GridType> Type;

  static std::shared_ptr<Type> create(GridType& grid, const int level)
  {
    assert(level >= 0);
    assert(level <= grid.maxLevel());
    return std::make_shared<Type>(grid, level);
  }
}; // struct LevelGridPartView< ..., false >


#endif // HAVE_DUNE_FEM


template <class SpaceType>
class GridPartView
{
  static_assert(std::is_base_of<SpaceInterface<typename SpaceType::Traits>, SpaceType>::value,
                "SpaceType has to be derived from SpaceInterface!");
  static const bool needs_grid_view = SpaceType::needs_grid_view;

public:
  typedef typename SpaceType::GridViewType::Grid GridType;
  typedef typename LeafGridPartView<GridType, needs_grid_view>::Type LeafGridViewType;
  typedef typename LevelGridPartView<GridType, needs_grid_view>::Type LevelGridViewType;

  static std::shared_ptr<LeafGridViewType> create_leaf(GridType& grid)
  {
    return LeafGridPartView<GridType, needs_grid_view>::create(grid);
  }

  static std::shared_ptr<LevelGridViewType> create_level(GridType& grid, const int level)
  {
    return LevelGridPartView<GridType, needs_grid_view>::create(grid, level);
  }
}; // struct GridPartView


} // namespace SpaceTools
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_TOOLS_HH
