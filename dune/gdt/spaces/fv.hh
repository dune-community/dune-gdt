// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2015, 2017)

#ifndef DUNE_GDT_SPACES_FV_HH
#define DUNE_GDT_SPACES_FV_HH

#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/layers.hh>

#include "fv/interface.hh"
#include "fv/default.hh"

#include <dune/xt/grid/gridprovider.hh>

namespace Dune {
namespace GDT {


template <class GridType,
          XT::Grid::Layers layer_type,
          Backends backend_type,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1>
class FvSpaceProvider
{
public:
  static const constexpr SpaceType space_type = SpaceType::fv;
  static const constexpr Backends space_backend = backend_type;
  static const constexpr XT::Grid::Layers grid_layer = layer_type;
  static const constexpr XT::Grid::Backends layer_backend = layer_from_backend<backend_type>::type;

  typedef typename XT::Grid::Layer<GridType, layer_type, layer_backend>::type GridLayerType;

private:
  template <class G, class R, size_t r, size_t rC, GDT::Backends b>
  struct SpaceChooser
  {
    static_assert(AlwaysFalse<G>::value, "No space available for this backend!");
  };

  template <class G, class R, size_t r, size_t rC>
  struct SpaceChooser<G, R, r, rC, GDT::Backends::gdt>
  {
    typedef GDT::FvSpace<GridLayerType, R, r, rC> Type;
  };

public:
  typedef typename SpaceChooser<GridType, RangeFieldType, dimRange, dimRangeCols, backend_type>::Type Type;
  typedef Type type;

  static Type create(GridLayerType grid_layer)
  {
    return Type(grid_layer);
  }

  template <class DdGridType>
  static Type create(XT::Grid::GridProvider<GridType, DdGridType>& grid_provider, const int level = 0)
  {
    return Type(grid_provider.template layer<layer_type, layer_backend>(level));
  }
}; // class FvSpaceProvider


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_HH
