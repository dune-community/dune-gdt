// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_DG_HH
#define DUNE_GDT_SPACES_DG_HH

#include <memory>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include "interface.hh"
#include "dg/default.hh"
#include "dg/dune-fem-wrapper.hh"
#include "../playground/spaces/dg/dune-pdelab-wrapper.hh"
#include "../playground/spaces/dg/dune-functions-wrapper.hh"
#include <dune/gdt/playground/spaces/block.hh>


namespace Dune {
namespace GDT {


template <class GridType,
          XT::Grid::Layers layer_type,
          Backends backend_type,
          int polOrder,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          XT::Grid::Backends grid_backend_type = layer_from_backend<backend_type>::type>
class DgSpaceProvider
{

public:
  static const constexpr SpaceType space_type = SpaceType::dg;
  static const constexpr Backends space_backend = backend_type;
  static const constexpr XT::Grid::Layers grid_layer = layer_type;
  static const constexpr XT::Grid::Backends layer_backend = grid_backend_type;

  typedef typename XT::Grid::Layer<GridType, layer_type, layer_backend>::type GridLayerType;

private:
  template <class G, int p, class R, size_t r, size_t rC, GDT::Backends b>
  struct SpaceChooser
  {
    static_assert(AlwaysFalse<G>::value, "No space available for this backend!");
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::Backends::fem>
  {
    typedef GDT::DuneFemDgSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::Backends::gdt>
  {
    static_assert(rC == 1, "");
    static_assert(r == 1, "");
    typedef GDT::DiscontinuousLagrangeSpace<GridLayerType, p, R> Type;
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::Backends::functions>
  {
    typedef GDT::DuneFunctionsDgSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::Backends::pdelab>
  {
    typedef GDT::DunePdelabDgSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

public:
  typedef typename SpaceChooser<GridType, polOrder, RangeFieldType, dimRange, dimRangeCols, backend_type>::Type Type;
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
}; // class DgSpaceProvider


template <class GridType,
          XT::Grid::Layers layer_type,
          Backends backend_type,
          int polOrder,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          XT::Grid::Backends grid_backend_type = layer_from_backend<backend_type>::type>
class BlockDgSpaceProvider
{
public:
  static const constexpr SpaceType space_type = SpaceType::block_dg;
  static const constexpr Backends space_backend = backend_type;
  static const constexpr XT::Grid::Layers grid_layer = layer_type;
  static const constexpr XT::Grid::Backends layer_backend = grid_backend_type;

  typedef typename XT::Grid::Layer<GridType, grid_layer, layer_backend>::type GridLayerType;

private:
  typedef typename DgSpaceProvider<GridType,
                                   grid_layer,
                                   space_backend,
                                   polOrder,
                                   RangeFieldType,
                                   dimRange,
                                   dimRangeCols,
                                   layer_backend>::Type LocalType;

public:
  typedef GDT::BlockSpace<LocalType> Type;
  typedef Type type;

  static Type create(GridLayerType /*grid_layer*/)
  {
    DUNE_THROW(NotImplemented, "Only usable to extract the correct type");
  }

  template <class DdGridType>
  static Type create(XT::Grid::GridProvider<GridType, DdGridType>& /*grid_provider*/, const int /*level*/ = 0)
  {
    DUNE_THROW(NotImplemented, "Only usable to extract the correct type");
  }
}; // class BlockDgSpaceProvider


} // namespace GDT
} // namespace Dune


#include "dg.lib.hh"


#endif // DUNE_GDT_SPACES_DG_HH
