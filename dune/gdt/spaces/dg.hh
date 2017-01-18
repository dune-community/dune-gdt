// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_DG_HH
#define DUNE_GDT_SPACES_DG_HH

#include <memory>

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include "interface.hh"
#include "dg/dune-fem-wrapper.hh"
#include "../playground/spaces/dg/dune-pdelab-wrapper.hh"


namespace Dune {
namespace GDT {


template <class GridType,
          XT::Grid::Layers layer_type,
          ChooseSpaceBackend backend_type,
          int polOrder,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1>
class DgSpaceProvider
{
  static const XT::Grid::Backends part_view_type = ChooseGridPartView<backend_type>::type;

public:
  typedef typename XT::Grid::Layer<GridType, layer_type, part_view_type>::type GridLayerType;

private:
  template <class G, int p, class R, size_t r, size_t rC, GDT::ChooseSpaceBackend b>
  struct SpaceChooser
  {
    static_assert(AlwaysFalse<G>::value, "No space available for this backend!");
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::fem>
  {
    typedef GDT::DuneFemDgSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::pdelab>
  {
    typedef GDT::DunePdelabDgSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

  typedef XT::Grid::GridProvider<GridType> GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface<GridType> MsGridProviderType;
#endif

public:
  typedef typename SpaceChooser<GridType, polOrder, RangeFieldType, dimRange, dimRangeCols, backend_type>::Type Type;
  typedef Type type;

  static Type create(GridLayerType grid_layer)
  {
    return Type(grid_layer);
  }

  static Type create(GridProviderType& grid_provider, const int level = 0)
  {
    return Type(grid_provider.template layer<layer_type, part_view_type>(level));
  }

#if HAVE_DUNE_GRID_MULTISCALE
  static Type create(const MsGridProviderType& grid_provider, const int level_or_subdomain = 0)
  {
    return Type(grid_provider.template layer<layer_type, part_view_type>(level_or_subdomain));
  }
#endif // HAVE_DUNE_GRID_MULTISCALE
}; // class DgSpaceProvider


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DG_HH
