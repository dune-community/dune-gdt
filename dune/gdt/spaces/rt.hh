// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2015)

#ifndef DUNE_GDT_SPACES_RT_HH
#define DUNE_GDT_SPACES_RT_HH

#include "rt/interface.hh"
#include "rt/dune-pdelab-wrapper.hh"

namespace Dune {
namespace GDT {


template <class GridType,
          XT::Grid::Layers layer_type,
          ChooseSpaceBackend backend_type,
          int polOrder,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1>
class RtSpaceProvider
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
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::pdelab>
  {
    typedef GDT::DunePdelabRtSpaceWrapper<GridLayerType, p, R, r, rC> Type;
  };

  typedef XT::Grid::GridProvider<GridType> GridProviderType;

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
}; // class RtSpaceProvider


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_HH
