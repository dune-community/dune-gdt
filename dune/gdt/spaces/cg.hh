// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CG_HH
#define DUNE_GDT_SPACES_CG_HH

#include <memory>

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider/interface.hh>

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/provider/interface.hh>
#endif

#include "interface.hh"
#include "cg/fem.hh"
#include "cg/pdelab.hh"


namespace Dune {
namespace GDT {
namespace Spaces {


template <class GridType, Stuff::Grid::ChooseLayer layer_type, ChooseSpaceBackend backend_type, int polOrder,
          class RangeFieldType, size_t dimRange, size_t dimRangeCols = 1>
class CGProvider
{
  static const Stuff::Grid::ChoosePartView part_view_type = ChooseGridPartView<backend_type>::type;

public:
  typedef typename Stuff::Grid::Layer<GridType, layer_type, part_view_type>::Type GridLayerType;

private:
  template <class G, int p, class R, size_t r, size_t rC, GDT::ChooseSpaceBackend b>
  struct SpaceChooser
  {
    static_assert(AlwaysFalse<G>::value, "No space available for this backend!");
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::fem>
  {
    typedef GDT::Spaces::CG::FemBased<GridLayerType, p, R, r> Type;
  };

  template <class G, int p, class R, size_t r, size_t rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::pdelab>
  {
    typedef GDT::Spaces::CG::PdelabBased<GridLayerType, p, R, r> Type;
  };

  typedef Stuff::Grid::ProviderInterface<GridType> GridProviderType;
#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface<GridType> MsGridProviderType;
#endif

public:
  typedef typename SpaceChooser<GridType, polOrder, RangeFieldType, dimRange, dimRangeCols, backend_type>::Type Type;

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
}; // class CGProvider


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_HH
