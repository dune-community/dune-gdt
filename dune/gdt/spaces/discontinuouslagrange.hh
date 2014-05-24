// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_HH
#define DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_HH

#include <memory>

#if HAVE_DUNE_GRID_MULTISCALE
# include <dune/grid/multiscale/provider/interface.hh>
#endif

#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider/interface.hh>

#include "interface.hh"
#include "discontinuouslagrange/fem-localfunctions.hh"


namespace Dune {
namespace GDT {
namespace Spaces {


template< class GridType, Stuff::Grid::ChooseLayer layer_type, ChooseSpaceBackend backend_type,
          int polOrder, class RangeFieldType, int dimRange, int dimRangeCols = 1 >
class DiscontinuousLagrangeProvider
{
  static const Stuff::Grid::ChoosePartView part_view_type = ChooseGridPartView< backend_type >::type;
public:
  typedef typename Stuff::Grid::Layer< GridType, layer_type, part_view_type >::Type GridLayerType;

private:
  template< class G, int p, class R, int r, int rC, GDT::ChooseSpaceBackend b >
  struct SpaceChooser
  {
    static_assert(AlwaysFalse< G >::value, "No space available for this backend!");
  };

  template< class G, int p, class R, int r, int rC >
  struct SpaceChooser< G, p, R, r, rC, GDT::ChooseSpaceBackend::fem_localfunctions >
  {
    typedef GDT::Spaces::DiscontinuousLagrange::FemLocalfunctionsBased< GridLayerType, p, R, r > Type;
  };

#if HAVE_DUNE_GRID_MULTISCALE
  typedef grid::Multiscale::ProviderInterface< GridType > GridProviderType;
#else
  typedef Stuff::Grid::ConstProviderInterface< GridType > GridProviderType;
#endif

public:
  typedef typename SpaceChooser< GridType, polOrder, RangeFieldType, dimRange, dimRangeCols, backend_type >::Type Type;

  static Type create(const std::shared_ptr< const GridLayerType > grid_layer)
  {
    return Type(grid_layer);
  }

  static Type create(const GridProviderType& grid_provider, const int level_or_subdomain = 0)
  {
    return Type(grid_provider.template layer< layer_type, part_view_type >(level_or_subdomain));
  }
}; // class DiscontinuousLagrangeProvider


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_HH
