// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_HH

#include <memory>

#include <dune/stuff/grid/partview.hh>
#include <dune/stuff/grid/provider/interface.hh>

#include "interface.hh"
#include "continuouslagrange/fem.hh"
#include "continuouslagrange/pdelab.hh"


namespace Dune {
namespace GDT {
namespace Spaces {


template <class GridType, Stuff::Grid::ChooseLayer layer_type, ChooseSpaceBackend backend_type, int polOrder,
          class RangeFieldType, int dimRange, int dimRangeCols = 1>
class ContinuousLagrangeProvider
{
  static const Stuff::Grid::ChoosePartView part_view_type = ChooseGridPartView<backend_type>::type;

public:
  typedef typename Stuff::Grid::Layer<GridType, layer_type, part_view_type>::Type GridLayerType;

private:
  template <class G, int p, class R, int r, int rC, GDT::ChooseSpaceBackend b>
  struct SpaceChooser
  {
    static_assert(AlwaysFalse<G>::value, "No space available for this backend!");
  };

  template <class G, int p, class R, int r, int rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::fem>
  {
    typedef GDT::Spaces::ContinuousLagrange::FemBased<GridLayerType, p, R, r> Type;
  };

  template <class G, int p, class R, int r, int rC>
  struct SpaceChooser<G, p, R, r, rC, GDT::ChooseSpaceBackend::pdelab>
  {
    typedef GDT::Spaces::ContinuousLagrange::PdelabBased<GridLayerType, p, R, r> Type;
  };

public:
  typedef typename SpaceChooser<GridType, polOrder, RangeFieldType, dimRange, dimRangeCols, backend_type>::Type Type;

  static Type create(const std::shared_ptr<const GridLayerType> grid_layer)
  {
    return Type(grid_layer);
  }

  static Type create(const Stuff::Grid::ConstProviderInterface<GridType>& grid_provider, const int level = 0)
  {
    return Type(grid_provider.template layer<layer_type, part_view_type>(level));
  }
}; // class ContinuousLagrangeProvider


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_HH
