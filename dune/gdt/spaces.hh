#ifndef DUNE_GDT_SPACES_HH
#define DUNE_GDT_SPACES_HH

#include "spaces/cg.hh"
#include "spaces/dg.hh"
#include "spaces/fv.hh"
#include "spaces/rt.hh"

namespace Dune {
namespace GDT {


template <class GridType,
          XT::Grid::Layers layer_type,
          SpaceType space_type,
          Backends backend_type,
          int polOrder,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1>
class SpaceProvider
{
  static_assert(AlwaysFalse<GridType>::value, "Please add a specialization for this space type!");
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC>
class SpaceProvider<G, layer, SpaceType::cg, backend, p, R, r, rC>
    : public CgSpaceProvider<G, layer, backend, p, R, r, rC>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC>
class SpaceProvider<G, layer, SpaceType::dg, backend, p, R, r, rC>
    : public DgSpaceProvider<G, layer, backend, p, R, r, rC>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC>
class SpaceProvider<G, layer, SpaceType::fv, backend, p, R, r, rC> : public FvSpaceProvider<G, layer, backend, R, r, rC>
{
  static_assert(p == 0, "There is no FV space with nonzero polOrder!");
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC>
class SpaceProvider<G, layer, SpaceType::rt, backend, p, R, r, rC>
    : public RtSpaceProvider<G, layer, backend, p, R, r, rC>
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_HH
