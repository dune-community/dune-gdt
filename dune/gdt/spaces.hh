// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

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
          size_t dimRangeCols = 1,
          XT::Grid::Backends grid_backend_type = layer_from_backend<backend_type>::type>
class SpaceProvider
{
  static_assert(AlwaysFalse<GridType>::value, "Please add a specialization for this space type!");
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::cg, backend, p, R, r, rC, g>
    : public CgSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::block_cg, backend, p, R, r, rC, g>
    : public BlockCgSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::dg, backend, p, R, r, rC, g>
    : public DgSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::block_dg, backend, p, R, r, rC, g>
    : public BlockDgSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::fv, backend, p, R, r, rC, g>
    : public FvSpaceProvider<G, layer, backend, R, r, rC, g>
{
  static_assert(p == 0, "There is no FV space with nonzero polOrder!");
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::block_fv, backend, p, R, r, rC, g>
    : public BlockFvSpaceProvider<G, layer, backend, R, r, rC, g>
{
  static_assert(p == 0, "There is no FV space with nonzero polOrder!");
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::rt, backend, p, R, r, rC, g>
    : public RtSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, class R, size_t r, size_t rC, XT::Grid::Backends g>
class SpaceProvider<G, layer, SpaceType::block_rt, backend, p, R, r, rC, g>
    : public BlockRtSpaceProvider<G, layer, backend, p, R, r, rC, g>
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_HH
