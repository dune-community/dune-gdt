// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_BOCHNER_HH
#define DUNE_GDT_SPACES_BOCHNER_HH

#include <dune/grid/onedgrid.hh>

#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class BochnerSpace
{
public:
  template <class... TemporalGridArgs>
  BochnerSpace(const SpaceInterface<GV, r, rC, R>& spatial_space, TemporalGridArgs&&... temporal_grid_args)
    : spatial_space_(spatial_space)
    , temporal_grid_(std::forward<TemporalGridArgs>(temporal_grid_args)...)
    , temporal_space_(temporal_grid_.leafGridView(), /*order=*/1)
  {
  }

  const ContinuousLagrangeSpace<typename OneDGrid::LeafGridView>& temporal_space() const
  {
    return temporal_space_;
  }

  const SpaceInterface<GV, r, rC, R>& spatial_space() const
  {
    return spatial_space_;
  }

private:
  const SpaceInterface<GV, r, rC, R>& spatial_space_;
  const OneDGrid temporal_grid_;
  const ContinuousLagrangeSpace<typename OneDGrid::LeafGridView> temporal_space_;
}; // class BochnerSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BOCHNER_HH
