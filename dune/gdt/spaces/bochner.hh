// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_SPACES_BOCHNER_HH
#define DUNE_GDT_SPACES_BOCHNER_HH

#include <dune/grid/onedgrid.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/ranges.hh>
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
    , time_interval_(std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::min()))
  {
    for (auto&& time_interval : elements(temporal_space_.grid_view())) {
      for (auto ii : XT::Common::value_range(time_interval.geometry().corners())) {
        const auto time_point = time_interval.geometry().corner(ii);
        time_interval_.first = std::min(time_interval_.first, time_point[0]);
        time_interval_.second = std::max(time_interval_.second, time_point[0]);
      }
    }
  }

  const ContinuousLagrangeSpace<typename OneDGrid::LeafGridView>& temporal_space() const
  {
    return temporal_space_;
  }

  const SpaceInterface<GV, r, rC, R>& spatial_space() const
  {
    return spatial_space_;
  }

  const std::pair<double, double>& time_interval() const
  {
    return time_interval_;
  }

  std::vector<double> time_points() const
  {
    auto temporal_basis = temporal_space_.basis().localize();
    DynamicVector<size_t> global_DoF_indices(temporal_space_.mapper().max_local_size());
    std::vector<double> points(temporal_space_.mapper().size(), 0.);
    for (auto&& time_interval : elements(temporal_space_.grid_view())) {
      temporal_basis->bind(time_interval);
      temporal_space_.mapper().global_indices(time_interval, global_DoF_indices);
      auto local_lagrange_points = temporal_basis->finite_element().lagrange_points();
      for (size_t ii = 0; ii < temporal_basis->size(); ++ii)
        points[global_DoF_indices[ii]] = time_interval.geometry().global(local_lagrange_points[ii])[0];
    }
    return points;
  } // ... time_points(...)

private:
  const SpaceInterface<GV, r, rC, R>& spatial_space_;
  const OneDGrid temporal_grid_;
  const ContinuousLagrangeSpace<typename OneDGrid::LeafGridView> temporal_space_;
  std::pair<double, double> time_interval_;
}; // class BochnerSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BOCHNER_HH
