// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_TOOLS_HYPERBOLIC_HH
#define DUNE_GDT_TOOLS_HYPERBOLIC_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Estimates dt via [Cockburn, Coquel, LeFloch, 1995]
 */
template <class GV,
          size_t m_as_size_t,
          class R,
          int m_as_int = m_as_size_t,
          typename = std::enable_if_t<ssize_t(m_as_size_t) == ssize_t(m_as_int), void>>
double estimate_dt_for_hyperbolic_system(
    const GridView<GV>& grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridView<GV>>, m_as_size_t, 1, R>& state,
    const XT::Functions::FunctionInterface<m_as_size_t, GridView<GV>::dimension, m_as_size_t, R>& flux,
    const std::pair<XT::Common::FieldVector<R, m_as_int>, XT::Common::FieldVector<R, m_as_int>>& boundary_data_range = {
        XT::Common::FieldVector<R, m_as_int>(std::numeric_limits<R>::max()),
        XT::Common::FieldVector<R, m_as_int>(std::numeric_limits<R>::min())})
{
  using D = typename GridView<GV>::ctype;
  static const constexpr size_t d = GridView<GV>::dimension;
  static const constexpr size_t m = m_as_size_t;
  // estimate data range
  auto data_minimum = boundary_data_range.first;
  auto data_maximum = boundary_data_range.second;
  auto local_state = state.local_function();
  for (auto&& element : elements(grid_view)) {
    local_state->bind(element);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.type(), local_state->order())) {
      const auto state_value = local_state->evaluate(quadrature_point.position());
      for (size_t ii = 0; ii < m; ++ii) {
        data_minimum[ii] = std::min(data_minimum[ii], state_value[ii]);
        data_maximum[ii] = std::max(data_maximum[ii], state_value[ii]);
      }
    }
  }
  // ensure distinct minima/maxima (otherwise the grid creation below will fail)
  for (size_t ii = 0; ii < m; ++ii)
    if (!(data_minimum[ii] < data_maximum[ii]))
      data_maximum[ii] = data_minimum[ii] + 1e-6 * data_minimum[ii];
  // estimate flux derivative range
  R max_flux_derivative = std::numeric_limits<R>::min();
  const auto flux_range_grid = XT::Grid::make_cube_grid<YaspGrid<m, EquidistantOffsetCoordinates<double, m>>>(
      data_minimum, data_maximum, XT::Common::FieldVector<unsigned int, m>(1));
  const auto flux_range = *flux_range_grid.leaf_view().template begin<0>();
  for (const auto& quadrature_point : QuadratureRules<R, m>::rule(flux_range.type(), flux.order())) {
    const auto df = flux.jacobian(flux_range.geometry().global(quadrature_point.position()));
    for (size_t ss = 0; ss < d; ++ss)
      max_flux_derivative = std::max(max_flux_derivative, df[ss].infinity_norm());
  }
  D perimeter_over_volume = std::numeric_limits<D>::min();
  for (auto&& element : elements(grid_view)) {
    D perimeter = 0;
    for (auto&& intersection : intersections(grid_view, element))
      perimeter += intersection.geometry().volume();
    perimeter_over_volume = std::max(perimeter_over_volume, perimeter / element.geometry().volume());
  }
  const auto dt = 1. / (perimeter_over_volume * max_flux_derivative);
  return dt;
} // ... estimate_dt_for_hyperbolic_system(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_HYPERBOLIC_HH
