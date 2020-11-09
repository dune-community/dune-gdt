// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_INTERPOLATIONS_HH
#define DUNE_GDT_INTERPOLATIONS_HH

#include "interpolations/boundary.hh"
#include "interpolations/default.hh"
#include "interpolations/raviart-thomas.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r, class R, class V, class GV, class IGV>
void interpolate(const XT::Functions::GridFunctionInterface<E, r, 1, R>& source,
                 DiscreteFunction<V, GV, r, 1, R>& target,
                 const GridView<IGV>& interpolation_grid_view,
                 const XT::Common::Parameter& param = {})
{
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GV>>, "");
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GridView<IGV>>>, "");

  if constexpr (r == GV::dimension) {
    if (target.space().type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target, interpolation_grid_view, param);
    else
      default_interpolation(source, target, interpolation_grid_view, param);
  } else
    default_interpolation(source, target, interpolation_grid_view, param);
} // ... interpolate(...)


template <class E, size_t r, class R, class V, class GV>
void interpolate(const XT::Functions::GridFunctionInterface<E, r, 1, R>& source,
                 DiscreteFunction<V, GV, r, 1, R>& target,
                 const XT::Common::Parameter& param = {})
{
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GV>>, "");

  if constexpr (r == GV::dimension) {
    if (target.space().type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target, param);
    else
      default_interpolation(source, target, param);
  } else
    default_interpolation(source, target, param);
} // ... interpolate(...)


template <class VectorType, class GV, size_t r, class R, class E, class IGV>
DiscreteFunction<VectorType, GV, r, 1, R> interpolate(const XT::Functions::GridFunctionInterface<E, r, 1, R>& source,
                                                      const SpaceInterface<GV, r, 1, R>& target_space,
                                                      const GridView<IGV>& interpolation_grid_view,
                                                      const XT::Common::Parameter& param = {})
{
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GV>>, "");
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GridView<IGV>>>, "");

  if constexpr (r == GV::dimension) {
    if (target_space.type() == SpaceType::raviart_thomas)
      return raviart_thomas_interpolation<VectorType>(source, target_space, interpolation_grid_view, param);
    else
      return default_interpolation<VectorType>(source, target_space, interpolation_grid_view, param);
  } else
    return default_interpolation<VectorType>(source, target_space, interpolation_grid_view, param);
} // ... interpolate(...)


template <class VectorType, class GV, size_t r, class R, class E>
DiscreteFunction<VectorType, GV, r, 1, R> interpolate(const XT::Functions::GridFunctionInterface<E, r, 1, R>& source,
                                                      const SpaceInterface<GV, r, 1, R>& target_space,
                                                      const XT::Common::Parameter& param = {})
{
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GV>>, "");

  if constexpr (r == GV::dimension) {
    if (target_space.type() == SpaceType::raviart_thomas)
      return raviart_thomas_interpolation<VectorType>(source, target_space, param);
    else
      return default_interpolation<VectorType>(source, target_space, param);
  } else
    return default_interpolation<VectorType>(source, target_space, param);
} // ... interpolate(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_HH
