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


template <class GV, size_t r, size_t rC, class R, class V, class IGVT>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            DiscreteFunction<V, GV, r, rC, R>& target,
            const GridView<IGVT>& interpolation_grid_view,
            const XT::Common::Parameter& param = {})
{
  DUNE_THROW_IF(target.space().type() == SpaceType::raviart_thomas && rC != 1 && r != GV::dimension,
                Exceptions::interpolation_error,
                "Raviart-Thomas interpolation does not make sense for these dimension!");
  if constexpr (rC == 1 && r == GV::dimension) {
    if (target.space().type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target, interpolation_grid_view, param);
    else
      default_interpolation(source, target, interpolation_grid_view, param);
  } else
    default_interpolation(source, target, interpolation_grid_view, param);
} // interpolate(...)


template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target,
                 const XT::Common::Parameter& param = {})
{
  DUNE_THROW_IF(target.space().type() == SpaceType::raviart_thomas && rC != 1 && r != GV::dimension,
                Exceptions::interpolation_error,
                "Raviart-Thomas interpolation does not make sense for these dimension!");
  if constexpr (rC == 1 && r == GV::dimension) {
    if (target.space().type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target, target.space().grid_view(), param);
    else
      default_interpolation(source, target, target.space().grid_view(), param);
  } else
    default_interpolation(source, target, target.space().grid_view(), param);
} // interpolate(...)


template <class VectorType, class GV, size_t r, size_t rC, class R, class IGVT>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const GridView<IGVT>& interpolation_grid_view,
            const XT::Common::Parameter& param = {})
{
  DUNE_THROW_IF(target_space.type() == SpaceType::raviart_thomas && rC != 1 && r != GV::dimension,
                Exceptions::interpolation_error,
                "Raviart-Thomas interpolation does not make sense for these dimension!");
  auto target_function = make_discrete_function<VectorType>(target_space);
  if constexpr (rC == 1 && r == GV::dimension) {
    if (target_space.type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target_function, interpolation_grid_view, param);
    else
      default_interpolation(source, target_function, interpolation_grid_view, param);
  } else
    default_interpolation(source, target_function, interpolation_grid_view, param);
  return target_function;
} // interpolate(...)


template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const XT::Common::Parameter& param = {})
{
  DUNE_THROW_IF(target_space.type() == SpaceType::raviart_thomas && rC != 1 && r != GV::dimension,
                Exceptions::interpolation_error,
                "Raviart-Thomas interpolation does not make sense for these dimension!");
  auto target_function = make_discrete_function<VectorType>(target_space);
  if constexpr (rC == 1 && r == GV::dimension) {
    if (target_space.type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target_function, target_space.grid_view(), param);
    else
      default_interpolation(source, target_function, target_space.grid_view(), param);
  } else
    default_interpolation(source, target_function, target_space.grid_view(), param);
  return target_function;
} // interpolate(...)


template <class GV, size_t r, size_t rC, class R>
DiscreteFunction<XT::LA::IstlDenseVector<R>, GV, r, rC, R>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const XT::Common::Parameter& param = {})
{
  DUNE_THROW_IF(target_space.type() == SpaceType::raviart_thomas && rC != 1 && r != GV::dimension,
                Exceptions::interpolation_error,
                "Raviart-Thomas interpolation does not make sense for these dimension!");
  auto target_function = make_discrete_function<XT::LA::IstlDenseVector<R>>(target_space);
  if constexpr (rC == 1 && r == GV::dimension) {
    if (target_space.type() == SpaceType::raviart_thomas)
      raviart_thomas_interpolation(source, target_function, target_space.grid_view(), param);
    else
      default_interpolation(source, target_function, target_space.grid_view(), param);
  } else
    default_interpolation(source, target_function, target_space.grid_view(), param);
  return target_function;
} // interpolate(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_HH
