// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_INTERPOLATIONS_HH
#define DUNE_GDT_INTERPOLATIONS_HH

#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/discretefunction/default.hh>


namespace Dune {
namespace GDT {


/**
 * \brief Interpolates a localizable function within the space of a given discrete function.
 *
 * Simply uses the interpolation() of the target spaces finite_element().
 *
 * \note This might not be optimal for all spaces. For instance, the polynomial order of source is not taken into
 *       account for local L^2-projection based interpolation. This is a limitation in dune-localfunctions and we need
 *       to completely replace the interpolation of the respective local finite element.
 */
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::LocalizableFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target)
{
  target.dofs().vector().set_all(0);
  auto local_dof_vector = target.dofs().localize();
  auto local_source = source.local_function();
  std::vector<R> local_dofs(target.space().mapper().max_local_size());
  for (auto&& element : elements(target.space().grid_view())) {
    local_source->bind(element);
    local_dof_vector.bind(element);
    const auto& fe = target.space().finite_element(element.geometry().type());
    fe.interpolation().interpolate(
        [&](const auto& xx) { return local_source->evaluate(xx)[0]; }, local_source->order(), local_dofs);
    for (size_t ii = 0; ii < local_dof_vector.size(); ++ii)
      local_dof_vector[ii] = local_dofs[ii];
  }
} // ... interpolate(...)


template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::LocalizableFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  interpolate(source, target_function);
  return target_function;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_HH
