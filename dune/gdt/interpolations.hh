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
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/default.hh>


namespace Dune {
namespace GDT {


// ### Variants for GridFunctionInterface


/**
 * \brief Interpolates a localizable function within a given space [most general variant].
 *
 * Simply uses the interpolation() of the target spaces finite_element().
 *
 *
 * \note This does not clear target.dofs().vector(). Thus, if interpolation_grid_view only covers a part of the domain
 *       of target.space().grid_view(), other contributions in target remain (which is on purpose).
 *
 * \note This might not be optimal for all spaces. For instance, the polynomial order of source is not taken into
 *       account for local L^2-projection based interpolation. This is a limitation in dune-localfunctions and we need
 *       to completely replace the interpolation of the respective local finite element.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            DiscreteFunction<V, GV, r, rC, R>& target,
            const GridView<IGV>& interpolation_grid_view)
{
  auto local_dof_vector = target.dofs().localize();
  auto local_source = source.local_function();
  std::vector<R> local_dofs(target.space().mapper().max_local_size());
  for (auto&& element : elements(interpolation_grid_view)) {
    local_source->bind(element);
    local_dof_vector.bind(element);
    const auto& fe = target.space().finite_element(element.geometry().type());
    fe.interpolation().interpolate(
        [&](const auto& xx) { return local_source->evaluate(xx); }, local_source->order(), local_dofs);
    for (size_t ii = 0; ii < local_dof_vector.size(); ++ii)
      local_dof_vector[ii] = local_dofs[ii];
  }
} // ... interpolate(...)


/**
 * \brief Interpolates a localizable function within a given space [uses target.space().grid_view() as
 *        interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target)
{
  interpolate(source, target, target.space().grid_view());
}


/**
 * \brief Interpolates a localizable function within a given space [creates a suitable target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value
                     && std::is_same<XT::Grid::extract_entity_t<GV>,
                                     typename IGV::Grid::template Codim<0>::Entity>::value,
                 DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const GridView<IGV>& interpolation_grid_view)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  interpolate(source, target_function, interpolation_grid_view);
  return target_function;
}


/**
 * \brief Interpolates a localizable function within a given space [creates a suitable target function, uses
 *        target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  interpolate(source, target_function);
  return target_function;
}


// ### Variants for FunctionInterface


/**
 * \brief Interpolates a function within a given space [most general variant].
 *
 * Simply calls as_grid_function<>() and redirects to the appropriate interpolate() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(const XT::Functions::FunctionInterface<GridView<IGV>::dimension, r, rC, R>& source,
            DiscreteFunction<V, GV, r, rC, R>& target,
            const GridView<IGV>& interpolation_grid_view)
{
  interpolate(source.as_grid_function(interpolation_grid_view), target, interpolation_grid_view);
}


/**
 * \brief Interpolates a function within a given space [uses target.space().grid_view() as
 *        interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target)
{
  interpolate(source, target, target.space().grid_view());
}


/**
 * \brief Interpolates a function within a given space [creates a suitable target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value
                     && std::is_same<XT::Grid::extract_entity_t<GV>,
                                     typename IGV::Grid::template Codim<0>::Entity>::value,
                 DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::FunctionInterface<GridView<IGV>::dimension, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const GridView<IGV>& interpolation_grid_view)
{
  return interpolate<VectorType>(
      source.as_grid_function(interpolation_grid_view), target_space, interpolation_grid_view);
}


/**
 * \brief Interpolates a function within a given space [creates a suitable target function, uses
 *        target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space)
{
  return interpolate<VectorType>(source, target_space, target_space.grid_view());
}


// ### Variants for LambdaFunction


/**
 * \brief Interpolates a function given as a lambda expression within a given space [most general variant].
 *
 * Simply creates a XT::Functions::LambdaFunction and redirects to the appropriate interpolate() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(
    const int source_order,
    const std::function<typename XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    DiscreteFunction<V, GV, r, rC, R>& target,
    const GridView<IGV>& interpolation_grid_view)
{
  interpolate(XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
              target,
              interpolation_grid_view);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [uses
 *        target.space().grid_view() as interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const int source_order,
                 const std::function<typename XT::Functions::LambdaFunction<GV::dimension, r, rC, R>::RangeReturnType(
                     const typename XT::Functions::LambdaFunction<GV::dimension, r, rC, R>::DomainType&,
                     const XT::Common::Parameter&)> source_evaluate_lambda,
                 DiscreteFunction<V, GV, r, rC, R>& target)
{
  interpolate(XT::Functions::LambdaFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda), target);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value
                     && std::is_same<XT::Grid::extract_entity_t<GV>,
                                     typename IGV::Grid::template Codim<0>::Entity>::value,
                 DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(
    const int source_order,
    const std::function<typename XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    const SpaceInterface<GV, r, rC, R>& target_space,
    const GridView<IGV>& interpolation_grid_view)
{
  return interpolate<VectorType>(
      XT::Functions::LambdaFunction<GridView<IGV>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
      target_space,
      interpolation_grid_view);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function, uses target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const int source_order,
            const std::function<typename XT::Functions::LambdaFunction<GV::dimension, r, rC, R>::RangeReturnType(
                const typename XT::Functions::LambdaFunction<GV::dimension, r, rC, R>::DomainType&,
                const XT::Common::Parameter&)> source_evaluate_lambda,
            const SpaceInterface<GV, r, rC, R>& target_space)
{
  return interpolate<VectorType>(
      XT::Functions::LambdaFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda), target_space);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_HH
