// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_INTERPOLATIONS_BOUNDARY_HH
#define DUNE_GDT_INTERPOLATIONS_BOUNDARY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/exceptions.hh>

#include "default.hh"

namespace Dune {
namespace GDT {


// ### Variants for GridFunctionInterface


/**
 * \brief Interpolates a localizable function restricted to certain boundary intersections within a given space [most
 *        general variant].
 *
 * Simply uses the interpolation() of the target spaces finite_element() and selects only those DoFs which lie on
 * intersections of correct boundary type.
 *
 * \note The same restrictions apply as for the default interpolations.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            DiscreteFunction<V, GV, r, rC, R>& target,
            const GridView<IGV>& interpolation_grid_view,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  using D = typename GridView<IGV>::ctype;
  static const constexpr int d = GridView<IGV>::dimension;
  auto local_dof_vector = target.dofs().localize();
  auto local_source = source.local_function();
  std::vector<R> local_dofs(target.space().mapper().max_local_size());
  for (auto&& element : elements(interpolation_grid_view)) {
    // first check if we should do something at all on this element
    size_t matching_boundary_intersections = 0;
    for (auto&& intersection : intersections(interpolation_grid_view, element))
      if (boundary_info.type(intersection) == target_boundary_type)
        ++matching_boundary_intersections;
    if (matching_boundary_intersections) {
      // some preparations
      local_source->bind(element);
      local_dof_vector.bind(element);
      const auto& fe = target.space().finite_element(element.geometry().type());
      // interpolate
      fe.interpolation().interpolate(
          [&](const auto& xx) { return local_source->evaluate(xx); }, local_source->order(), local_dofs);
      const auto& reference_element = ReferenceElements<D, d>::general(element.geometry().type());
      // but keep only those DoFs associated with the intersection, therefore
      // * determine these DoFs
      std::set<size_t> local_target_boundary_DoFs;
      const auto local_key_indices = fe.coefficients().local_key_indices();
      for (auto&& intersection : intersections(interpolation_grid_view, element)) {
        if (boundary_info.type(intersection) == target_boundary_type) {
          const auto intersection_index = intersection.indexInInside();
          for (const auto& local_DoF : local_key_indices[1][intersection_index])
            local_target_boundary_DoFs.insert(local_DoF);
          for (int cc = 2; cc <= d; ++cc) {
            for (int ii = 0; ii < reference_element.size(intersection_index, 1, cc); ++ii) {
              const auto subentity_id = reference_element.subEntity(intersection_index, 1, ii, cc);
              for (const auto& local_DoF : local_key_indices[cc][subentity_id])
                local_target_boundary_DoFs.insert(local_DoF);
            }
          }
        }
      }
      DUNE_THROW_IF(local_target_boundary_DoFs.size() == 0,
                    Exceptions::interpolation_error,
                    "The finite element is not suitable: although the element has an intersection on the boundary, no "
                    "DoFs are associated with it!");
      // * and use only those
      for (const auto& local_DoF_id : local_target_boundary_DoFs)
        local_dof_vector[local_DoF_id] = local_dofs[local_DoF_id];
    }
  }
} // ... interpolate(...)


/**
 * \brief Interpolates a localizable function restricted to certain boundary intersections within a given space [uses
 *        target.space().grid_view() as interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target,
                 const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
                 const XT::Grid::BoundaryType& target_boundary_type)
{
  interpolate(source, target, target.space().grid_view(), boundary_info, target_boundary_type);
}


/**
 * \brief Interpolates a localizable function restricted to certain boundary intersections within a given space [creates
 *        a suitable target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const GridView<IGV>& interpolation_grid_view,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  interpolate(source, target_function, interpolation_grid_view, boundary_info, target_boundary_type);
  return target_function;
}


/**
 * \brief Interpolates a localizable function restricted to certain boundary intersections within a given space [creates
 *        a suitable target function, uses target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  interpolate(source, target_function, boundary_info, target_boundary_type);
  return target_function;
}


// ### Variants for FunctionInterface


/**
 * \brief Interpolates a function restricted to certain boundary intersections within a given space [most general
 *        variant].
 *
 * Simply calls as_grid_function<>() and redirects to the appropriate interpolate() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(const XT::Functions::FunctionInterface<GridView<IGV>::dimension, r, rC, R>& source,
            DiscreteFunction<V, GV, r, rC, R>& target,
            const GridView<IGV>& interpolation_grid_view,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  interpolate(source.as_grid_function(interpolation_grid_view),
              target,
              interpolation_grid_view,
              boundary_info,
              target_boundary_type);
}


/**
 * \brief Interpolates a function restricted to certain boundary intersections within a given space [uses
 *        target.space().grid_view() as interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
                 DiscreteFunction<V, GV, r, rC, R>& target,
                 const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
                 const XT::Grid::BoundaryType& target_boundary_type)
{
  interpolate(source, target, target.space().grid_view(), boundary_info, target_boundary_type);
}


/**
 * \brief Interpolates a function restricted to certain boundary intersections within a given space [creates a suitable
 *        target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::FunctionInterface<GridView<IGV>::dimension, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const GridView<IGV>& interpolation_grid_view,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  return interpolate<VectorType>(source.as_grid_function(interpolation_grid_view),
                                 target_space,
                                 interpolation_grid_view,
                                 boundary_info,
                                 target_boundary_type);
}


/**
 * \brief Interpolates a function restricted to certain boundary intersections within a given space [creates a suitable
 *        target function, uses target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  return interpolate<VectorType>(source, target_space, target_space.grid_view(), boundary_info, target_boundary_type);
}


// ### Variants for GenericFunction


/**
 * \brief Interpolates a function given as a lambda expression restricted to certain boundary intersections within a
 *        given space [most general variant].
 *
 * Simply creates a XT::Functions::GenericFunction and redirects to the appropriate interpolate() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
                 void>
interpolate(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    DiscreteFunction<V, GV, r, rC, R>& target,
    const GridView<IGV>& interpolation_grid_view,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
    const XT::Grid::BoundaryType& target_boundary_type)
{
  interpolate(XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
              target,
              interpolation_grid_view,
              boundary_info,
              target_boundary_type);
}


/**
 * \brief Interpolates a function given as a lambda expression restricted to certain boundary intersections within a
 *        given space [uses target.space().grid_view() as interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void interpolate(const int source_order,
                 const std::function<typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::RangeReturnType(
                     const typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::DomainType&,
                     const XT::Common::Parameter&)> source_evaluate_lambda,
                 DiscreteFunction<V, GV, r, rC, R>& target,
                 const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
                 const XT::Grid::BoundaryType& target_boundary_type)
{
  interpolate(XT::Functions::GenericFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda),
              target,
              boundary_info,
              target_boundary_type);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGV>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGV::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    const SpaceInterface<GV, r, rC, R>& target_space,
    const GridView<IGV>& interpolation_grid_view,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridView<IGV>>>& boundary_info,
    const XT::Grid::BoundaryType& target_boundary_type)
{
  return interpolate<VectorType>(
      XT::Functions::GenericFunction<GridView<IGV>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
      target_space,
      interpolation_grid_view,
      boundary_info,
      target_boundary_type);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function, uses target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
interpolate(const int source_order,
            const std::function<typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::RangeReturnType(
                const typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::DomainType&,
                const XT::Common::Parameter&)> source_evaluate_lambda,
            const SpaceInterface<GV, r, rC, R>& target_space,
            const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GV>>& boundary_info,
            const XT::Grid::BoundaryType& target_boundary_type)
{
  return interpolate<VectorType>(
      XT::Functions::GenericFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda),
      target_space,
      boundary_info,
      target_boundary_type);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_BOUNDARY_HH
