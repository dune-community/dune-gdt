// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_INTERPOLATIONS_RAVIART_THOMAS_HH
#define DUNE_GDT_INTERPOLATIONS_RAVIART_THOMAS_HH

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

namespace Dune {
namespace GDT {


/**
 * \note source is averaged on intersections
 */
template <class E, size_t d, class R, class V, class GV, class IGV>
void raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                                  DiscreteFunction<V, GV, d, 1, R>& target,
                                  const GridView<IGV>& interpolation_grid_view,
                                  const XT::Common::Parameter& param = {})
{
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GV>>, "");
  static_assert(std::is_same_v<E, XT::Grid::extract_entity_t<GridView<IGV>>>, "");
  static_assert(d == GV::dimension, "");
  XT::Common::DefaultLogger logger("raviart_thomas_interpolation");
  LOG(debug) << "(source=" << &source << ", target=" << &target
             << ", interpolation_grid_view=" << &interpolation_grid_view << ", param=" << param << ")" << std::endl;
  if (target.space().type() != SpaceType::raviart_thomas)
    LOG(warn) << "target.space().type() is " << target.space().type() << ", not raviart_thomas! Continuing anyway ..."
              << std::endl;
  using D = typename GridView<IGV>::ctype;
  // some preparations
  const FiniteVolumeMapper<GridView<IGV>> element_mapper(interpolation_grid_view);
  auto local_target = target.local_discrete_function();
  auto local_source_element = source.local_function();
  auto local_source_neighbor = source.local_function();
  auto rt_basis = target.space().basis().localize();
  for (auto&& element : elements(interpolation_grid_view)) {
    local_target->bind(element);
    local_source_element->bind(element);
    rt_basis->bind(element);
    const auto& rt_fe = rt_basis->finite_element();
    // prepare
    const size_t sz = rt_basis->size();
    std::vector<bool> local_key_was_handled(sz, false);
    // determine the face dofs, therefore walk the intersections
    for (auto&& intersection : intersections(interpolation_grid_view, element)) {
      const auto intersection_to_local_key_map = rt_fe.coefficients().local_key_indices(1);
      const auto intersection_index = intersection.indexInInside();
      const auto& local_keys_assosiated_with_intersection = intersection_to_local_key_map[intersection_index];
      if (local_keys_assosiated_with_intersection.size() > 0) {
        const auto intersection_fe =
            make_local_orthonormal_finite_element<D, d - 1, R>(intersection.type(), rt_fe.order());
        const auto& intersection_Pk_basis = intersection_fe->basis();
        DUNE_THROW_IF(intersection_Pk_basis.size() != local_keys_assosiated_with_intersection.size(),
                      Exceptions::interpolation_error,
                      "intersection_Pk_basis.size() = " << intersection_Pk_basis.size()
                                                        << "\n   local_keys_assosiated_with_intersection.size() = "
                                                        << local_keys_assosiated_with_intersection.size());
        XT::LA::CommonDenseMatrix<R> lhs(
            local_keys_assosiated_with_intersection.size(), intersection_Pk_basis.size(), 0);
        XT::LA::CommonDenseVector<R> rhs(intersection_Pk_basis.size(), 0);
        bool there_are_intersection_dofs_to_determine = false;
        if (intersection.neighbor()) {
          const auto neighbor = intersection.outside();
          // only look at each intersection once
          if (element_mapper.global_index(element, 0) < element_mapper.global_index(neighbor, 0)) {
            there_are_intersection_dofs_to_determine = true;
            local_source_neighbor->bind(neighbor);
            // do a face quadrature, average source
            for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                     intersection.type(),
                     std::max(rt_basis->order() + intersection_Pk_basis.order(),
                              std::max(local_source_element->order(param), local_source_neighbor->order(param))
                                  + intersection_Pk_basis.order()))) {
              const auto point_on_reference_intersection = quadrature_point.position();
              const auto point_in_reference_element =
                  intersection.geometryInInside().global(point_on_reference_intersection);
              const auto point_in_reference_neighbor =
                  intersection.geometryInOutside().global(point_on_reference_intersection);
              const auto quadrature_weight = quadrature_point.weight();
              const auto normal = intersection.unitOuterNormal(point_on_reference_intersection);
              const auto integration_factor =
                  intersection.geometry().integrationElement(point_on_reference_intersection);
              const auto rt_basis_values = rt_basis->evaluate_set(point_in_reference_element);
              const auto intersection_Pk_basis_values = intersection_Pk_basis.evaluate(point_on_reference_intersection);
              const auto local_source_values = (local_source_element->evaluate(point_in_reference_element, param)
                                                + local_source_neighbor->evaluate(point_in_reference_neighbor, param))
                                               * 0.5;
              for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
                const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
                for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                  lhs.add_to_entry(ii,
                                   jj,
                                   quadrature_weight * integration_factor * (rt_basis_values[local_key_index] * normal)
                                       * intersection_Pk_basis_values[jj]);
              }
              for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                rhs[jj] += quadrature_weight * integration_factor * (local_source_values * normal)
                           * intersection_Pk_basis_values[jj];
            }
          } else {
            there_are_intersection_dofs_to_determine = false;
            // we still need to mark the local DoFs as handled, since they were determined earlier when looking at this
            // intersection from the other side
            for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
              const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
              assert(!local_key_was_handled[local_key_index]);
              local_key_was_handled[local_key_index] = true;
            }
          }
        } else {
          there_are_intersection_dofs_to_determine = true;
          // do a face quadrature
          for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                   intersection.type(),
                   std::max(rt_basis->order() + intersection_Pk_basis.order(),
                            local_source_element->order(param) + intersection_Pk_basis.order()))) {
            const auto point_on_reference_intersection = quadrature_point.position();
            const auto point_in_reference_element =
                intersection.geometryInInside().global(point_on_reference_intersection);
            const auto quadrature_weight = quadrature_point.weight();
            const auto normal = intersection.unitOuterNormal(point_on_reference_intersection);
            const auto integration_factor = intersection.geometry().integrationElement(point_on_reference_intersection);
            const auto rt_basis_values = rt_basis->evaluate_set(point_in_reference_element);
            const auto intersection_Pk_basis_values = intersection_Pk_basis.evaluate(point_on_reference_intersection);
            const auto local_source_values = local_source_element->evaluate(point_in_reference_element, param);
            for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
              const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
              for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                lhs.add_to_entry(ii,
                                 jj,
                                 quadrature_weight * integration_factor * (rt_basis_values[local_key_index] * normal)
                                     * intersection_Pk_basis_values[jj]);
            }
            for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
              rhs[jj] += quadrature_weight * integration_factor * (local_source_values * normal)
                         * intersection_Pk_basis_values[jj];
          }
        }
        if (there_are_intersection_dofs_to_determine) {
          XT::LA::CommonDenseVector<R> intersection_dofs(local_keys_assosiated_with_intersection.size(), 0);
          try {
            intersection_dofs = XT::LA::solve(lhs, rhs);
          } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
            DUNE_THROW(Exceptions::interpolation_error,
                       "Failed to solve for DoFs associated with intersection "
                           << intersection_index << ", this was the original error:\n   " << ee.what());
          }
          for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
            const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
            assert(!local_key_was_handled[local_key_index]);
            local_target->dofs()[local_key_index] = intersection_dofs[ii];
            local_key_was_handled[local_key_index] = true;
          }
        }
      }
    }
    // determine the volume dofs
    const auto element_to_local_key_map = rt_fe.coefficients().local_key_indices(0);
    const auto& local_keys_assosiated_with_element = element_to_local_key_map[0];
    if (local_keys_assosiated_with_element.size() > 0) {
      DUNE_THROW_IF(rt_basis->order(param) < 1,
                    Exceptions::interpolation_error,
                    "DoFs associated with the element only make sense for orders >= 1!");
      const auto element_fe = make_local_orthonormal_finite_element<D, d, R, d>(element.type(), rt_fe.order() - 1);
      const auto& element_Pkminus1_basis = element_fe->basis();
      DUNE_THROW_IF(element_Pkminus1_basis.size() != local_keys_assosiated_with_element.size(),
                    Exceptions::interpolation_error,
                    "element_Pkminus1_basis.size() = " << element_Pkminus1_basis.size()
                                                       << "\n   local_keys_assosiated_with_element.size() = "
                                                       << local_keys_assosiated_with_element.size());
      XT::LA::CommonDenseMatrix<R> lhs(local_keys_assosiated_with_element.size(), element_Pkminus1_basis.size(), 0);
      XT::LA::CommonDenseVector<R> rhs(element_Pkminus1_basis.size(), 0);
      // do a volume quadrature
      for (auto&& quadrature_point :
           QuadratureRules<D, d>::rule(element.type(),
                                       std::max(rt_basis->order() + element_Pkminus1_basis.order(),
                                                local_source_element->order(param) + element_Pkminus1_basis.order()))) {
        const auto point_in_reference_element = quadrature_point.position();
        const auto quadrature_weight = quadrature_point.weight();
        const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
        const auto rt_basis_values = rt_basis->evaluate_set(point_in_reference_element);
        const auto element_Pkminus1_basis_values = element_Pkminus1_basis.evaluate(point_in_reference_element);
        const auto local_source_values = local_source_element->evaluate(point_in_reference_element, param);
        for (size_t ii = 0; ii < local_keys_assosiated_with_element.size(); ++ii) {
          const size_t local_key_index = local_keys_assosiated_with_element[ii];
          for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
            lhs.add_to_entry(ii,
                             jj,
                             quadrature_weight * integration_factor
                                 * (rt_basis_values[local_key_index] * element_Pkminus1_basis_values[jj]));
        }
        for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
          rhs[jj] += quadrature_weight * integration_factor * (local_source_values * element_Pkminus1_basis_values[jj]);
      }
      XT::LA::CommonDenseVector<R> element_dofs(local_keys_assosiated_with_element.size(), 0);
      try {
        element_dofs = XT::LA::solve(lhs, rhs);
      } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
        DUNE_THROW(Exceptions::interpolation_error,
                   "Failed to solve for volume DoFs, this was the original error:\n   " << ee.what());
      }
      for (size_t ii = 0; ii < local_keys_assosiated_with_element.size(); ++ii) {
        const size_t local_key_index = local_keys_assosiated_with_element[ii];
        assert(!local_key_was_handled[local_key_index]);
        local_target->dofs()[local_key_index] = element_dofs[ii];
        local_key_was_handled[local_key_index] = true;
      }
    }
    // final checks that there are no other dofs left
    for (size_t ii = 0; ii < sz; ++ii)
      DUNE_THROW_IF(!local_key_was_handled[ii],
                    Exceptions::interpolation_error,
                    "The following DoF is neither associated with an intersection, nor with the element!"
                        << "\n   local DoF index: " << ii
                        << "\n   associated local_key: " << rt_fe.coefficients().local_key(ii));
  }
} // ... raviart_thomas_interpolation()


template <class E, size_t d, class R, class V, class GV>
void raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                                  const DiscreteFunction<V, GV, d, 1, R>& target,
                                  const XT::Common::Parameter& param = {})
{
  raviart_thomas_interpolation(source, target, target.space().grid_view(), param);
}


template <class VectorType, // <- has to be specified manually
          class GV,
          size_t d,
          class R,
          class E,
          class IGV>
DiscreteFunction<VectorType, GV, d, 1, R>
raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                             const SpaceInterface<GV, d, 1, R>& target_space,
                             const GridView<IGV>& interpolation_grid_view,
                             const XT::Common::Parameter& param = {})
{
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  DiscreteFunction<VectorType, GV, d, 1, R> target(target_space);
  raviart_thomas_interpolation(source, target, interpolation_grid_view, param);
  return target;
}

template <class GV, size_t d, class R, class E, class IGV>
auto raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                                  const SpaceInterface<GV, d, 1, R>& target_space,
                                  const GridView<IGV>& interpolation_grid_view,
                                  const XT::Common::Parameter& param = {})
{
  return raviart_thomas_interpolation<XT::LA::IstlDenseVector<R>>(source, target_space, interpolation_grid_view, param);
}


template <class VectorType, // <- has to be specified manually
          class GV,
          size_t d,
          class R,
          class E>
auto raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                                  const SpaceInterface<GV, d, 1, R>& target_space,
                                  const XT::Common::Parameter& param = {})
{
  return raviart_thomas_interpolation<VectorType>(source, target_space, target_space.grid_view(), param);
}

template <class GV, size_t d, class R, class E>
auto raviart_thomas_interpolation(const XT::Functions::GridFunctionInterface<E, d, 1, R>& source,
                                  const SpaceInterface<GV, d, 1, R>& target_space,
                                  const XT::Common::Parameter& param = {})
{
  return raviart_thomas_interpolation<XT::LA::IstlDenseVector<R>>(source, target_space, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_RAVIART_THOMAS_HH
