// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH

#include <cmath>
#include <functional>
#include <type_traits>

#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class SpaceType, class GridLayerType>
class LocalAdvectionDgOperator;


namespace internal {


template <class SpaceType, class GridLayerType>
class LocalAdvectionDgOperatorTraits
{
  static_assert(is_space<SpaceType>::value, "");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");

public:
  using derived_type = LocalAdvectionDgOperator<SpaceType, GridLayerType>;
};


} // namespace internal


template <class SpaceType, class GridLayerType = typename SpaceType::GridLayerType>
class LocalAdvectionDgOperator
    : public LocalOperatorInterface<internal::LocalAdvectionDgOperatorTraits<SpaceType, GridLayerType>>
{
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = SpaceType::dimRangeCols;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>;

public:
  using FluxType = XT::Functions::LocalizableFluxFunctionInterface<E, D, d, StateType, 0, R, d>;
  using NumericalFluxType = std::function<R(const typename FluxType::LocalfunctionType&,
                                            const typename StateType::RangeType&,
                                            const typename StateType::RangeType&,
                                            const typename StateType::DomainType&)>;

  LocalAdvectionDgOperator(const GridLayerType& grid_layer, const FluxType& flux, NumericalFluxType numerical_flux)
    : grid_layer_(grid_layer)
    , flux_(flux)
    , numerical_flux_(numerical_flux)
  {
  }

  template <class VectorType>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    // prepare local functions
    const auto& basis = local_range.basis();
    const auto local_flux = flux_.local_function(entity);
    const auto local_source_en = source.local_function(entity);
    // prepare local containers
    XT::LA::CommonDenseMatrix<R> local_mass_matrix(basis.size(), basis.size(), 0.);
    XT::LA::CommonDenseVector<R> local_advection_op(basis.size(), 0.);
    // create volume quadrature
    const auto volume_integrand_order =
        std::max(2 * basis.order(),
                 (local_flux->order() * local_source_en->order()) + (std::max(ssize_t(0), ssize_t(basis.order()) - 1)));
    // loop over all quadrature points
    for (const auto& quadrature_point :
         QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate everything
      const auto basis_values = basis.evaluate(xx);
      const auto basis_jacobians = basis.jacobian(xx);
      const auto source_value = local_source_en->evaluate(xx);
      const auto flux_value = local_flux->evaluate(xx, source_value);
      // compute mass matrix
      for (size_t ii = 0; ii < basis.size(); ++ii)
        for (size_t jj = 0; jj < basis.size(); ++jj)
          local_mass_matrix.add_to_entry(
              ii, jj, integration_factor * quadrature_weight * (basis_values[ii] * basis_values[jj]));
      // compute volume part of advection operator
      for (size_t jj = 0; jj < basis.size(); ++jj)
        local_advection_op[jj] += integration_factor * quadrature_weight * -1. * (basis_jacobians[jj][0] * flux_value);
    }
    // walk the neighbors
    const auto intersection_it_end = grid_layer_.iend(entity);
    for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
      const auto& intersection = *intersection_it;
      if (intersection.neighbor()) {
        const auto neighbor = intersection.outside();
        const auto local_source_ne = source.local_function(neighbor);
        const auto face_integrand_order =
            basis.order() + local_flux->order() * std::max(local_source_en->order(), local_source_ne->order());
        for (const auto& quadrature_point :
             QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
          const auto x_intersection = quadrature_point.position();
          const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
          const auto quadrature_weight = quadrature_point.weight();
          const auto normal = intersection.unitOuterNormal(x_intersection);
          const auto x_entity = intersection.geometryInInside().global(x_intersection);
          const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
          // evaluate everything
          const auto basis_values = basis.evaluate(x_entity);
          const auto g = numerical_flux_(
              *local_flux, local_source_en->evaluate(x_entity), local_source_ne->evaluate(x_neighbor), normal);
          for (size_t jj = 0; jj < basis.size(); ++jj) {
            local_advection_op[jj] += integration_factor * quadrature_weight * (g * basis_values[jj]);
          }
        }
      }
    }
    XT::LA::CommonDenseVector<R> local_DoF_vector(basis.size(), 0.);
    XT::LA::make_solver(local_mass_matrix).apply(local_advection_op, local_DoF_vector);
    for (size_t ii = 0; ii < basis.size(); ++ii)
      local_range.vector().set(ii, local_DoF_vector[ii]);
  } // ... apply(...)

private:
  const GridLayerType& grid_layer_;
  const FluxType& flux_;
  const NumericalFluxType numerical_flux_;
}; // class LocalAdvectionDgOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
