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

#include <type_traits>

#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>

#include "advection-fv.hh"

namespace Dune {
namespace GDT {


// forward
template <class SpaceType>
class LocalAdvectionDgVolumeOperator;

template <class SpaceType>
class LocalAdvectionDgCouplingOperator;

template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomExtrapolation;

template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux;


namespace internal {


template <class SpaceType>
class LocalAdvectionDgVolumeOperatorTraits
{
  static_assert(is_space<SpaceType>::value, "");

public:
  using derived_type = LocalAdvectionDgVolumeOperator<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionDgCouplingOperatorTraits
{
  static_assert(is_space<SpaceType>::value, "");

public:
  using derived_type = LocalAdvectionDgCouplingOperator<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomExtrapolationTraits
{
  static_assert(is_space<SpaceType>::value, "");

public:
  using derived_type = LocalAdvectionDgBoundaryOperatorByCustomExtrapolation<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxTraits
{
  static_assert(is_space<SpaceType>::value, "");

public:
  using derived_type = LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux<SpaceType>;
};


} // namespace internal


template <class SpaceType>
class LocalAdvectionDgVolumeOperator
    : public LocalOperatorInterface<internal::LocalAdvectionDgVolumeOperatorTraits<SpaceType>>
{
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;

  LocalAdvectionDgVolumeOperator(const NumericalFluxType& numerical_flux)
    : numerical_flux_(numerical_flux)
  {
  }

  template <class VectorType>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range) const
  {
    const auto& entity = local_range.entity();
    // prepare local functions
    const auto& basis = local_range.basis();
    const auto& flux = numerical_flux_.flux();
    const auto local_source = source.local_function(entity);
    // create volume quadrature
    const auto volume_integrand_order =
        (flux.order() * local_source->order()) + (std::max(ssize_t(0), ssize_t(basis.order()) - 1));
    // loop over all quadrature points
    for (const auto& quadrature_point :
         QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate everything
      const auto basis_jacobians = basis.jacobian(xx);
      const auto source_value = local_source->evaluate(xx);
      const auto flux_value = flux.evaluate(xx, source_value);
      for (size_t ii = 0; ii < basis.size(); ++ii)
        local_range.vector().add(ii, integration_factor * quadrature_weight * -1. * (flux_value * basis_jacobians[ii]));
    }


    //    const auto& entity = local_range_entity.entity();
    //    const auto& neighbor = local_range_neighbor.entity();
    //    const auto u_inside = source.local_discrete_function(entity);
    //    const auto u_outside = source.local_discrete_function(neighbor);
    //    const auto& basis_entity = local_range_entity.basis();
    //    const auto& basis_neighbor = local_range_neighbor.basis();
    //    const auto face_integrand_order =
    //        std::max(basis_entity.order(), basis_neighbor.order())
    //        + numerical_flux_.flux().order() * std::max(u_inside->order(), u_outside->order());
    //    for (const auto& quadrature_point :
    //         QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
    //      const auto x_intersection = quadrature_point.position();
    //      const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
    //      const auto quadrature_weight = quadrature_point.weight();
    //      const auto normal = intersection.unitOuterNormal(x_intersection);
    //      const auto x_entity = intersection.geometryInInside().global(x_intersection);
    //      const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
    //      // evaluate everything
    //      const auto basis_entity_values = basis_entity.evaluate(x_entity);
    //      const auto basis_neighbor_values = basis_neighbor.evaluate(x_neighbor);
    //      const auto g = numerical_flux_.apply(u_inside->evaluate(x_entity), u_outside->evaluate(x_neighbor), normal,
    //      mu);
    //      for (size_t ii = 0; ii < basis_entity.size(); ++ii)
    //        local_range_entity.vector().add(ii, integration_factor * quadrature_weight * (g *
    //        basis_entity_values[ii]));
    //      for (size_t ii = 0; ii < basis_neighbor.size(); ++ii)
    //        local_range_neighbor.vector().add(
    //            ii, integration_factor * quadrature_weight * -1.0 * (g * basis_neighbor_values[ii]));
    //    }
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;

  //  template <class VectorType>
  //  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
  //             LocalDiscreteFunction<SpaceType, VectorType>& local_range) const
  //  {
  //    const auto& entity = local_range.entity();
  //    // prepare local functions
  //    const auto& basis = local_range.basis();
  //    const auto local_flux = flux_.local_function(entity);
  //    const auto local_source_en = source.local_function(entity);
  //    // prepare local containers
  //    XT::LA::CommonDenseMatrix<R> local_mass_matrix(basis.size(), basis.size(), 0.);
  //    XT::LA::CommonDenseVector<R> local_advection_op(basis.size(), 0.);
  //    // create volume quadrature
  //    const auto volume_integrand_order =
  //        std::max(2 * basis.order(),
  //                 (local_flux->order() * local_source_en->order()) + (std::max(ssize_t(0), ssize_t(basis.order()) -
  //                 1)));
  //    // loop over all quadrature points
  //    for (const auto& quadrature_point :
  //         QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
  //      const auto xx = quadrature_point.position();
  //      // integration factors
  //      const auto integration_factor = entity.geometry().integrationElement(xx);
  //      const auto quadrature_weight = quadrature_point.weight();
  //      // evaluate everything
  //      const auto basis_values = basis.evaluate(xx);
  //      const auto basis_jacobians = basis.jacobian(xx);
  //      const auto source_value = local_source_en->evaluate(xx);
  //      const auto flux_value = local_flux->evaluate(xx, source_value);
  //      // compute mass matrix
  //      for (size_t ii = 0; ii < basis.size(); ++ii)
  //        for (size_t jj = 0; jj < basis.size(); ++jj)
  //          local_mass_matrix.add_to_entry(
  //              ii, jj, integration_factor * quadrature_weight * (basis_values[ii] * basis_values[jj]));
  //      // compute volume part of advection operator
  //      for (size_t jj = 0; jj < basis.size(); ++jj)
  //        local_advection_op[jj] += integration_factor * quadrature_weight * -1. * (basis_jacobians[jj][0] *
  //        flux_value);
  //    }
  //    // walk the neighbors
  //    const auto intersection_it_end = grid_layer_.iend(entity);
  //    for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
  //    ++intersection_it) {
  //      const auto& intersection = *intersection_it;
  //      if (intersection.neighbor()) {
  //        const auto neighbor = intersection.outside();
  //        const auto local_source_ne = source.local_function(neighbor);
  //        const auto face_integrand_order =
  //            basis.order() + local_flux->order() * std::max(local_source_en->order(), local_source_ne->order());
  //        for (const auto& quadrature_point :
  //             QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
  //          const auto x_intersection = quadrature_point.position();
  //          const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
  //          const auto quadrature_weight = quadrature_point.weight();
  //          const auto normal = intersection.unitOuterNormal(x_intersection);
  //          const auto x_entity = intersection.geometryInInside().global(x_intersection);
  //          const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
  //          // evaluate everything
  //          const auto basis_values = basis.evaluate(x_entity);
  //          const auto g = numerical_flux_(
  //              *local_flux, local_source_en->evaluate(x_entity), local_source_ne->evaluate(x_neighbor), normal);
  //          for (size_t jj = 0; jj < basis.size(); ++jj) {
  //            local_advection_op[jj] += integration_factor * quadrature_weight * (g * basis_values[jj]);
  //          }
  //        }
  //      }
  //    }
  //    XT::LA::CommonDenseVector<R> local_DoF_vector(basis.size(), 0.);
  //    XT::LA::make_solver(local_mass_matrix).apply(local_advection_op, local_DoF_vector);
  //    for (size_t ii = 0; ii < basis.size(); ++ii)
  //      local_range.vector().set(ii, local_DoF_vector[ii]);
  //  } // ... apply(...)

  // private:
  //  const GridLayerType& grid_layer_;
  //  const FluxType& flux_;
  //  const NumericalFluxType numerical_flux_;
}; // class LocalAdvectionDgVolumeOperator


template <class SpaceType>
class LocalAdvectionDgCouplingOperator
    : public LocalCouplingOperatorInterface<internal::LocalAdvectionDgCouplingOperatorTraits<SpaceType>>
{
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;

  LocalAdvectionDgCouplingOperator(const NumericalFluxType& numerical_flux)
    : numerical_flux_(numerical_flux)
  {
  }

  const XT::Common::ParameterType& parameter_type() const override final
  {
    return numerical_flux_.parameter_type();
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor,
             const XT::Common::Parameter& mu = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range_entity.entity();
    const auto& neighbor = local_range_neighbor.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto u_outside = source.local_discrete_function(neighbor);
    const auto& basis_entity = local_range_entity.basis();
    const auto& basis_neighbor = local_range_neighbor.basis();
    const auto face_integrand_order =
        std::max(basis_entity.order(), basis_neighbor.order())
        + numerical_flux_.flux().order() * std::max(u_inside->order(), u_outside->order());
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
      const auto x_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(x_intersection);
      const auto x_entity = intersection.geometryInInside().global(x_intersection);
      const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
      // evaluate everything
      const auto basis_entity_values = basis_entity.evaluate(x_entity);
      const auto basis_neighbor_values = basis_neighbor.evaluate(x_neighbor);
      const auto g = numerical_flux_.apply(u_inside->evaluate(x_entity), u_outside->evaluate(x_neighbor), normal, mu);
      for (size_t ii = 0; ii < basis_entity.size(); ++ii)
        local_range_entity.vector().add(ii, integration_factor * quadrature_weight * (g * basis_entity_values[ii]));
      for (size_t ii = 0; ii < basis_neighbor.size(); ++ii)
        local_range_neighbor.vector().add(
            ii, integration_factor * quadrature_weight * -1.0 * (g * basis_neighbor_values[ii]));
    }
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;

  //  template <class VectorType>
  //  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
  //             LocalDiscreteFunction<SpaceType, VectorType>& local_range) const
  //  {
  //    const auto& entity = local_range.entity();
  //    // prepare local functions
  //    const auto& basis = local_range.basis();
  //    const auto local_flux = flux_.local_function(entity);
  //    const auto local_source_en = source.local_function(entity);
  //    // prepare local containers
  //    XT::LA::CommonDenseMatrix<R> local_mass_matrix(basis.size(), basis.size(), 0.);
  //    XT::LA::CommonDenseVector<R> local_advection_op(basis.size(), 0.);
  //    // create volume quadrature
  //    const auto volume_integrand_order =
  //        std::max(2 * basis.order(),
  //                 (local_flux->order() * local_source_en->order()) + (std::max(ssize_t(0), ssize_t(basis.order()) -
  //                 1)));
  //    // loop over all quadrature points
  //    for (const auto& quadrature_point :
  //         QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
  //      const auto xx = quadrature_point.position();
  //      // integration factors
  //      const auto integration_factor = entity.geometry().integrationElement(xx);
  //      const auto quadrature_weight = quadrature_point.weight();
  //      // evaluate everything
  //      const auto basis_values = basis.evaluate(xx);
  //      const auto basis_jacobians = basis.jacobian(xx);
  //      const auto source_value = local_source_en->evaluate(xx);
  //      const auto flux_value = local_flux->evaluate(xx, source_value);
  //      // compute mass matrix
  //      for (size_t ii = 0; ii < basis.size(); ++ii)
  //        for (size_t jj = 0; jj < basis.size(); ++jj)
  //          local_mass_matrix.add_to_entry(
  //              ii, jj, integration_factor * quadrature_weight * (basis_values[ii] * basis_values[jj]));
  //      // compute volume part of advection operator
  //      for (size_t jj = 0; jj < basis.size(); ++jj)
  //        local_advection_op[jj] += integration_factor * quadrature_weight * -1. * (basis_jacobians[jj][0] *
  //        flux_value);
  //    }
  //    // walk the neighbors
  //    const auto intersection_it_end = grid_layer_.iend(entity);
  //    for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
  //    ++intersection_it) {
  //      const auto& intersection = *intersection_it;
  //      if (intersection.neighbor()) {
  //        const auto neighbor = intersection.outside();
  //        const auto local_source_ne = source.local_function(neighbor);
  //        const auto face_integrand_order =
  //            basis.order() + local_flux->order() * std::max(local_source_en->order(), local_source_ne->order());
  //        for (const auto& quadrature_point :
  //             QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
  //          const auto x_intersection = quadrature_point.position();
  //          const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
  //          const auto quadrature_weight = quadrature_point.weight();
  //          const auto normal = intersection.unitOuterNormal(x_intersection);
  //          const auto x_entity = intersection.geometryInInside().global(x_intersection);
  //          const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
  //          // evaluate everything
  //          const auto basis_values = basis.evaluate(x_entity);
  //          const auto g = numerical_flux_(
  //              *local_flux, local_source_en->evaluate(x_entity), local_source_ne->evaluate(x_neighbor), normal);
  //          for (size_t jj = 0; jj < basis.size(); ++jj) {
  //            local_advection_op[jj] += integration_factor * quadrature_weight * (g * basis_values[jj]);
  //          }
  //        }
  //      }
  //    }
  //    XT::LA::CommonDenseVector<R> local_DoF_vector(basis.size(), 0.);
  //    XT::LA::make_solver(local_mass_matrix).apply(local_advection_op, local_DoF_vector);
  //    for (size_t ii = 0; ii < basis.size(); ++ii)
  //      local_range.vector().set(ii, local_DoF_vector[ii]);
  //  } // ... apply(...)

  // private:
  //  const GridLayerType& grid_layer_;
  //  const FluxType& flux_;
  //  const NumericalFluxType numerical_flux_;
}; // class LocalAdvectionDgCouplingOperator


template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomExtrapolation
    : public LocalBoundaryOperatorInterface<internal::
                                                LocalAdvectionDgBoundaryOperatorByCustomExtrapolationTraits<SpaceType>>
{
  using BaseType =
      LocalBoundaryOperatorInterface<internal::LocalAdvectionDgBoundaryOperatorByCustomExtrapolationTraits<SpaceType>>;
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using IntersectionType = XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>;
  using IntersectionPointType = FieldVector<D, d - 1>;
  using RangeType = typename StateType::RangeType;
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;
  using FluxType = typename NumericalFluxType::FluxType;
  using LambdaType = std::function<RangeType(const IntersectionType&,
                                             const IntersectionPointType&,
                                             const FluxType&,
                                             const RangeType&,
                                             const XT::Common::Parameter&)>;

  LocalAdvectionDgBoundaryOperatorByCustomExtrapolation(
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_treatment_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(boundary_treatment_param_type)
    , numerical_flux_(numerical_flux)
    , boundary_treatment_(boundary_treatment_lambda)
  {
    if (numerical_flux_.parameter_type() != boundary_treatment_param_type)
      DUNE_THROW(NotImplemented, "todo: merge boundary_treatment_param_type and numerical_flux.parameter_type()!");
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range,
             const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range.entity();
    const auto& basis = local_range.basis();
    const auto u_inside = source.local_discrete_function(entity);
    const auto face_integrand_order = basis.order() + numerical_flux_.flux().order() * u_inside->order();
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
      const auto x_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(x_intersection);
      const auto x_entity = intersection.geometryInInside().global(x_intersection);
      const auto basis_values = basis.evaluate(x_entity);
      const auto u = u_inside->evaluate(x_entity);
      const auto v = boundary_treatment_(intersection, x_intersection, numerical_flux_.flux(), u, param);
      const auto g = numerical_flux_.apply(u, v, normal, param);
      for (size_t ii = 0; ii < basis.size(); ++ii)
        local_range.vector().add(ii, integration_factor * quadrature_weight * (g * basis_values[ii]));
    }
    //    const auto x_intersection = ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0);
    //    const auto normal = intersection.unitOuterNormal(x_intersection);
    //    // copy the local DoF vector to matching FieldVector
    //    typename StateType::RangeType u;
    //    if (u.size() != u_inside->vector().size())
    //      DUNE_THROW(Exceptions::local_operator_error,
    //                 "u.size() = " << u.size() << "\n   u_inside->vector().size() = " << u_inside->vector().size());
    //    for (size_t ii = 0; ii < u.size(); ++ii)
    //      u[ii] = u_inside->vector().get(ii);
    //    const auto v = boundary_treatment_(intersection, x_intersection, numerical_flux_.flux(), u, param);
    //    const auto g = numerical_flux_.apply(u, v, normal, /*mu*/ {});
    //    const auto h = entity.geometry().volume();
    //    for (size_t ii = 0; ii < r; ++ii)
    //      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;
  const LambdaType boundary_treatment_;
}; // class LocalAdvectionDgBoundaryOperatorByCustomExtrapolation


template <class SpaceType>
class LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux
    : public LocalBoundaryOperatorInterface<internal::
                                                LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxTraits<SpaceType>>
{
  using BaseType =
      LocalBoundaryOperatorInterface<internal::LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxTraits<SpaceType>>;
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using DomainType = typename StateType::DomainType;
  using RangeType = typename StateType::RangeType;
  using LambdaType = std::function<RangeType(const RangeType&, const DomainType&, const XT::Common::Parameter&)>;

  LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux(
      LambdaType boundary_numerical_flux_lambda,
      const int order,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(boundary_treatment_param_type)
    , boundary_numerical_flux_lambda_(boundary_numerical_flux_lambda)
    , order_(order >= 0 ? order : 0)
  {
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range,
             const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range.entity();
    const auto& basis = local_range.basis();
    const auto u_inside = source.local_discrete_function(entity);
    const auto face_integrand_order = basis.order() + order_ * u_inside->order();
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
      const auto x_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(x_intersection);
      const auto x_entity = intersection.geometryInInside().global(x_intersection);
      const auto basis_values = basis.evaluate(x_entity);
      const auto g = boundary_numerical_flux_lambda_(u_inside->evaluate(x_entity), normal, param);
      for (size_t ii = 0; ii < basis.size(); ++ii)
        local_range.vector().add(ii, integration_factor * quadrature_weight * (g * basis_values[ii]));
    }
    //    const auto x_intersection = ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0);
    //    const auto normal = intersection.unitOuterNormal(x_intersection);
    //    // copy the local DoF vector to matching FieldVector
    //    typename StateType::RangeType u;
    //    if (u.size() != u_inside->vector().size())
    //      DUNE_THROW(Exceptions::local_operator_error,
    //                 "u.size() = " << u.size() << "\n   u_inside->vector().size() = " << u_inside->vector().size());
    //    for (size_t ii = 0; ii < u.size(); ++ii)
    //      u[ii] = u_inside->vector().get(ii);
    //    const auto g = boundary_numerical_flux_lambda_(u, normal, param);
    //    const auto h = entity.geometry().volume();
    //    for (size_t ii = 0; ii < r; ++ii)
    //      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  } // ... apply(...)

private:
  const LambdaType boundary_numerical_flux_lambda_;
  const int order_;
}; // class LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
