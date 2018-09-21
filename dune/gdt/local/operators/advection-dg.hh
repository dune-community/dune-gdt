// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH

#include <functional>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/parameter.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/dof-vector.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/local/numerical-fluxes/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \note See also LocalElementOperatorInterface for a description of the template arguments.
 *
 * \sa LocalElementOperatorInterface
 */
template <class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionDgVolumeOperator : public LocalElementOperatorInterface<SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionDgVolumeOperator<SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalElementOperatorInterface<SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::SourceType;
  using typename BaseType::LocalRangeType;

  using FluxType = XT::Functions::FunctionInterface<m, d, m, RF>;

  LocalAdvectionDgVolumeOperator(const FluxType& flux)
    : BaseType(flux.parameter_type())
    , flux_(flux)
  {
  }

  LocalAdvectionDgVolumeOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , flux_(other.flux_)
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const SourceType& source,
             LocalRangeType& local_range,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range.element();
    const auto& basis = local_range.basis();
    const auto local_source = source.local_discrete_function(element);
    const auto local_source_order = local_source->order(param);
    const auto local_basis_order = basis.order(param);
    const auto integrand_order = flux_.order(param) * local_source_order + std::max(local_basis_order - 1, 0);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      // prepare
      const auto point_in_reference_element = quadrature_point.position();
      const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate
      basis.jacobians(point_in_reference_element, basis_jacobians_, param);
      const auto source_value = local_source->evaluate(point_in_reference_element, param);
      const auto flux_value = flux_.evaluate(source_value, param);
      // compute
      for (size_t ii = 0; ii < basis.size(param); ++ii)
        local_range.dofs()[ii] += integration_factor * quadrature_weight * -1. * (flux_value * basis_jacobians_[ii]);
    }
  } // ... apply(...)

private:
  const FluxType& flux_;
  mutable std::vector<typename LocalRangeType::LocalBasisType::DerivativeRangeType> basis_jacobians_;
}; // class LocalAdvectionDgVolumeOperator


/**
 * \note See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 */
template <class I,
          class SV,
          class SGV,
          size_t m = 1,
          class SR = double,
          class RR = SR,
          class IRGV = SGV,
          class IRV = SV,
          class ORR = RR,
          class ORGV = IRGV,
          class ORV = IRV>
class LocalAdvectionDgCouplingOperator
    : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>
{
  using ThisType = LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SR, RR, IRGV, IRV, ORR, ORGV, ORV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::IntersectionType;
  using typename BaseType::SourceType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;

  using NumericalFluxType = NumericalFluxInterface<d, m, RR>;

  LocalAdvectionDgCouplingOperator(const NumericalFluxType& numerical_flux, bool compute_outside = true)
    : BaseType(numerical_flux.parameter_type())
    , numerical_flux_(numerical_flux.copy())
    , compute_outside_(compute_outside)
  {
  }

  LocalAdvectionDgCouplingOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_flux_(other.numerical_flux_->copy())
    , compute_outside_(other.compute_outside_)
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& local_range_outside,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& inside_element = local_range_inside.element();
    const auto& outside_element = local_range_outside.element();
    const auto& inside_basis = local_range_inside.basis();
    const auto& outside_basis = local_range_outside.basis();
    const auto u = source.local_discrete_function(inside_element);
    const auto v = source.local_discrete_function(outside_element);
    const auto integrand_order = std::max(inside_basis.order(param), outside_basis.order(param))
                                 + numerical_flux_->flux().order(param) * std::max(u->order(param), v->order(param));
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.geometry().type(), integrand_order)) {
      // prepare
      const auto point_in_reference_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(point_in_reference_intersection);
      const auto point_in_inside_reference_element =
          intersection.geometryInInside().global(point_in_reference_intersection);
      const auto point_in_outside_reference_element =
          intersection.geometryInOutside().global(point_in_reference_intersection);
      // evaluate
      inside_basis.evaluate(point_in_inside_reference_element, inside_basis_values_);
      if (compute_outside_)
        outside_basis.evaluate(point_in_outside_reference_element, outside_basis_values_);
      const auto g = numerical_flux_->apply(u->evaluate(point_in_inside_reference_element),
                                            v->evaluate(point_in_outside_reference_element),
                                            normal,
                                            param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        local_range_inside.dofs()[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
      if (compute_outside_)
        for (size_t ii = 0; ii < outside_basis.size(param); ++ii)
          local_range_outside.dofs()[ii] -= integration_factor * quadrature_weight * (g * outside_basis_values_[ii]);
    }
  } // ... apply(...)

private:
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  const bool compute_outside_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
  mutable std::vector<typename LocalOutsideRangeType::LocalBasisType::RangeType> outside_basis_values_;
}; // class LocalAdvectionDgCouplingOperator


/**
 * \note See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 */
template <class I, class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator
    : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::IntersectionType;
  using typename BaseType::SourceType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;

  using StateDomainType = FieldVector<typename SGV::ctype, SGV::dimension>;
  using StateRangeType = typename XT::Functions::RangeTypeSelector<SF, m, 1>::type;
  using LambdaType = std::function<StateRangeType(
      const StateRangeType& /*u*/, const StateDomainType& /*n*/, const XT::Common::Parameter& /*param*/)>;

  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(
      LambdaType numerical_boundary_flux,
      const int numerical_flux_order,
      const XT::Common::ParameterType& numerical_flux_param_type = {})
    : BaseType(numerical_flux_param_type)
    , numerical_boundary_flux_(numerical_boundary_flux)
    , numerical_flux_order_(numerical_flux_order)
  {
  }

  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_boundary_flux_(other.numerical_boundary_flux_)
    , numerical_flux_order_(other.numerical_flux_order_)
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range_inside.element();
    const auto& inside_basis = local_range_inside.basis();
    const auto u = source.local_discrete_function(element);
    const auto integrand_order = inside_basis.order(param) + numerical_flux_order_ * u->order(param);
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.geometry().type(), integrand_order)) {
      // prepare
      const auto point_in_reference_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(point_in_reference_intersection);
      const auto point_in_inside_reference_element =
          intersection.geometryInInside().global(point_in_reference_intersection);
      // evaluate
      inside_basis.evaluate(point_in_inside_reference_element, inside_basis_values_);
      const auto g = numerical_boundary_flux_(u->evaluate(point_in_inside_reference_element), normal, param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        local_range_inside.dofs()[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
    }
  } // ... apply(...)

private:
  const LambdaType numerical_boundary_flux_;
  const int numerical_flux_order_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
}; // class LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator


/**
 * \note See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 */
template <class I, class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator
    : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::IntersectionType;
  using typename BaseType::SourceType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;

  using D = typename IntersectionType::ctype;
  using NumericalFluxType = NumericalFluxInterface<d, m, RF>;
  using FluxType = typename NumericalFluxType::FluxType;
  using StateRangeType = typename XT::Functions::RangeTypeSelector<RF, m, 1>::type;
  using LambdaType =
      std::function<StateRangeType(const IntersectionType& /*intersection*/,
                                   const FieldVector<D, d - 1>& /*xx_in_reference_intersection_coordinates*/,
                                   const FluxType& /*flux*/,
                                   const StateRangeType& /*u*/,
                                   const XT::Common::Parameter& /*param*/)>;

  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_extrapolation_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(numerical_flux.parameter_type() + boundary_treatment_param_type)
    , numerical_flux_(numerical_flux.copy())
    , extrapolate_(boundary_extrapolation_lambda)
  {
  }

  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_flux_(other.numerical_flux_->copy())
    , extrapolate_(other.extrapolate_)
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range_inside.element();
    const auto& inside_basis = local_range_inside.basis();
    const auto local_source = source.local_discrete_function(element);
    const auto integrand_order =
        inside_basis.order(param) + numerical_flux_->flux().order(param) * local_source->order(param);
    for (const auto& quadrature_point :
         QuadratureRules<D, d - 1>::rule(intersection.geometry().type(), integrand_order)) {
      // prepare
      const auto point_in_reference_intersection = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      const auto normal = intersection.unitOuterNormal(point_in_reference_intersection);
      const auto point_in_inside_reference_element =
          intersection.geometryInInside().global(point_in_reference_intersection);
      // evaluate
      inside_basis.evaluate(point_in_inside_reference_element, inside_basis_values_);
      const auto u = local_source->evaluate(point_in_inside_reference_element);
      const auto v = extrapolate_(intersection, point_in_reference_intersection, numerical_flux_->flux(), u, param);
      const auto g = numerical_flux_->apply(u, v, normal, param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        local_range_inside.dofs()[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
    }
  } // ... apply(...)

private:
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  const LambdaType extrapolate_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
}; // class LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
