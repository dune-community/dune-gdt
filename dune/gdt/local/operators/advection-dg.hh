// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH

#include <functional>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/intersection.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/dof-vector.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/local/numerical-fluxes/interface.hh>
#include <dune/gdt/tools/local-mass-matrix.hh>

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
  using typename BaseType::E;
  using typename BaseType::LocalRangeType;
  using typename BaseType::LocalSourceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceType;

  using FluxType = XT::Functions::FluxFunctionInterface<E, m, d, m, RF>;
  using LocalMassMatrixProviderType = LocalMassMatrixProvider<RGV, m, 1, RF>;

  LocalAdvectionDgVolumeOperator(const SourceType& source, const FluxType& flux)
    : BaseType(source, flux.parameter_type())
    , flux_(flux)
    , local_flux_(flux_.local_function())
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgVolumeOperator(const SourceType& source,
                                 const LocalMassMatrixProviderType& local_mass_matrices,
                                 const FluxType& flux)
    : BaseType(source, flux.parameter_type())
    , flux_(flux)
    , local_flux_(flux_.local_function())
    , local_mass_matrices_(local_mass_matrices)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgVolumeOperator(const SourceSpaceType& source_space,
                                 const SV& source_vector,
                                 const LocalMassMatrixProviderType& local_mass_matrices,
                                 const FluxType& flux)
    : BaseType(source_space, source_vector, flux.parameter_type())
    , flux_(flux)
    , local_flux_(flux_.local_function())
    , local_mass_matrices_(local_mass_matrices)
  {}

  LocalAdvectionDgVolumeOperator(const ThisType& other)
    : BaseType(other)
    , flux_(other.flux_)
    , local_flux_(flux_.local_function())
    , local_mass_matrices_(other.local_mass_matrices_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range.element();
    const auto& basis = local_range.basis();
    local_dofs_.resize(basis.size(param));
    local_dofs_ *= 0.;
    local_source_->bind(element);
    const auto local_source_order = local_source_->order(param);
    const auto local_basis_order = basis.order(param);
    local_flux_->bind(element);
    const auto integrand_order = local_flux_->order(param) * local_source_order + std::max(local_basis_order - 1, 0);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      // prepare
      const auto point_in_reference_element = quadrature_point.position();
      const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate
      basis.jacobians(point_in_reference_element, basis_jacobians_, param);
      const auto source_value = local_source_->evaluate(point_in_reference_element, param);
      const auto flux_value = local_flux_->evaluate(point_in_reference_element, source_value, param);
      // compute
      for (size_t ii = 0; ii < basis.size(param); ++ii)
        local_dofs_[ii] += integration_factor * quadrature_weight * -1. * (flux_value * basis_jacobians_[ii]);
    }
    // apply local mass matrix, if required (not optimal, uses a temporary)
    if (local_mass_matrices_.valid())
      local_dofs_ = local_mass_matrices_.access().local_mass_matrix_inverse(element) * local_dofs_;
    // add to local range
    for (size_t ii = 0; ii < basis.size(param); ++ii)
      local_range.dofs()[ii] += local_dofs_[ii];
  } // ... apply(...)

private:
  using BaseType::local_source_;
  const FluxType& flux_;
  std::unique_ptr<typename FluxType::LocalFunctionType> local_flux_;
  const XT::Common::ConstStorageProvider<LocalMassMatrixProviderType> local_mass_matrices_;
  mutable std::vector<typename LocalRangeType::LocalBasisType::DerivativeRangeType> basis_jacobians_;
  mutable XT::LA::CommonDenseVector<RF> local_dofs_;
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
          class IRR = SR,
          class IRGV = SGV,
          class IRV = SV,
          class ORR = IRR,
          class ORGV = IRGV,
          class ORV = IRV>
class LocalAdvectionDgCouplingOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, IRR, IRGV, IRV, ORGV, ORV>
{
  using ThisType = LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SR, IRR, IRGV, IRV, ORR, ORGV, ORV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, IRR, IRGV, IRV, ORGV, ORV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::LocalSourceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceType;

  using NumericalFluxType = NumericalFluxInterface<I, d, m, IRR>;
  using LocalMassMatrixProviderType = LocalMassMatrixProvider<IRGV, m, 1, IRR>;

  LocalAdvectionDgCouplingOperator(const SourceType& source,
                                   const NumericalFluxType& numerical_flux,
                                   bool compute_outside = true)
    : BaseType(source, numerical_flux.parameter_type())
    , local_source_outside__(source_.access().local_function())
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , compute_outside_(compute_outside)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgCouplingOperator(const SourceType& source,
                                   const LocalMassMatrixProviderType& local_mass_matrices,
                                   const NumericalFluxType& numerical_flux,
                                   bool compute_outside = true)
    : BaseType(source, numerical_flux.parameter_type())
    , local_source_outside__(source_.access().local_function())
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , compute_outside_(compute_outside)
    , local_mass_matrices_(local_mass_matrices)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgCouplingOperator(const SourceSpaceType& source_space,
                                   const SV& source_vector,
                                   const LocalMassMatrixProviderType& local_mass_matrices,
                                   const NumericalFluxType& numerical_flux,
                                   bool compute_outside = true)
    : BaseType(source_space, source_vector, numerical_flux.parameter_type())
    , local_source_outside__(source_.access().local_function())
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , compute_outside_(compute_outside)
    , local_mass_matrices_(local_mass_matrices)
  {}

  LocalAdvectionDgCouplingOperator(const ThisType& other)
    : BaseType(other)
    , local_source_outside__(source_.access().local_function())
    , numerical_flux_(other.numerical_flux_->copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , compute_outside_(other.compute_outside_)
    , local_mass_matrices_(other.local_mass_matrices_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& local_range_outside,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& inside_element = local_range_inside.element();
    const auto& outside_element = local_range_outside.element();
    const auto& inside_basis = local_range_inside.basis();
    const auto& outside_basis = local_range_outside.basis();
    inside_local_dofs_.resize(inside_basis.size(param));
    outside_local_dofs_.resize(outside_basis.size(param));
    inside_local_dofs_ *= 0.;
    outside_local_dofs_ *= 0.;
    numerical_flux_->bind(intersection);
    local_flux_->bind(intersection.inside());
    const auto inside_flux_order = local_flux_->order(param);
    local_flux_->bind(intersection.outside());
    const auto outside_flux_order = local_flux_->order(param);
    local_source_->bind(inside_element);
    local_source_outside__->bind(outside_element);
    const auto u_order = local_source_->order(param);
    const auto v_order = local_source_outside__->order(param);
    const auto integrand_order = std::max(inside_basis.order(param), outside_basis.order(param))
                                 + std::max(inside_flux_order * u_order, outside_flux_order * v_order);
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
      const auto u_val = local_source_->evaluate(point_in_inside_reference_element);
      const auto v_val = local_source_outside__->evaluate(point_in_outside_reference_element);
      const auto g = numerical_flux_->apply(point_in_reference_intersection, u_val, v_val, normal, param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        inside_local_dofs_[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
      if (compute_outside_)
        for (size_t ii = 0; ii < outside_basis.size(param); ++ii)
          outside_local_dofs_[ii] -= integration_factor * quadrature_weight * (g * outside_basis_values_[ii]);
    }
    // apply local mass matrix, if required (not optimal, uses a temporary)
    if (local_mass_matrices_.valid())
      inside_local_dofs_ = local_mass_matrices_.access().local_mass_matrix_inverse(inside_element) * inside_local_dofs_;
    if (compute_outside_ && local_mass_matrices_.valid())
      outside_local_dofs_ =
          local_mass_matrices_.access().local_mass_matrix_inverse(outside_element) * outside_local_dofs_;
    // add to local range
    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
      local_range_inside.dofs()[ii] += inside_local_dofs_[ii];
    if (compute_outside_)
      for (size_t ii = 0; ii < outside_basis.size(param); ++ii)
        local_range_outside.dofs()[ii] += outside_local_dofs_[ii];
  } // ... apply(...)

private:
  using BaseType::local_source_;
  using BaseType::source_;
  std::unique_ptr<LocalSourceType> local_source_outside__;
  const std::unique_ptr<NumericalFluxType> numerical_flux_;
  std::unique_ptr<typename NumericalFluxType::FluxType::LocalFunctionType> local_flux_;
  const bool compute_outside_;
  const XT::Common::ConstStorageProvider<LocalMassMatrixProviderType> local_mass_matrices_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
  mutable std::vector<typename LocalOutsideRangeType::LocalBasisType::RangeType> outside_basis_values_;
  mutable XT::LA::CommonDenseVector<IRR> inside_local_dofs_;
  mutable XT::LA::CommonDenseVector<ORR> outside_local_dofs_;
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
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceType;

  using StateDomainType = FieldVector<typename SGV::ctype, SGV::dimension>;
  using StateRangeType = typename XT::Functions::RangeTypeSelector<SF, m, 1>::type;
  using LambdaType = std::function<StateRangeType(
      const StateRangeType& /*u*/, const StateDomainType& /*n*/, const XT::Common::Parameter& /*param*/)>;
  using LocalMassMatrixProviderType = LocalMassMatrixProvider<RGV, m, 1, RF>;

  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(
      const SourceType& source,
      LambdaType numerical_boundary_flux,
      const int numerical_flux_order,
      const XT::Common::ParameterType& numerical_flux_param_type = {})
    : BaseType(source, numerical_flux_param_type)
    , numerical_boundary_flux_(numerical_boundary_flux)
    , numerical_flux_order_(numerical_flux_order)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(
      const SourceType& source,
      const LocalMassMatrixProviderType& local_mass_matrices,
      LambdaType numerical_boundary_flux,
      const int numerical_flux_order,
      const XT::Common::ParameterType& numerical_flux_param_type = {})
    : BaseType(source, numerical_flux_param_type)
    , numerical_boundary_flux_(numerical_boundary_flux)
    , numerical_flux_order_(numerical_flux_order)
    , local_mass_matrices_(local_mass_matrices)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(
      const SourceSpaceType& source_space,
      const SV& source_vector,
      const LocalMassMatrixProviderType& local_mass_matrices,
      LambdaType numerical_boundary_flux,
      const int numerical_flux_order,
      const XT::Common::ParameterType& numerical_flux_param_type = {})
    : BaseType(source_space, source_vector, numerical_flux_param_type)
    , numerical_boundary_flux_(numerical_boundary_flux)
    , numerical_flux_order_(numerical_flux_order)
    , local_mass_matrices_(local_mass_matrices)
  {}


  LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator(const ThisType& other)
    : BaseType(other)
    , numerical_boundary_flux_(other.numerical_boundary_flux_)
    , numerical_flux_order_(other.numerical_flux_order_)
    , local_mass_matrices_(other.local_mass_matrices_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range_inside.element();
    const auto& inside_basis = local_range_inside.basis();
    inside_local_dofs_.resize(inside_basis.size(param));
    inside_local_dofs_ *= 0.;
    local_source_->bind(element);
    const auto integrand_order = inside_basis.order(param) + numerical_flux_order_ * local_source_->order(param);
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
      const auto g =
          numerical_boundary_flux_(local_source_->evaluate(point_in_inside_reference_element), normal, param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        inside_local_dofs_[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
    }
    // apply local mass matrix, if required (not optimal, uses a temporary)
    if (local_mass_matrices_.valid())
      inside_local_dofs_ = local_mass_matrices_.access().local_mass_matrix_inverse(element) * inside_local_dofs_;
    // add to local range
    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
      local_range_inside.dofs()[ii] += inside_local_dofs_[ii];
  } // ... apply(...)

private:
  using BaseType::local_source_;
  const LambdaType numerical_boundary_flux_;
  const int numerical_flux_order_;
  const XT::Common::ConstStorageProvider<LocalMassMatrixProviderType> local_mass_matrices_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
  mutable XT::LA::CommonDenseVector<RF> inside_local_dofs_;
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
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceType;

  using D = typename IntersectionType::ctype;
  using NumericalFluxType = NumericalFluxInterface<I, d, m, RF>;
  using FluxType = typename NumericalFluxType::FluxType;
  using StateRangeType = typename XT::Functions::RangeTypeSelector<RF, m, 1>::type;
  using LambdaType =
      std::function<StateRangeType(const IntersectionType& /*intersection*/,
                                   const FieldVector<D, d - 1>& /*xx_in_reference_intersection_coordinates*/,
                                   const FluxType& /*flux*/,
                                   const StateRangeType& /*u*/,
                                   const XT::Common::Parameter& /*param*/)>;
  using LocalMassMatrixProviderType = LocalMassMatrixProvider<RGV, m, 1, RF>;

  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(
      const SourceType& source,
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_extrapolation_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(source, numerical_flux.parameter_type() + boundary_treatment_param_type)
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , extrapolate_(boundary_extrapolation_lambda)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(
      const SourceType& source,
      const LocalMassMatrixProviderType& local_mass_matrices,
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_extrapolation_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(source, numerical_flux.parameter_type() + boundary_treatment_param_type)
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , extrapolate_(boundary_extrapolation_lambda)
    , local_mass_matrices_(local_mass_matrices)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(
      const SourceSpaceType& source_space,
      const SV& source_vector,
      const LocalMassMatrixProviderType& local_mass_matrices,
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_extrapolation_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(source_space, source_vector, numerical_flux.parameter_type() + boundary_treatment_param_type)
    , numerical_flux_(numerical_flux.copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , extrapolate_(boundary_extrapolation_lambda)
    , local_mass_matrices_(local_mass_matrices)
  {}

  LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator(const ThisType& other)
    : BaseType(other)
    , numerical_flux_(other.numerical_flux_->copy())
    , local_flux_(numerical_flux_->flux().local_function())
    , extrapolate_(other.extrapolate_)
    , local_mass_matrices_(other.local_mass_matrices_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range_inside.element();
    const auto& inside_basis = local_range_inside.basis();
    inside_local_dofs_.resize(inside_basis.size(param));
    inside_local_dofs_ *= 0.;
    numerical_flux_->bind(intersection);
    local_flux_->bind(intersection.inside());
    local_source_->bind(element);
    const auto integrand_order = inside_basis.order(param) + local_flux_->order(param) * local_source_->order(param);
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
      const auto u = local_source_->evaluate(point_in_inside_reference_element);
      const auto v = extrapolate_(intersection, point_in_reference_intersection, numerical_flux_->flux(), u, param);
      const auto g = numerical_flux_->apply(point_in_reference_intersection, u, v, normal, param);
      // compute
      for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
        inside_local_dofs_[ii] += integration_factor * quadrature_weight * (g * inside_basis_values_[ii]);
    }
    // apply local mass matrix, if required (not optimal, uses a temporary)
    if (local_mass_matrices_.valid())
      inside_local_dofs_ = local_mass_matrices_.access().local_mass_matrix_inverse(element) * inside_local_dofs_;
    // add to local range
    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
      local_range_inside.dofs()[ii] += inside_local_dofs_[ii];
  } // ... apply(...)

private:
  using BaseType::local_source_;
  const std::unique_ptr<NumericalFluxType> numerical_flux_;
  std::unique_ptr<typename NumericalFluxType::FluxType::LocalFunctionType> local_flux_;
  const LambdaType extrapolate_;
  const XT::Common::ConstStorageProvider<LocalMassMatrixProviderType> local_mass_matrices_;
  mutable std::vector<typename LocalInsideRangeType::LocalBasisType::RangeType> inside_basis_values_;
  mutable XT::LA::CommonDenseVector<RF> inside_local_dofs_;
}; // class LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator


template <class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionDgArtificialViscosityShockCapturingOperator
  : public LocalElementOperatorInterface<SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionDgArtificialViscosityShockCapturingOperator<SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalElementOperatorInterface<SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::LocalRangeType;
  using typename BaseType::LocalSourceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceType;

  using FluxType = XT::Functions::FunctionInterface<m, d, m, RF>;
  using LocalMassMatrixProviderType = LocalMassMatrixProvider<RGV, m, 1, RF>;

  LocalAdvectionDgArtificialViscosityShockCapturingOperator(const SourceType& source,
                                                            const SGV& assembly_grid_view,
                                                            const double& nu_1 = 0.2,
                                                            const double& alpha_1 = 1.0,
                                                            const size_t index = 0)
    : BaseType(source)
    , local_source_outside__(source_.access().local_function())
    , assembly_grid_view_(assembly_grid_view)
    , nu_1_(nu_1)
    , alpha_1_(alpha_1)
    , index_(index)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgArtificialViscosityShockCapturingOperator(const SourceType& source,
                                                            const LocalMassMatrixProviderType& local_mass_matrices,
                                                            const SGV& assembly_grid_view,
                                                            const double& nu_1 = 0.2,
                                                            const double& alpha_1 = 1.0,
                                                            const size_t index = 0)
    : BaseType(source)
    , local_source_outside__(source_.access().local_function())
    , assembly_grid_view_(assembly_grid_view)
    , nu_1_(nu_1)
    , alpha_1_(alpha_1)
    , index_(index)
    , local_mass_matrices_(local_mass_matrices)
  {}

  /// Applies the inverse of the local mass matrix.
  LocalAdvectionDgArtificialViscosityShockCapturingOperator(


      const SourceSpaceType& source_space,
      const SV& source_vector,

      const LocalMassMatrixProviderType& local_mass_matrices,
      const SGV& assembly_grid_view,
      const double& nu_1 = 0.2,
      const double& alpha_1 = 1.0,
      const size_t index = 0)
    : BaseType(source_space, source_vector)
    , local_source_outside__(source_.access().local_function())
    , assembly_grid_view_(assembly_grid_view)
    , nu_1_(nu_1)
    , alpha_1_(alpha_1)
    , index_(index)
    , local_mass_matrices_(local_mass_matrices)
  {}

  LocalAdvectionDgArtificialViscosityShockCapturingOperator(const ThisType& other)
    : BaseType(other)
    , local_source_outside__(source_.access().local_function())
    , assembly_grid_view_(other.assembly_grid_view_)
    , nu_1_(other.nu_1_)
    , alpha_1_(other.alpha_1_)
    , index_(other.index_)
    , local_mass_matrices_(other.local_mass_matrices_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const override final
  {
    const auto& element = local_range.element();
    const auto& basis = local_range.basis();
    local_dofs_.resize(basis.size(param));
    local_dofs_ *= 0.;
    local_source_->bind(element);
    if (local_source_->order(param) <= 0 || basis.order(param) <= 0)
      return;
    // compute jump indicator (8.176)
    double element_jump_indicator = 0;
    double element_boundary_without_domain_boundary = (d == 1) ? 1. : 0.;
    for (auto&& intersection : intersections(assembly_grid_view_, element)) {
      if (intersection.neighbor() && !intersection.boundary()) {
        if (d > 1)
          element_boundary_without_domain_boundary += XT::Grid::diameter(intersection);
        const auto neighbor = intersection.outside();
        local_source_outside__->bind(neighbor);
        const auto integration_order =
            std::pow(std::max(local_source_->order(param), local_source_outside__->order(param)), 2);
        for (auto&& quadrature_point :
             QuadratureRules<D, d - 1>::rule(intersection.geometry().type(), integration_order)) {
          const auto point_in_reference_intersection = quadrature_point.position();
          const auto point_in_reference_element =
              intersection.geometryInInside().global(point_in_reference_intersection);
          const auto point_in_reference_neighbor =
              intersection.geometryInOutside().global(point_in_reference_intersection);
          const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
          const auto quadrature_weight = quadrature_point.weight();
          const auto value_on_element = local_source_->evaluate(point_in_reference_element, param)[index_];
          const auto value_on_neighbor = local_source_outside__->evaluate(point_in_reference_neighbor, param)[index_];
          element_jump_indicator +=
              integration_factor * quadrature_weight * std::pow(value_on_element - value_on_neighbor, 2);
        }
      }
      element_jump_indicator /= element_boundary_without_domain_boundary * element.geometry().volume();
    }
    // compute smoothed discrete jump indicator (8.180)
    double smoothed_discrete_jump_indicator = 0;
    const double xi_min = 0.5;
    const double xi_max = 1.5;
    if (element_jump_indicator < xi_min)
      smoothed_discrete_jump_indicator = 0;
    else if (!(element_jump_indicator < xi_max))
      smoothed_discrete_jump_indicator = 1;
    else
      smoothed_discrete_jump_indicator =
          0.5 * std::sin(M_PI * (element_jump_indicator - (xi_max - xi_min)) / (2 * (xi_max - xi_min))) + 0.5;
    // evaluate artificial viscosity form (8.183)
    if (smoothed_discrete_jump_indicator > 0) {
      const auto h = element.geometry().volume();
      for (const auto& quadrature_point : QuadratureRules<D, d>::rule(
               element.type(), std::max(0, local_source_->order(param) - 1) + std::max(0, basis.order(param) - 1))) {
        const auto point_in_reference_element = quadrature_point.position();
        const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
        const auto quadrature_weight = quadrature_point.weight();
        const auto source_jacobian = local_source_->jacobian(point_in_reference_element, param);
        basis.jacobians(point_in_reference_element, basis_jacobians_, param);
        // compute beta_h
        for (size_t ii = 0; ii < basis.size(param); ++ii)
          local_dofs_[ii] += integration_factor * quadrature_weight * nu_1_ * std::pow(h, alpha_1_)
                             * smoothed_discrete_jump_indicator * (source_jacobian[0] * basis_jacobians_[ii][0]);
      }
    }
    // apply local mass matrix, if required (not optimal, uses a temporary)
    if (local_mass_matrices_.valid())
      local_dofs_ = local_mass_matrices_.access().local_mass_matrix_inverse(element) * local_dofs_;
    // add to local range
    for (size_t ii = 0; ii < basis.size(param); ++ii)
      local_range.dofs()[ii] += local_dofs_[ii];
  } // ... apply(...)

private:
  using BaseType::local_source_;
  using BaseType::source_;
  std::unique_ptr<LocalSourceType> local_source_outside__;
  const SGV& assembly_grid_view_;
  const double nu_1_;
  const double alpha_1_;
  const size_t index_;
  const XT::Common::ConstStorageProvider<LocalMassMatrixProviderType> local_mass_matrices_;
  mutable std::vector<typename LocalRangeType::LocalBasisType::DerivativeRangeType> basis_jacobians_;
  mutable XT::LA::CommonDenseVector<RF> local_dofs_;
}; // class LocalAdvectionDgArtificialViscosityShockCapturingOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_DG_HH
