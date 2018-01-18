// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_OPERATORS_ADVECTION_FV_HH

#include <functional>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/onedgrid.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/local/operators/lambda.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// forward
template <class DF, class GL>
class AdvectionFvOperator;


namespace internal {


template <class DF, class GL>
class AdvectionFvOperatorTraits
{
public:
  using derived_type = AdvectionFvOperator<DF, GL>;
  using FieldType = double;
  using JacobianType = int;
};


} // namespace internal


template <class DF, class GL = typename DF::SpaceType::GridLayerType>
class AdvectionFvOperator : public OperatorInterface<internal::AdvectionFvOperatorTraits<DF, GL>>
{
  static_assert(is_discrete_function<DF>::value, "");
  using SpaceType = typename DF::SpaceType;
  static_assert(SpaceType::polOrder == 0, "Use AdvectionDgOperator instead!");
  static_assert(XT::Grid::is_layer<GL>::value, "");
  using ThisType = AdvectionFvOperator<DF, GL>;
  using BaseType = OperatorInterface<internal::AdvectionFvOperatorTraits<DF, GL>>;
  using LocalCouplingOperatorType = LocalAdvectionFvCouplingOperator<SpaceType>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::JacobianType;
  using NumericalFluxType = typename LocalCouplingOperatorType::NumericalFluxType;
  using LocalAdvectionFvBoundaryOperatorByCustomExtrapolationType =
      LocalAdvectionFvBoundaryOperatorByCustomExtrapolation<SpaceType>;
  using LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxType =
      LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux<SpaceType>;

  AdvectionFvOperator(
      const GL& grid_layer,
      const NumericalFluxType& numerical_flx,
      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new XT::Grid::ApplyOn::NoIntersections<GL>())
    : BaseType(numerical_flx.parameter_type())
    , grid_layer_(grid_layer)
    , numerical_flux_(numerical_flx)
    , local_coupling_operator_(numerical_flux_.access())
    , periodicity_exception_(periodicity_exception)
  {
  }

  AdvectionFvOperator(
      const GL& grid_layer,
      NumericalFluxType*&& numerical_flx_ptr,
      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new XT::Grid::ApplyOn::NoIntersections<GL>())
    : BaseType(numerical_flx_ptr->parameter_type())
    , grid_layer_(grid_layer)
    , numerical_flux_(std::move(numerical_flx_ptr))
    , local_coupling_operator_(numerical_flux_.access())
    , periodicity_exception_(periodicity_exception)
  {
  }

  AdvectionFvOperator(
      const GL& grid_layer,
      const typename NumericalFluxType::FluxType& flux,
      std::function<typename NumericalFluxType::RangeType(const typename NumericalFluxType::RangeType&,
                                                          const typename NumericalFluxType::RangeType&,
                                                          const typename NumericalFluxType::DomainType&,
                                                          const XT::Common::Parameter&)> numerical_flux_lambda,
      const XT::Common::ParameterType& numerical_flux_parameter_type = {},
      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new XT::Grid::ApplyOn::NoIntersections<GL>())
    : BaseType(numerical_flux_parameter_type)
    , grid_layer_(grid_layer)
    , numerical_flux_(new NumericalLambdaFlux<DF>(flux, numerical_flux_lambda, numerical_flux_parameter_type))
    , local_coupling_operator_(numerical_flux_.access())
    , periodicity_exception_(periodicity_exception)
  {
  }

  AdvectionFvOperator(const ThisType& other) = delete;
  AdvectionFvOperator(ThisType&& source) = delete;

  const NumericalFluxType& numerical_flux() const
  {
    return numerical_flux_.access();
  }

  void append(typename LocalAdvectionFvBoundaryOperatorByCustomExtrapolationType::LambdaType boundary_treatment_lambda,
              XT::Grid::ApplyOn::WhichIntersection<GL>*&& filter,
              const XT::Common::ParameterType& param_type = {})
  {
    boundary_treatment_by_extrapolation_operators_.emplace_back(
        new LocalAdvectionFvBoundaryOperatorByCustomExtrapolationType(
            numerical_flux(), boundary_treatment_lambda, param_type),
        std::move(filter));
  }

  void
  append(typename LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxType::LambdaType boundary_numerical_flux_lambda,
         XT::Grid::ApplyOn::WhichIntersection<GL>*&& filter,
         const XT::Common::ParameterType& param_type = {})
  {
    boundary_treatment_by_boundary_flux_operators_.emplace_back(
        new LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxType(boundary_numerical_flux_lambda, param_type),
        std::move(filter));
  }

  void apply(const DF& source, DF& range, const XT::Common::Parameter& param = {}) const
  {
    if (!source.vector().valid())
      DUNE_THROW(InvalidStateException, "source contains inf or nan!");
    range.vector() *= 0.;
    if (!range.vector().valid())
      DUNE_THROW(InvalidStateException, "range contains inf or nan!");
    LocalizableOperatorBase<GL, DF, DF> walker(grid_layer_, source, range);
    walker.append(
        local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GL>(), param, "fv_coupling_inner");
    walker.append(local_coupling_operator_,
                  XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GL>() && !(*periodicity_exception_),
                  param,
                  "fv_coupling_periodic");
    for (const auto& boundary_treatment_op_and_filter : boundary_treatment_by_extrapolation_operators_) {
      const auto& op = boundary_treatment_op_and_filter.first.access();
      const auto& filter = *boundary_treatment_op_and_filter.second;
      walker.append(op, filter.copy(), param, "fv_boundary_extrapolation");
    }
    for (const auto& lambda_boundary_op_and_filter : boundary_treatment_by_boundary_flux_operators_) {
      const auto& op = lambda_boundary_op_and_filter.first.access();
      const auto& filter = *lambda_boundary_op_and_filter.second;
      walker.append(op, filter.copy(), param, "fv_boundary_custom_flux");
    }
    walker.walk();
    if (!range.vector().valid())
      DUNE_THROW(InvalidStateException, "range contains inf or nan!");
  } // ... apply(...)

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/,
                   const SourceType& /*source*/,
                   const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return FieldType();
  }

  template <class SourceType>
  JacobianType jacobian(const SourceType& /*source*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return JacobianType();
  }

  template <class SourceType>
  void
  jacobian(const SourceType& /*source*/, JacobianType& /*jac*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/,
                     SourceType& /*source*/,
                     const XT::Common::Configuration& /*opts*/,
                     const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "");
    return std::vector<std::string>();
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "");
    return XT::Common::Configuration();
  }

private:
  const GL& grid_layer_;
  const XT::Common::ConstStorageProvider<NumericalFluxType> numerical_flux_;
  const LocalCouplingOperatorType local_coupling_operator_;
  std::list<std::pair<XT::Common::ConstStorageProvider<LocalAdvectionFvBoundaryOperatorByCustomExtrapolationType>,
                      std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>>>>
      boundary_treatment_by_extrapolation_operators_;
  std::list<std::pair<XT::Common::ConstStorageProvider<LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxType>,
                      std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>>>>
      boundary_treatment_by_boundary_flux_operators_;
  std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>> periodicity_exception_;
}; // class AdvectionFvOperator


/**
 * \brief Estimates dt via [Cockburn, Coquel, LeFloch, 1995]
 */
template <class GL, class E, class D, size_t d, class R>
double estimate_dt_for_scalar_advection_equation(
    const GL& grid_layer,
    const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1, 1>& initial_data,
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1, 1>, 0, R, d, 1>&
            flux,
    const std::pair<R, R>& boundary_data_range = {std::numeric_limits<R>::max(), std::numeric_limits<R>::min()})
{
  R data_minimum = boundary_data_range.first;
  R data_maximum = boundary_data_range.second;
  for (auto&& entity : elements(grid_layer)) {
    const auto u0_local = initial_data.local_function(entity);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(entity.type(), u0_local->order())) {
      const auto value = u0_local->evaluate(quadrature_point.position());
      data_minimum = std::min(data_minimum, value[0]);
      data_maximum = std::max(data_maximum, value[0]);
    }
  }
  R max_flux_derivative = std::numeric_limits<R>::min();
  OneDGrid max_flux_grid(1, data_minimum, data_maximum);
  const auto max_flux_interval = *max_flux_grid.leafGridView().template begin<0>();
  for (const auto& quadrature_point : QuadratureRules<R, 1>::rule(max_flux_interval.type(), flux.order())) {
    const auto df = flux.partial_u({}, max_flux_interval.geometry().global(quadrature_point.position()));
    max_flux_derivative = std::max(max_flux_derivative, df.infinity_norm());
  }
  D perimeter_over_volume = std::numeric_limits<D>::min();
  for (auto&& entity : elements(grid_layer)) {
    D perimeter = 0;
    for (auto&& intersection : intersections(grid_layer, entity))
      perimeter += intersection.geometry().volume();
    perimeter_over_volume = std::max(perimeter_over_volume, perimeter / entity.geometry().volume());
  }
  const auto dt = 1. / (perimeter_over_volume * max_flux_derivative);
  return dt;
} // ... estimate_dt_for_scalar_advection_equation(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_FV_HH
