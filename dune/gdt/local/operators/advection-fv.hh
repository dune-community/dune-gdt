// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH

#include <functional>

#include <dune/xt/common/parameter.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/dof-vector.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/local/numerical-fluxes/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \note Presumes that the basis evaluates to 1.
 *
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
class LocalAdvectionFvCouplingOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>
{
  using ThisType = LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SR, RR, IRGV, IRV, ORR, ORGV, ORV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>;

public:
  using BaseType::d;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceType;
  using StateType = typename XT::Functions::RangeTypeSelector<SR, m, 1>::type;
  using NumericalFluxType = NumericalFluxInterface<I, d, m, RR>;
  using LocalIntersectionCoords = typename NumericalFluxType::LocalIntersectionCoords;

  LocalAdvectionFvCouplingOperator(const NumericalFluxType& numerical_flux)
    : BaseType(numerical_flux.parameter_type())
    , numerical_flux_(numerical_flux.copy())
  {}

  LocalAdvectionFvCouplingOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_flux_(other.numerical_flux_->copy())
  {}

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
    DUNE_THROW_IF((local_range_inside.space().type() != SpaceType::finite_volume)
                      || (local_range_outside.space().type() != SpaceType::finite_volume),
                  Exceptions::operator_error,
                  "Use LocalAdvectionDgCouplingOperator instead!");
    evaluate_inside_state(source, intersection, u_);
    evaluate_outside_state(source, intersection, v_);
    const auto normal = intersection.centerUnitOuterNormal();
    numerical_flux_->bind(intersection);
    if (numerical_flux_->x_dependent())
      x_in_intersection_coords_ = intersection.geometry().local(intersection.geometry().center());
    const auto g = numerical_flux_->apply(x_in_intersection_coords_, u_, v_, normal, param);
    const auto h_intersection = intersection.geometry().volume();
    const auto h_inside_element = intersection.inside().geometry().volume();
    const auto h_outside_element = intersection.outside().geometry().volume();
    for (size_t ii = 0; ii < m; ++ii) {
      local_range_inside.dofs()[ii] += (g[ii] * h_intersection) / h_inside_element;
      local_range_outside.dofs()[ii] -= (g[ii] * h_intersection) / h_outside_element;
    }
  } // ... apply(...)

  static void evaluate_inside_state(const SourceType& source, const IntersectionType& intersection, StateType& u)
  {
    const auto local_source = source.local_discrete_function(intersection.inside());
    if (source.space().type() != SpaceType::finite_volume) {
      mutex_.lock();
      u = local_source->evaluate(intersection.geometryInInside().center());
      mutex_.unlock();
    } else {
      for (size_t ii = 0; ii < m; ++ii)
        u[ii] = local_source->dofs()[ii];
    }
  }

  static void evaluate_outside_state(const SourceType& source, const IntersectionType& intersection, StateType& v)
  {
    const auto local_source = source.local_discrete_function(intersection.outside());
    if (source.space().type() != SpaceType::finite_volume) {
      mutex_.lock();
      v = local_source->evaluate(intersection.geometryInOutside().center());
      mutex_.unlock();
    } else {
      for (size_t ii = 0; ii < m; ++ii)
        v[ii] = local_source->dofs()[ii];
    }
  }

private:
  std::unique_ptr<NumericalFluxType> numerical_flux_;
  mutable LocalIntersectionCoords x_in_intersection_coords_;
  mutable StateType u_;
  mutable StateType v_;
  static std::mutex mutex_;
}; // class LocalAdvectionFvCouplingOperator

template <class I,
          class SV,
          class SGV,
          size_t m,
          class SR,
          class RR,
          class IRGV,
          class IRV,
          class ORR,
          class ORGV,
          class ORV>
std::mutex LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SR, RR, IRGV, IRV, ORR, ORGV, ORV>::mutex_;


template <class I, class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;
  using CouplingOperatorType = LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceType;

  using StateDomainType = FieldVector<typename SGV::ctype, SGV::dimension>;
  using StateType = typename CouplingOperatorType::StateType;
  using LambdaType = std::function<StateType(
      const StateType& /*u*/, const StateDomainType& /*n*/, const XT::Common::Parameter& /*param*/)>;

  LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator(
      LambdaType numerical_boundary_flux_lambda, const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(boundary_treatment_param_type)
    , numerical_boundary_flux_(numerical_boundary_flux_lambda)
  {}

  LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_boundary_flux_(other.numerical_boundary_flux_)
  {}

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
    DUNE_THROW_IF((source.space().type() != SpaceType::finite_volume)
                      || (local_range_inside.space().type() != SpaceType::finite_volume),
                  Exceptions::operator_error,
                  "Use LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator instead!");
    const auto& element = local_range_inside.element();
    CouplingOperatorType::evaluate_inside_state(source, intersection, u_);
    const auto normal = intersection.centerUnitOuterNormal();
    const auto g = numerical_boundary_flux_(u_, normal, param);
    const auto h_intersection = intersection.geometry().volume();
    const auto h_element = element.geometry().volume();
    for (size_t ii = 0; ii < m; ++ii)
      local_range_inside.dofs()[ii] += (g[ii] * h_intersection) / h_element;
  } // ... apply(...)

private:
  const LambdaType numerical_boundary_flux_;
  mutable StateType u_;
}; // class LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator


template <class I, class SV, class SGV, size_t m = 1, class SF = double, class RF = SF, class RGV = SGV, class RV = SV>
class LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>
{
  using ThisType = LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SF, m, 1, RF, RGV, RV>;
  using CouplingOperatorType = LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>;

public:
  using BaseType::d;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceType;

  using D = typename IntersectionType::ctype;
  using NumericalFluxType = NumericalFluxInterface<I, d, m, RF>;
  using LocalIntersectionCoords = typename NumericalFluxType::LocalIntersectionCoords;
  using FluxType = typename NumericalFluxType::FluxType;
  using StateType = typename CouplingOperatorType::StateType;
  using LambdaType = std::function<StateType(const IntersectionType& /*intersection*/,
                                             const FieldVector<D, d - 1>& /*xx_in_reference_intersection_coordinates*/,
                                             const FluxType& /*flux*/,
                                             const StateType& /*u*/,
                                             const XT::Common::Parameter& /*param*/)>;

  LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator(
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_extrapolation_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(numerical_flux.parameter_type() + boundary_treatment_param_type)
    , numerical_flux_(numerical_flux.copy())
    , extrapolate_(boundary_extrapolation_lambda)
  {}

  LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_flux_(other.numerical_flux_->copy())
    , extrapolate_(other.extrapolate_)
  {}

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
    DUNE_THROW_IF(local_range_inside.space().type() != SpaceType::finite_volume,
                  Exceptions::operator_error,
                  "Use LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator instead!");
    numerical_flux_->bind(intersection);
    if (numerical_flux_->x_dependent())
      x_in_intersection_coords_ = intersection.geometry().local(intersection.geometry().center());
    CouplingOperatorType::evaluate_inside_state(source, intersection, u_);
    v_ = extrapolate_(intersection,
                      ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0),
                      numerical_flux_->flux(),
                      u_,
                      param);
    const auto normal = intersection.centerUnitOuterNormal();
    const auto g = numerical_flux_->apply(x_in_intersection_coords_, u_, v_, normal, param);
    const auto h_intersection = intersection.geometry().volume();
    const auto h_element = intersection.inside().geometry().volume();
    for (size_t ii = 0; ii < m; ++ii)
      local_range_inside.dofs()[ii] += (g[ii] * h_intersection) / h_element;
  } // ... apply(...)

private:
  std::unique_ptr<NumericalFluxType> numerical_flux_;
  const LambdaType extrapolate_;
  mutable LocalIntersectionCoords x_in_intersection_coords_;
  mutable StateType u_;
  mutable StateType v_;
}; // class LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
