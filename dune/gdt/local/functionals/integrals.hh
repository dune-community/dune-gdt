// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2016, 2018)
//   René Milk       (2017)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/local/integrands/generic.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r = 1, size_t rC = 1, class R = double, class F = R>
class LocalElementIntegralFunctional : public LocalElementFunctionalInterface<E, r, rC, R, F>
{
  using ThisType = LocalElementIntegralFunctional;
  using BaseType = LocalElementFunctionalInterface<E, r, rC, R, F>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::LocalBasisType;
  using IntegrandType = LocalUnaryElementIntegrandInterface<E, r, rC, R>;
  using GenericIntegrand = GenericLocalUnaryElementIntegrand<E, r, rC, R>;

  LocalElementIntegralFunctional(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy_as_unary_element_integrand())
    , over_integrate_(over_integrate)
  {}

  LocalElementIntegralFunctional(
      typename GenericIntegrand::GenericOrderFunctionType order_function,
      typename GenericIntegrand::GenericEvaluateFunctionType evaluate_function,
      typename GenericIntegrand::GenericPostBindFunctionType post_bind_function = [](const E&) {},
      const XT::Common::ParameterType& param_type = {},
      const int over_integrate = 0)
    : BaseType(param_type)
    , integrand_(
          GenericIntegrand(order_function, evaluate_function, post_bind_function).copy_as_unary_element_integrand())
    , over_integrate_(over_integrate)
  {}

  LocalElementIntegralFunctional(const ThisType& other)
    : BaseType(other)
    , integrand_(other.integrand_->copy_as_unary_element_integrand())
    , over_integrate_(other.over_integrate_)
  {}

  LocalElementIntegralFunctional(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  void apply(const LocalBasisType& basis,
             DynamicVector<F>& result,
             const XT::Common::Parameter& param = {}) const override final
  {
    // prepare integand
    const auto& element = basis.element();
    integrand_->bind(element);
    // prepare storage
    const auto size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    result *= 0;
    // loop over all quadrature points
    const auto integrand_order = integrand_->order(basis, param) + over_integrate_;
    for (auto&& quadrature_point : QuadratureRules<D, d>::rule(element.type(), integrand_order)) {
      const auto point_in_reference_element = quadrature_point.position();
      // integration factors
      const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_->evaluate(basis, point_in_reference_element, integrand_values_, param);
      assert(integrand_values_.size() >= size && "This must not happen!");
      // compute integral
      for (size_t ii = 0; ii < size; ++ii)
        result[ii] += integrand_values_[ii] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  mutable std::unique_ptr<IntegrandType> integrand_;
  const int over_integrate_;
  mutable DynamicVector<F> integrand_values_;
}; // class LocalElementIntegralFunctional


template <class I, size_t r = 1, size_t rC = 1, class R = double, class F = R>
class LocalIntersectionIntegralFunctional : public LocalIntersectionFunctionalInterface<I, r, rC, R, F>
{
  using ThisType = LocalIntersectionIntegralFunctional;
  using BaseType = LocalIntersectionFunctionalInterface<I, r, rC, R, F>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalBasisType;
  using IntegrandType = LocalUnaryIntersectionIntegrandInterface<I, r, rC, R>;
  //  using GenericIntegrand = GenericLocalUnaryIntersectionIntegrand<I, r, rC, R>;

  LocalIntersectionIntegralFunctional(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy_as_unary_intersection_integrand())
    , over_integrate_(over_integrate)
  {}

  //  LocalIntersectionIntegralFunctional(typename GenericIntegrand::GenericOrderFunctionType order_function,
  //                                      typename GenericIntegrand::GenericEvaluateFunctionType evaluate_function,
  //                                      const XT::Common::ParameterType& param_type = {},
  //                                      const int over_integrate = 0)
  //    : BaseType(param_type)
  //    , integrand_(GenericIntegrand(order_function, evaluate_function).copy())
  //    , over_integrate_(over_integrate)
  //  {}

  LocalIntersectionIntegralFunctional(const ThisType& other)
    : BaseType(other)
    , integrand_(other.integrand_->copy_as_unary_intersection_integrand())
    , over_integrate_(other.over_integrate_)
  {}

  LocalIntersectionIntegralFunctional(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  bool inside() const override final
  {
    return integrand_->inside();
  }

  using BaseType::apply;

  void apply(const IntersectionType& intersection,
             const LocalBasisType& test_basis,
             DynamicVector<F>& result,
             const XT::Common::Parameter& param = {}) const override final
  {
    // prepare integand
    integrand_->bind(intersection);
    // prepare storage
    const size_t size = test_basis.size(param);
    if (result.size() < size)
      result.resize(size);
    result *= 0;
    // loop over all quadrature points
    const auto integrand_order = integrand_->order(test_basis, param) + over_integrate_;
    for (const auto& quadrature_point : QuadratureRules<D, d - 1>::rule(intersection.type(), integrand_order)) {
      const auto point_in_reference_intersection = quadrature_point.position();
      // integration factors
      const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_->evaluate(test_basis, point_in_reference_intersection, integrand_values_, param);
      assert(integrand_values_.size() >= size && "This must not happen!");
      // compute integral
      for (size_t ii = 0; ii < size; ++ii)
        result[ii] += integrand_values_[ii] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  mutable std::unique_ptr<IntegrandType> integrand_;
  const int over_integrate_;
  mutable DynamicVector<F> integrand_values_;
}; // class LocalIntersectionIntegralFunctional


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
