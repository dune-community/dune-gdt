// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r = 1, size_t rC = 1, class R = double, class F = R>
class LocalElementIntegralFunctional : public LocalElementFunctionalInterface<E, r, rC, R, F>
{
  using ThisType = LocalElementIntegralFunctional<E, r, rC, R, F>;
  using BaseType = LocalElementFunctionalInterface<E, r, rC, R, F>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::LocalBasisType;
  using IntegrandType = LocalUnaryElementIntegrandInterface<E, r, rC, R>;

  LocalElementIntegralFunctional(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy())
    , over_integrate_(over_integrate)
  {
  }

  LocalElementIntegralFunctional(const ThisType& other)
    : BaseType(other.parameter_type())
    , integrand_(other.integrand_->copy())
    , over_integrate_(other.over_integrate_)
  {
  }

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
    const auto& element = basis.entity();
    integrand_->bind(element);
    // prepare storage
    const auto size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    result *= 0;
    // loop over all quadrature points
    const auto integration_order = integrand_->order(basis, param) + over_integrate_;
    for (auto&& quadrature_point : QuadratureRules<D, d>::rule(element.type(), integration_order)) {
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


#if 0
template <class UnaryEvaluationType,
          class TestBase,
          class Intersection,
          class Field = typename TestBase::RangeFieldType>
class LocalFaceIntegralFunctional : public LocalFaceFunctionalInterface<TestBase, Intersection, Field>
{
  static_assert(is_unary_face_integrand<UnaryEvaluationType>::value, "");
  typedef LocalFaceIntegralFunctional<UnaryEvaluationType, TestBase, Intersection, Field> ThisType;
  typedef LocalFaceFunctionalInterface<TestBase, Intersection, Field> BaseType;

  typedef typename TestBase::DomainFieldType D;
  static const size_t d = TestBase::dimDomain;

public:
  using typename BaseType::TestBaseType;
  using typename BaseType::IntersectionType;
  using typename BaseType::FieldType;

  template <class... Args>
  explicit LocalFaceIntegralFunctional(Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  explicit LocalFaceIntegralFunctional(const int over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(boost::numeric_cast<size_t>(over_integrate))
  {
  }

  template <class... Args>
  explicit LocalFaceIntegralFunctional(const size_t over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  LocalFaceIntegralFunctional(const ThisType& other) = default;
  LocalFaceIntegralFunctional(ThisType&& source) = default;

  using BaseType::apply;

  void apply(const TestBaseType& test_base,
             const IntersectionType& intersection,
             DynamicVector<FieldType>& ret) const override final
  {
    const auto& entity = test_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base) + over_integrate_;
    const auto& quadrature =
        QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t size = test_base.size();
    ret *= 0.0;
    assert(ret.size() >= size);
    DynamicVector<FieldType> evaluation_result(size, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = intersection.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, intersection, xx, evaluation_result);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii)
        ret[ii] += evaluation_result[ii] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const UnaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalFaceIntegralFunctional
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
