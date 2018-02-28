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

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class UnaryEvaluationType, class TestBase, class Field = typename TestBase::RangeFieldType>
class LocalVolumeIntegralFunctional : public LocalVolumeFunctionalInterface<TestBase, Field>
{
  static_assert(is_unary_volume_integrand<UnaryEvaluationType>::value, "");
  typedef LocalVolumeIntegralFunctional<UnaryEvaluationType, TestBase, Field> ThisType;
  typedef LocalVolumeFunctionalInterface<TestBase, Field> BaseType;

  typedef typename TestBase::DomainFieldType D;
  static const size_t d = TestBase::dimDomain;

public:
  using typename BaseType::TestBaseType;
  using typename BaseType::FieldType;

  template <class... Args>
  explicit LocalVolumeIntegralFunctional(Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  explicit LocalVolumeIntegralFunctional(const int over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(boost::numeric_cast<size_t>(over_integrate))
  {
  }

  template <class... Args>
  explicit LocalVolumeIntegralFunctional(const size_t over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  LocalVolumeIntegralFunctional(const ThisType& other) = default;
  LocalVolumeIntegralFunctional(ThisType&& source) = default;

  using BaseType::apply;

  void apply(const TestBaseType& test_base, DynamicVector<FieldType>& ret) const override final
  {
    const auto& entity = test_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base) + over_integrate_;
    const auto& quadrature = QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t size = test_base.size();
    ret *= 0.0;
    assert(ret.size() >= size);
    DynamicVector<FieldType> evaluation_result(size, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, xx, evaluation_result);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii)
        ret[ii] += evaluation_result[ii] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const UnaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalVolumeIntegralFunctional


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
