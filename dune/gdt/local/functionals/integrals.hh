// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class UnaryEvaluationImp>
class LocalVolumeIntegralFunctional;

template <class UnaryEvaluationImp>
class LocalFaceIntegralFunctional;


namespace internal {


template <class UnaryEvaluationImp>
class LocalVolumeIntegralFunctionalTraits
{
  static_assert(is_unary_volume_integrand<UnaryEvaluationImp>::value,
                "UnaryEvaluationImp has to be derived from LocalVolumeIntegrandInterface< ..., 1 >!");

public:
  typedef LocalVolumeIntegralFunctional<UnaryEvaluationImp> derived_type;
};


template <class UnaryEvaluationImp>
class LocalFaceIntegralFunctionalTraits
{
  static_assert(is_unary_face_integrand<UnaryEvaluationImp>::value,
                "UnaryEvaluationImp has to be derived from LocalFaceIntegrandInterface< ..., 1 >!");

public:
  typedef LocalFaceIntegralFunctional<UnaryEvaluationImp> derived_type;
};


} // namespace internal


template <class UnaryEvaluationType>
class LocalVolumeIntegralFunctional
    : public LocalVolumeFunctionalInterface<internal::LocalVolumeIntegralFunctionalTraits<UnaryEvaluationType>>
{
  typedef LocalVolumeIntegralFunctional<UnaryEvaluationType> ThisType;
  typedef LocalVolumeFunctionalInterface<internal::LocalVolumeIntegralFunctionalTraits<UnaryEvaluationType>> BaseType;

public:
  typedef internal::LocalVolumeIntegralFunctionalTraits<UnaryEvaluationType> Traits;

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
  LocalVolumeIntegralFunctional(ThisType&& source)     = default;

  using BaseType::apply;

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  void apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_base,
             Dune::DynamicVector<R>& ret) const
  {
    const auto& entity         = test_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base) + over_integrate_;
    const auto& quadrature = QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t size = test_base.size();
    ret *= 0.0;
    assert(ret.size() >= size);
    Dune::DynamicVector<R> evaluation_result(size, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight  = quadrature_point.weight();
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


template <class UnaryEvaluationType>
class LocalFaceIntegralFunctional
    : public LocalFaceFunctionalInterface<internal::LocalFaceIntegralFunctionalTraits<UnaryEvaluationType>>
{
  typedef LocalFaceIntegralFunctional<UnaryEvaluationType> ThisType;
  typedef LocalFaceFunctionalInterface<internal::LocalFaceIntegralFunctionalTraits<UnaryEvaluationType>> BaseType;

public:
  typedef internal::LocalFaceIntegralFunctionalTraits<UnaryEvaluationType> Traits;

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
  LocalFaceIntegralFunctional(ThisType&& source)     = default;

  using BaseType::apply;

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  void apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_base,
             const IntersectionType& intersection, Dune::DynamicVector<R>& ret) const
  {
    const auto& entity         = test_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base) + over_integrate_;
    const auto& quadrature =
        QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t size = test_base.size();
    ret *= 0.0;
    assert(ret.size() >= size);
    Dune::DynamicVector<R> evaluation_result(size, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = intersection.geometry().integrationElement(xx);
      const auto quadrature_weight  = quadrature_point.weight();
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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTEGRALS_HH
