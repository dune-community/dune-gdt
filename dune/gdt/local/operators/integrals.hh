// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH

#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/xt/common/matrix.hh>

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards, to be used in the traits
template <class BinaryEvaluationType>
class LocalVolumeIntegralOperator;

template <class QuaternaryFaceIntegrandTypeType>
class LocalCouplingIntegralOperator;

template <class BinaryEvaluationType>
class LocalBoundaryIntegralOperator;


namespace internal {


template <class BinaryEvaluationType>
class LocalVolumeIntegralOperatorTraits
{
  static_assert(is_binary_volume_integrand<BinaryEvaluationType>::value,
                "BinaryEvaluationType has to be derived from LocalVolumeIntegrandInterface< ..., 2 >!");

public:
  typedef LocalVolumeIntegralOperator<BinaryEvaluationType> derived_type;
};


template <class QuaternaryFaceIntegrandTypeType>
class LocalCouplingIntegralOperatorTraits
{
  static_assert(is_quaternary_face_integrand<QuaternaryFaceIntegrandTypeType>::value,
                "QuaternaryFaceIntegrandTypeType has to be derived from LocalFaceIntegrandInterface< ..., 4 >!");

public:
  typedef LocalCouplingIntegralOperator<QuaternaryFaceIntegrandTypeType> derived_type;
};


template <class BinaryEvaluationType>
class LocalBoundaryIntegralOperatorTraits
{
  static_assert(is_binary_face_integrand<BinaryEvaluationType>::value,
                "BinaryEvaluationType has to be derived from LocalFaceIntegrandInterface< ..., 2 >!");

public:
  typedef LocalBoundaryIntegralOperator<BinaryEvaluationType> derived_type;
};


} // namespace internal


template <class BinaryEvaluationType>
class LocalVolumeIntegralOperator
    : public LocalVolumeTwoFormInterface<internal::LocalVolumeIntegralOperatorTraits<BinaryEvaluationType>>
{
  typedef LocalVolumeIntegralOperator<BinaryEvaluationType> ThisType;
  typedef LocalVolumeTwoFormInterface<internal::LocalVolumeIntegralOperatorTraits<BinaryEvaluationType>> BaseType;

public:
  typedef internal::LocalVolumeIntegralOperatorTraits<BinaryEvaluationType> Traits;

  template <class... Args>
  explicit LocalVolumeIntegralOperator(Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  explicit LocalVolumeIntegralOperator(const int over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(boost::numeric_cast<size_t>(over_integrate))
  {
  }

  template <class... Args>
  explicit LocalVolumeIntegralOperator(const size_t over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  LocalVolumeIntegralOperator(const ThisType& other) = default;
  LocalVolumeIntegralOperator(ThisType&& source)     = default;

  using BaseType::apply2;

  template <class E, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void apply2(const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base,
              const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base,
              Dune::DynamicMatrix<R>& ret) const
  {
    const auto& entity         = ansatz_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, ansatz_base, test_base) + over_integrate_;
    const auto& quadrature = QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    Dune::DynamicMatrix<R> evaluation_result(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight  = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, ansatz_base, xx, evaluation_result);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        const auto& evaluation_result_row = evaluation_result[ii];
        auto& ret_row                     = ret[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          ret_row[jj] += evaluation_result_row[jj] * integration_factor * quadrature_weight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalVolumeIntegralOperator


template <class QuaternaryFaceIntegrandTypeType>
class LocalCouplingIntegralOperator
    : public LocalCouplingTwoFormInterface<internal::
                                               LocalCouplingIntegralOperatorTraits<QuaternaryFaceIntegrandTypeType>>
{
public:
  typedef internal::LocalCouplingIntegralOperatorTraits<QuaternaryFaceIntegrandTypeType> Traits;

  template <class... Args>
  explicit LocalCouplingIntegralOperator(Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  explicit LocalCouplingIntegralOperator(const int over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(boost::numeric_cast<size_t>(over_integrate))
  {
  }

  template <class... Args>
  explicit LocalCouplingIntegralOperator(const size_t over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  template <class E, class N, class IntersectionType, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA,
            size_t rCA>
  void apply2(const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base_en,
              const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base_en,
              const Stuff::LocalfunctionSetInterface<N, D, d, R, rT, rCT>& test_base_ne,
              const Stuff::LocalfunctionSetInterface<N, D, d, R, rA, rCA>& ansatz_base_ne,
              const IntersectionType& intersection, Dune::DynamicMatrix<R>& ret_en_en,
              Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
              Dune::DynamicMatrix<R>& ret_ne_en) const
  {
    // local inducing function
    const auto& entity            = test_base_en.entity();
    const auto& neighbor          = test_base_ne.entity();
    const auto local_functions_en = integrand_.localFunctions(entity);
    const auto local_functions_ne = integrand_.localFunctions(neighbor);
    // quadrature
    const size_t integrand_order =
        integrand_.order(
            local_functions_en, local_functions_ne, test_base_en, ansatz_base_en, test_base_ne, ansatz_base_ne)
        + over_integrate_;
    const auto& quadrature =
        QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(integrand_order));
    // check matrices
    ret_en_en *= 0.0;
    ret_ne_ne *= 0.0;
    ret_en_ne *= 0.0;
    ret_ne_en *= 0.0;
    const size_t rows_en = test_base_en.size();
    const size_t cols_en = ansatz_base_en.size();
    const size_t rows_ne = test_base_ne.size();
    const size_t cols_ne = ansatz_base_ne.size();
    assert(ret_en_en.rows() >= rows_en);
    assert(ret_en_en.cols() >= cols_en);
    assert(ret_ne_ne.rows() >= rows_ne);
    assert(ret_ne_ne.cols() >= cols_ne);
    assert(ret_en_ne.rows() >= rows_en);
    assert(ret_en_ne.cols() >= cols_ne);
    assert(ret_ne_en.rows() >= rows_en);
    assert(ret_ne_en.cols() >= cols_en);
    Dune::DynamicMatrix<R> evaluation_result_en_en(
        rows_en, cols_en, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> evaluation_result_ne_ne(
        rows_ne, cols_ne, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> evaluation_result_en_ne(
        rows_en, cols_ne, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> evaluation_result_ne_en(
        rows_ne, cols_en, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx                 = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(xx);
      const auto quadrature_weight  = quadrature_point.weight();
      // evaluate local
      integrand_.evaluate(local_functions_en,
                          local_functions_ne,
                          test_base_en,
                          ansatz_base_en,
                          test_base_ne,
                          ansatz_base_ne,
                          intersection,
                          xx,
                          evaluation_result_en_en,
                          evaluation_result_ne_ne,
                          evaluation_result_en_ne,
                          evaluation_result_ne_en);
      // compute integrals
      // loop over all entity test basis functions
      for (size_t ii = 0; ii < rows_en; ++ii) {
        auto& ret_en_en_row                     = ret_en_en[ii];
        auto& ret_en_ne_row                     = ret_en_ne[ii];
        const auto& evaluation_result_en_en_row = evaluation_result_en_en[ii];
        const auto& evaluation_result_en_ne_row = evaluation_result_en_ne[ii];
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < cols_en; ++jj) {
          ret_en_en_row[jj] += evaluation_result_en_en_row[jj] * integration_factor * quadrature_weight;
        }
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < cols_ne; ++jj) {
          ret_en_ne_row[jj] += evaluation_result_en_ne_row[jj] * integration_factor * quadrature_weight;
        }
      }
      // loop over all neighbor test basis functions
      for (size_t ii = 0; ii < rows_ne; ++ii) {
        auto& ret_ne_ne_row                     = ret_ne_ne[ii];
        auto& ret_ne_en_row                     = ret_ne_en[ii];
        const auto& evaluation_result_ne_ne_row = evaluation_result_ne_ne[ii];
        const auto& evaluation_result_ne_en_row = evaluation_result_ne_en[ii];
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < cols_ne; ++jj) {
          ret_ne_ne_row[jj] += evaluation_result_ne_ne_row[jj] * integration_factor * quadrature_weight;
        }
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < cols_en; ++jj) {
          ret_ne_en_row[jj] += evaluation_result_ne_en_row[jj] * integration_factor * quadrature_weight;
        }
      }
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const QuaternaryFaceIntegrandTypeType integrand_;
  const size_t over_integrate_;
}; // class LocalCouplingIntegralOperator


template <class BinaryFaceIntegrandType>
class LocalBoundaryIntegralOperator
    : public LocalBoundaryTwoFormInterface<internal::LocalBoundaryIntegralOperatorTraits<BinaryFaceIntegrandType>>
{
public:
  typedef internal::LocalBoundaryIntegralOperatorTraits<BinaryFaceIntegrandType> Traits;

  template <class... Args>
  LocalBoundaryIntegralOperator(Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  LocalBoundaryIntegralOperator(const int over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(boost::numeric_cast<size_t>(over_integrate))
  {
  }

  template <class... Args>
  LocalBoundaryIntegralOperator(const size_t over_integrate, Args&&... args)
    : integrand_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  template <class E, class IntersectionType, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void apply2(const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& test_base,
              const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatz_base,
              const IntersectionType& intersection, Dune::DynamicMatrix<R>& ret) const
  {
    // local inducing function
    const auto& entity         = test_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const auto integrand_order = integrand_.order(local_functions, test_base, ansatz_base) + over_integrate_;
    const auto& quadrature =
        QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    ret *= 0.0;
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    Dune::DynamicMatrix<R> evaluation_result(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx              = quadrature_point.position();
      const R integration_factor = intersection.geometry().integrationElement(xx);
      const R quadrature_weight  = quadrature_point.weight();
      // evaluate local
      integrand_.evaluate(local_functions, test_base, ansatz_base, intersection, xx, evaluation_result);
      // compute integral
      assert(evaluation_result.rows() >= rows);
      assert(evaluation_result.cols() >= cols);
      // loop over all test basis functions
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& ret_row                     = ret[ii];
        const auto& evaluation_result_row = evaluation_result[ii];
        // loop over all ansatz basis functions
        for (size_t jj = 0; jj < cols; ++jj) {
          ret_row[jj] += evaluation_result_row[jj] * integration_factor * quadrature_weight;
        }
      }
    }
  } // ... apply(...)

private:
  const BinaryFaceIntegrandType integrand_;
  const size_t over_integrate_;
}; // class LocalBoundaryIntegralOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH
