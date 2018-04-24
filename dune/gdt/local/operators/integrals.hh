// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH

#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/matrix.hh>

#include <dune/xt/functions/interfaces.hh>

#include <dune/xt/la/container/common/matrix/dense.hh>

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class BinaryEvaluationType,
          class TestBase,
          class AnsatzBase = TestBase,
          class Field = typename TestBase::RangeFieldType>
class LocalVolumeIntegralOperator : public LocalVolumeTwoFormInterface<TestBase, AnsatzBase, Field>
{
  static_assert(is_binary_volume_integrand<BinaryEvaluationType>::value, "");
  static_assert(std::is_same<typename TestBase::EntityType, typename AnsatzBase::EntityType>::value, "");
  static_assert(std::is_same<typename TestBase::DomainFieldType, typename AnsatzBase::DomainFieldType>::value, "");
  static_assert(TestBase::dimDomain == AnsatzBase::dimDomain, "");

  typedef LocalVolumeIntegralOperator<BinaryEvaluationType, TestBase, AnsatzBase, Field> ThisType;
  typedef LocalVolumeTwoFormInterface<TestBase, AnsatzBase, Field> BaseType;

  typedef typename TestBase::DomainFieldType D;
  static const size_t d = TestBase::dimDomain;

public:
  using typename BaseType::TestBaseType;
  using typename BaseType::AnsatzBaseType;
  using typename BaseType::FieldType;

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
  LocalVolumeIntegralOperator(ThisType&& source) = default;

  using BaseType::apply2;

  // copied from apply2 below
  // TODO: fix properly (use CommonDenseMatrix instead of DynamicMatrix everywhere in dune-gdt?)
  void apply2(const TestBaseType& test_base,
              const AnsatzBaseType& ansatz_base,
              XT::LA::CommonDenseMatrix<FieldType>& ret) const
  {
    const auto& entity = ansatz_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base, ansatz_base) + over_integrate_;
    const auto& quadrature = QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    DynamicMatrix<FieldType> integrand_eval(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, ansatz_base, xx, integrand_eval);
      // compute integral
      ret.axpy(integration_factor * quadrature_weight, integrand_eval);
    } // loop over all quadrature points
  } // ... apply2(...)

  void
  apply2(const TestBaseType& test_base, const AnsatzBaseType& ansatz_base, DynamicMatrix<FieldType>& ret) const override
  {
    const auto& entity = ansatz_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, test_base, ansatz_base) + over_integrate_;
    const auto& quadrature = QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
    // prepare storage
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    DynamicMatrix<FieldType> integrand_eval(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, ansatz_base, xx, integrand_eval);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        const auto& integrand_eval_row = integrand_eval[ii];
        auto& ret_row = ret[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          ret_row[jj] += integrand_eval_row[jj] * integration_factor * quadrature_weight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalVolumeIntegralOperator


template <class QuaternaryFaceIntegrandType,
          class TestBaseEntity,
          class Intersection,
          class AnsatzBaseEntity = TestBaseEntity,
          class TestBaseNeighbor = TestBaseEntity,
          class AnsatzBaseNeighbor = AnsatzBaseEntity,
          class Field = typename TestBaseEntity::RangeFieldType>
class LocalCouplingIntegralOperator : public LocalCouplingTwoFormInterface<TestBaseEntity,
                                                                           Intersection,
                                                                           AnsatzBaseEntity,
                                                                           TestBaseNeighbor,
                                                                           AnsatzBaseNeighbor,
                                                                           Field>
{
  static_assert(is_quaternary_face_integrand<QuaternaryFaceIntegrandType>::value, "");
  static_assert(std::is_same<typename TestBaseEntity::EntityType, typename AnsatzBaseEntity::EntityType>::value, "");
  static_assert(std::is_same<typename TestBaseEntity::EntityType, typename TestBaseNeighbor::EntityType>::value, "");
  static_assert(std::is_same<typename TestBaseEntity::EntityType, typename AnsatzBaseNeighbor::EntityType>::value, "");
  static_assert(
      std::is_same<typename TestBaseEntity::DomainFieldType, typename AnsatzBaseEntity::DomainFieldType>::value, "");
  static_assert(
      std::is_same<typename TestBaseEntity::DomainFieldType, typename TestBaseNeighbor::DomainFieldType>::value, "");
  static_assert(
      std::is_same<typename TestBaseEntity::DomainFieldType, typename AnsatzBaseNeighbor::DomainFieldType>::value, "");
  static_assert(TestBaseEntity::dimDomain == AnsatzBaseEntity::dimDomain, "");
  static_assert(TestBaseEntity::dimDomain == TestBaseNeighbor::dimDomain, "");
  static_assert(TestBaseEntity::dimDomain == AnsatzBaseNeighbor::dimDomain, "");

  typedef LocalCouplingIntegralOperator<QuaternaryFaceIntegrandType,
                                        TestBaseEntity,
                                        Intersection,
                                        AnsatzBaseEntity,
                                        TestBaseNeighbor,
                                        AnsatzBaseNeighbor,
                                        Field>
      ThisType;
  typedef LocalCouplingTwoFormInterface<TestBaseEntity,
                                        Intersection,
                                        AnsatzBaseEntity,
                                        TestBaseNeighbor,
                                        AnsatzBaseNeighbor,
                                        Field>
      BaseType;
  typedef typename TestBaseEntity::DomainFieldType D;
  static const size_t d = TestBaseEntity::dimDomain;

public:
  using typename BaseType::TestBaseEntityType;
  using typename BaseType::AnsatzBaseEntityType;
  using typename BaseType::TestBaseNeighborType;
  using typename BaseType::AnsatzBaseNeighborType;
  using typename BaseType::IntersectionType;
  using typename BaseType::FieldType;

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

  LocalCouplingIntegralOperator(const ThisType& other) = default;
  LocalCouplingIntegralOperator(ThisType&& source) = default;

  void apply2(const TestBaseEntityType& test_base_en,
              const AnsatzBaseEntityType& ansatz_base_en,
              const TestBaseNeighborType& test_base_ne,
              const AnsatzBaseNeighborType& ansatz_base_ne,
              const IntersectionType& intersection,
              DynamicMatrix<FieldType>& ret_en_en,
              DynamicMatrix<FieldType>& ret_ne_ne,
              DynamicMatrix<FieldType>& ret_en_ne,
              DynamicMatrix<FieldType>& ret_ne_en) const override
  {
    // local inducing function
    const auto& entity = test_base_en.entity();
    const auto& neighbor = test_base_ne.entity();
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
    // \todo: make mutable member, after SMP refactor
    DynamicMatrix<FieldType> integrand_eval_en_en(rows_en, cols_en, 0.);
    DynamicMatrix<FieldType> integrand_eval_ne_ne(rows_ne, cols_ne, 0.);
    DynamicMatrix<FieldType> integrand_eval_en_ne(rows_en, cols_ne, 0.);
    DynamicMatrix<FieldType> integrand_eval_ne_en(rows_ne, cols_en, 0.);
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate local
      integrand_.evaluate(local_functions_en,
                          local_functions_ne,
                          test_base_en,
                          ansatz_base_en,
                          test_base_ne,
                          ansatz_base_ne,
                          intersection,
                          xx,
                          integrand_eval_en_en,
                          integrand_eval_ne_ne,
                          integrand_eval_en_ne,
                          integrand_eval_ne_en);
      // compute integrals
      // loop over all entity test basis functions
      for (size_t ii = 0; ii < rows_en; ++ii) {
        auto& ret_en_en_row = ret_en_en[ii];
        auto& ret_en_ne_row = ret_en_ne[ii];
        const auto& integrand_eval_en_en_row = integrand_eval_en_en[ii];
        const auto& integrand_eval_en_ne_row = integrand_eval_en_ne[ii];
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < cols_en; ++jj) {
          ret_en_en_row[jj] += integrand_eval_en_en_row[jj] * integration_factor * quadrature_weight;
        }
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < cols_ne; ++jj) {
          ret_en_ne_row[jj] += integrand_eval_en_ne_row[jj] * integration_factor * quadrature_weight;
        }
      }
      // loop over all neighbor test basis functions
      for (size_t ii = 0; ii < rows_ne; ++ii) {
        auto& ret_ne_ne_row = ret_ne_ne[ii];
        auto& ret_ne_en_row = ret_ne_en[ii];
        const auto& integrand_eval_ne_ne_row = integrand_eval_ne_ne[ii];
        const auto& integrand_eval_ne_en_row = integrand_eval_ne_en[ii];
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < cols_ne; ++jj) {
          ret_ne_ne_row[jj] += integrand_eval_ne_ne_row[jj] * integration_factor * quadrature_weight;
        }
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < cols_en; ++jj) {
          ret_ne_en_row[jj] += integrand_eval_ne_en_row[jj] * integration_factor * quadrature_weight;
        }
      }
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const QuaternaryFaceIntegrandType integrand_;
  const size_t over_integrate_;
}; // class LocalCouplingIntegralOperator


template <class BinaryFaceIntegrandType,
          class TestBase,
          class Intersection,
          class AnsatzBase = TestBase,
          class Field = typename TestBase::RangeFieldType>
class LocalBoundaryIntegralOperator : public LocalBoundaryTwoFormInterface<TestBase, Intersection, AnsatzBase, Field>
{
  static_assert(is_binary_face_integrand<BinaryFaceIntegrandType>::value, "");
  static_assert(std::is_same<typename TestBase::EntityType, typename AnsatzBase::EntityType>::value, "");
  static_assert(std::is_same<typename TestBase::DomainFieldType, typename AnsatzBase::DomainFieldType>::value, "");
  static_assert(TestBase::dimDomain == AnsatzBase::dimDomain, "");

  typedef LocalBoundaryIntegralOperator<BinaryFaceIntegrandType, TestBase, Intersection, AnsatzBase, Field> ThisType;
  typedef LocalBoundaryTwoFormInterface<TestBase, Intersection, AnsatzBase, Field> BaseType;
  typedef typename TestBase::DomainFieldType D;
  static const size_t d = TestBase::dimDomain;

public:
  using typename BaseType::TestBaseType;
  using typename BaseType::AnsatzBaseType;
  using typename BaseType::IntersectionType;
  using typename BaseType::FieldType;

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

  void apply2(const TestBaseType& test_base,
              const AnsatzBaseType& ansatz_base,
              const IntersectionType& intersection,
              DynamicMatrix<FieldType>& ret) const override
  {
    // local inducing function
    const auto& entity = test_base.entity();
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
    DynamicMatrix<FieldType> integrand_eval(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      const auto integration_factor = intersection.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate local
      integrand_.evaluate(local_functions, test_base, ansatz_base, intersection, xx, integrand_eval);
      // compute integral
      assert(integrand_eval.rows() >= rows);
      assert(integrand_eval.cols() >= cols);
      // loop over all test basis functions
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& ret_row = ret[ii];
        const auto& integrand_eval_row = integrand_eval[ii];
        // loop over all ansatz basis functions
        for (size_t jj = 0; jj < cols; ++jj) {
          ret_row[jj] += integrand_eval_row[jj] * integration_factor * quadrature_weight;
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
