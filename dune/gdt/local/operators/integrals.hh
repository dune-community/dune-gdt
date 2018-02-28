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

#include <dune/geometry/quadraturerules.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * For an explanation of the template arguments \sa LocalElementTwoFormInterface
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalElementIntegralOperator : public LocalElementTwoFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = LocalElementIntegralOperator<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;
  using BaseType = LocalElementTwoFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::LocalTestBasisType;
  using typename BaseType::LocalAnsatzBasisType;

  using IntegrandType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

  LocalElementIntegralOperator(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy())
    , over_integrate_(over_integrate)
  {
  }

  LocalElementIntegralOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , integrand_(other.integrand_->copy())
    , over_integrate_(other.over_integrate_)
  {
  }

  LocalElementIntegralOperator(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const LocalTestBasisType& test_basis,
              const LocalAnsatzBasisType& ansatz_basis,
              DynamicMatrix<F>& result,
              const XT::Common::Parameter& param = {}) const override
  {
    // prepare integand
    const auto& element = ansatz_basis.entity();
    assert(test_basis.entity() == element && "This must not happen!");
    integrand_->bind(element);
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // loop over all quadrature points
    const size_t integrand_order = integrand_->order(test_basis, ansatz_basis) + over_integrate_;
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.type(), integrand_order)) {
      const auto point_in_reference_element = quadrature_point.position();
      // integration factors
      const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_->evaluate(test_basis, ansatz_basis, point_in_reference_element, integrand_values_, param);
      assert(integrand_values_.rows() >= rows && "This must not happen!");
      assert(integrand_values_.cols() >= cols && "This must not happen!");
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += integrand_values_[ii][jj] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  mutable std::unique_ptr<IntegrandType> integrand_;
  const int over_integrate_;
  mutable DynamicMatrix<F> integrand_values_;
}; // class LocalElementIntegralOperator


#if 0
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
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTEGRALS_HH
