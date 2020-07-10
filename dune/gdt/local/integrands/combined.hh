// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_COMBINED_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_COMBINED_HH

#include <dune/xt/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r, size_t rC, class R, class F>
class LocalUnaryElementIntegrandSum : public LocalUnaryElementIntegrandInterface<E, r, rC, R, F>
{
  using ThisType = LocalUnaryElementIntegrandSum;
  using BaseType = LocalUnaryElementIntegrandInterface<E, r, rC, R, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalTestBasisType;

  LocalUnaryElementIntegrandSum(const BaseType& left, const BaseType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left.copy_as_unary_element_integrand().release())
    , right_(right.copy_as_unary_element_integrand().release())
  {}

  LocalUnaryElementIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy_as_unary_element_integrand().release())
    , right_(other.right_.access().copy_as_unary_element_integrand().release())
  {}

  LocalUnaryElementIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_unary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& elmnt) override final
  {
    left_.access().bind(elmnt);
    right_.access().bind(elmnt);
  }

public:
  int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(left_.access().order(basis, param), right_.access().order(basis, param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& basis,
                const DomainType& point_in_reference_element,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Each integrand clears its storage, so we let the left one write into ...
    left_.access().evaluate(basis, point_in_reference_element, result, param);
    // ..., the right one into ...
    right_.access().evaluate(basis, point_in_reference_element, right_result_, param);
    // ... and simply add them up (cannot use += here, vectors might have different sizes).
    const size_t size = basis.size(param);
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] += right_result_[ii];
  } // ... evaluate(...)

private:
  XT::Common::StorageProvider<BaseType> left_;
  XT::Common::StorageProvider<BaseType> right_;
  mutable DynamicVector<F> right_result_;
}; // class LocalUnaryElementIntegrandSum


template <class I, size_t r, size_t rC, class R, class F>
class LocalUnaryIntersectionIntegrandSum : public LocalUnaryIntersectionIntegrandInterface<I, r, rC, R, F>
{
  using ThisType = LocalUnaryIntersectionIntegrandSum<I, r, rC, R, F>;
  using BaseType = LocalUnaryIntersectionIntegrandInterface<I, r, rC, R, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalTestBasisType;

  LocalUnaryIntersectionIntegrandSum(const BaseType& left, const BaseType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left.copy_as_unary_intersection_integrand().release())
    , right_(right.copy_as_unary_intersection_integrand().release())
  {
    DUNE_THROW_IF(left.inside() != right.inside(),
                  Exceptions::integrand_error,
                  "left.inside() = " << left.inside() << "\n   right.inside() = " << right.inside());
  }

  LocalUnaryIntersectionIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy_as_unary_intersection_integrand().release())
    , right_(other.right_.access().copy_as_unary_intersection_integrand().release())
  {}

  LocalUnaryIntersectionIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_unary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    left_.access().bind(intrsctn);
    right_.access().bind(intrsctn);
  }

public:
  bool inside() const override final
  {
    return left_.access().inside();
  }

  int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(left_.access().order(basis, param), right_.access().order(basis, param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& basis,
                const DomainType& point_in_reference_intersection,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Each integrand clears its storage, so we let the left one write into ...
    left_.access().evaluate(basis, point_in_reference_intersection, result, param);
    // ..., the right one into ...
    right_.access().evaluate(basis, point_in_reference_intersection, right_result_, param);
    // ... and simply add them up (cannot use += here, vectors might have different sizes).
    const size_t size = basis.size(param);
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] += right_result_[ii];
  } // ... evaluate(...)

private:
  XT::Common::StorageProvider<BaseType> left_;
  XT::Common::StorageProvider<BaseType> right_;
  mutable DynamicVector<F> right_result_;
}; // class LocalUnaryIntersectionIntegrandSum


template <class E, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryElementIntegrandSum : public LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using ThisType = LocalBinaryElementIntegrandSum;

public:
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  LocalBinaryElementIntegrandSum(const BaseType& left, const BaseType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left.copy_as_binary_element_integrand().release())
    , right_(right.copy_as_binary_element_integrand().release())
  {}

  LocalBinaryElementIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy_as_binary_element_integrand().release())
    , right_(other.right_.access().copy_as_binary_element_integrand().release())
  {}

  LocalBinaryElementIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& elmnt) override final
  {
    left_.access().bind(elmnt);
    right_.access().bind(elmnt);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(left_.access().order(test_basis, ansatz_basis, param),
                    right_.access().order(test_basis, ansatz_basis, param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_element,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Each integrand clears its storage, so we let the left one write into ...
    left_.access().evaluate(test_basis, ansatz_basis, point_in_reference_element, result, param);
    // ..., the right one into ...
    right_.access().evaluate(test_basis, ansatz_basis, point_in_reference_element, right_result_, param);
    // ... and simply add them up (cannot use += here, matrices might have different sizes).
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += right_result_[ii][jj];
  } // ... evaluate(...)

private:
  XT::Common::StorageProvider<BaseType> left_;
  XT::Common::StorageProvider<BaseType> right_;
  mutable DynamicMatrix<F> right_result_;
}; // class LocalBinaryElementIntegrandSum


template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryIntersectionIntegrandSum
  : public LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using ThisType = LocalBinaryIntersectionIntegrandSum<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  LocalBinaryIntersectionIntegrandSum(const BaseType& left, const BaseType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left.copy_as_binary_intersection_integrand().release())
    , right_(right.copy_as_binary_intersection_integrand().release())
  {
    DUNE_THROW_IF(left.inside() != right.inside(),
                  Exceptions::integrand_error,
                  "left.inside() = " << left.inside() << "\n   right.inside() = " << right.inside());
  }

  LocalBinaryIntersectionIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy_as_binary_intersection_integrand().release())
    , right_(other.right_.access().copy_as_binary_intersection_integrand().release())
  {}

  LocalBinaryIntersectionIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    left_.access().bind(intrsctn);
    right_.access().bind(intrsctn);
  }

public:
  bool inside() const override final
  {
    return left_.access().inside();
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(left_.access().order(test_basis, ansatz_basis, param),
                    right_.access().order(test_basis, ansatz_basis, param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Each integrand clears its storage, so we let the left one write into ...
    left_.access().evaluate(test_basis, ansatz_basis, point_in_reference_intersection, result, param);
    // ..., the right one into ...
    right_.access().evaluate(test_basis, ansatz_basis, point_in_reference_intersection, right_result_, param);
    // ... and simply add them up (cannot use += here, matrices might have different sizes).
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += right_result_[ii][jj];
  } // ... evaluate(...)

private:
  XT::Common::StorageProvider<BaseType> left_;
  XT::Common::StorageProvider<BaseType> right_;
  mutable DynamicMatrix<F> right_result_;
}; // class LocalBinaryIntersectionIntegrandSum


template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalQuaternaryIntersectionIntegrandSum
  : public LocalQuaternaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using ThisType = LocalQuaternaryIntersectionIntegrandSum<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  LocalQuaternaryIntersectionIntegrandSum(const BaseType& left, const BaseType& right)
    : BaseType(left.parameter_type() + right.parameter_type())
    , left_(left.copy_as_quaternary_intersection_integrand().release())
    , right_(right.copy_as_quaternary_intersection_integrand().release())
  {}

  LocalQuaternaryIntersectionIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy_as_quaternary_intersection_integrand().release())
    , right_(other.right_.access().copy_as_quaternary_intersection_integrand().release())
  {}

  LocalQuaternaryIntersectionIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    left_.access().bind(intrsctn);
    right_.access().bind(intrsctn);
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(
        left_.access().order(test_basis_inside, ansatz_basis_inside, test_basis_outside, ansatz_basis_outside, param),
        right_.access().order(test_basis_inside, ansatz_basis_inside, test_basis_outside, ansatz_basis_outside, param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis_inside,
                const LocalAnsatzBasisType& ansatz_basis_inside,
                const LocalTestBasisType& test_basis_outside,
                const LocalAnsatzBasisType& ansatz_basis_outside,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result_in_in,
                DynamicMatrix<F>& result_in_out,
                DynamicMatrix<F>& result_out_in,
                DynamicMatrix<F>& result_out_out,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Each integrand clears its storage, so we let the left one write into ...
    left_.access().evaluate(test_basis_inside,
                            ansatz_basis_inside,
                            test_basis_outside,
                            ansatz_basis_outside,
                            point_in_reference_intersection,
                            result_in_in,
                            result_in_out,
                            result_out_in,
                            result_out_out,
                            param);
    // ..., the right one into ...
    right_.access().evaluate(test_basis_inside,
                             ansatz_basis_inside,
                             test_basis_outside,
                             ansatz_basis_outside,
                             point_in_reference_intersection,
                             result_in_in_,
                             result_in_out_,
                             result_out_in_,
                             result_out_out_,
                             param);
    // ... and simply add them up (cannot use += here, matrices might be larger).
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    for (size_t ii = 0; ii < rows_in; ++ii)
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_in_in[ii][jj] += result_in_in_[ii][jj];
    for (size_t ii = 0; ii < rows_in; ++ii)
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_in_out[ii][jj] += result_in_out_[ii][jj];
    for (size_t ii = 0; ii < rows_out; ++ii)
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_out_in[ii][jj] += result_out_in_[ii][jj];
    for (size_t ii = 0; ii < rows_out; ++ii)
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_out_out[ii][jj] += result_out_out_[ii][jj];
  } // ... evaluate(...)

private:
  XT::Common::StorageProvider<BaseType> left_;
  XT::Common::StorageProvider<BaseType> right_;
  mutable DynamicMatrix<F> result_in_in_;
  mutable DynamicMatrix<F> result_in_out_;
  mutable DynamicMatrix<F> result_out_in_;
  mutable DynamicMatrix<F> result_out_out_;
}; // class LocalQuaternaryIntersectionIntegrandSum


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_COMBINED_HH
