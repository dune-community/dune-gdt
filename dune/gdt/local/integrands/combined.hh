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


/// \todo add LocalUnaryElementIntegrandSum
/// \todo add operator+ to LocalUnaryElementIntegrandInterface
/// \sa LocalQuaternaryIntersectionIntegrandSum
// template <class E, size_t r, size_t rC, class R, class F>
// class LocalUnaryElementIntegrandSum
//  : public XT::Common::ParametricInterface
//  , public XT::Grid::ElementBoundObject<Element>
//{
//  static_assert(XT::Grid::is_entity<Element>::value, "");

//  using ThisType = LocalUnaryElementIntegrandInterface<Element, range_dim, range_dim_cols, RangeField, Field>;

// public:
//  using E = Element;
//  using D = typename Element::Geometry::ctype;
//  static const constexpr size_t d = E::dimension;
//  using F = Field;

//  using R = RangeField;
//  static const constexpr size_t r = range_dim;
//  static const constexpr size_t rC = range_dim_cols;

//  using typename XT::Grid::ElementBoundObject<Element>::ElementType;
//  using DomainType = FieldVector<D, d>;
//  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

//  LocalUnaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
//    : XT::Common::ParametricInterface(param_type)
//  {}

//  virtual ~LocalUnaryElementIntegrandInterface() = default;

//  virtual std::unique_ptr<ThisType> copy() const = 0;

// protected:
//  virtual void post_bind(const IntersectionType& intrsctn) = 0;

// public:

//  /**
//   * Returns the polynomial order of the integrand, given the basis.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * Computes the evaluation of this integrand at the given point for each function in the basis.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual void evaluate(const LocalBasisType& basis,
//                        const DomainType& point_in_reference_element,
//                        DynamicVector<F>& result,
//                        const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * This method is provided for convenience and should not be used within library code.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual DynamicVector<F> evaluate(const LocalBasisType& basis,
//                                    const DomainType& point_in_reference_element,
//                                    const XT::Common::Parameter& param = {}) const
//  {
//    DynamicVector<F> result(basis.size(param));
//    evaluate(basis, point_in_reference_element, result, param);
//    return result;
//  }
//}; // class LocalUnaryElementIntegrandSum


/// \todo add LocalBinaryElementIntegrandSum
/// \todo add operator+ to LocalBinaryElementIntegrandInterface
/// \sa LocalQuaternaryIntersectionIntegrandSum
// template <class E, size_t t_r, size_t t_RC, class TF, class F, size_t a_r, size_t a_rC, class AF>
// class LocalBinaryElementIntegrandSum
//  : public XT::Common::ParametricInterface
//  , public XT::Grid::ElementBoundObject<Element>
//{
//  static_assert(XT::Grid::is_entity<Element>::value, "");

//  using ThisType = LocalBinaryElementIntegrandInterface<Element,
//                                                        test_range_dim,
//                                                        test_range_dim_cols,
//                                                        TestRangeField,
//                                                        Field,
//                                                        ansatz_range_dim,
//                                                        ansatz_range_dim_cols,
//                                                        AnsatzRangeField>;

// public:
//  using E = Element;
//  using D = typename Element::Geometry::ctype;
//  static const constexpr size_t d = E::dimension;
//  using F = Field;

//  using TR = TestRangeField;
//  static const constexpr size_t t_r = test_range_dim;
//  static const constexpr size_t t_rC = test_range_dim_cols;

//  using AR = AnsatzRangeField;
//  static const constexpr size_t a_r = ansatz_range_dim;
//  static const constexpr size_t a_rC = ansatz_range_dim_cols;

//  using typename XT::Grid::ElementBoundObject<Element>::ElementType;
//  using DomainType = FieldVector<D, d>;
//  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, t_r, t_rC, TR>;
//  using LocalAnsatzBasisType = XT::Functions::ElementFunctionSetInterface<E, a_r, a_rC, AR>;

//  LocalBinaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
//    : XT::Common::ParametricInterface(param_type)
//  {}

//  virtual ~LocalBinaryElementIntegrandInterface() = default;

//  virtual std::unique_ptr<ThisType> copy() const = 0;

// protected:
//  virtual void post_bind(const IntersectionType& intrsctn) = 0;

// public:
//  /**
//   * Returns the polynomial order of the integrand, given the bases.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual int order(const LocalTestBasisType& test_basis,
//                    const LocalAnsatzBasisType& ansatz_basis,
//                    const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * Computes the evaluation of this integrand at the given point for each combination of functions from the two
//   bases.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual void evaluate(const LocalTestBasisType& test_basis,
//                        const LocalAnsatzBasisType& ansatz_basis,
//                        const DomainType& point_in_reference_element,
//                        DynamicMatrix<F>& result,
//                        const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * This method is provided for convenience and should not be used within library code.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual DynamicMatrix<F> evaluate(const LocalTestBasisType& test_basis,
//                                    const LocalAnsatzBasisType& ansatz_basis,
//                                    const DomainType& point_in_reference_element,
//                                    const XT::Common::Parameter& param = {}) const
//  {
//    DynamicMatrix<F> result(test_basis.size(param), ansatz_basis.size(param), 0);
//    evaluate(test_basis, ansatz_basis, point_in_reference_element, result, param);
//    return result;
//  }
//}; // class LocalBinaryElementIntegrandInterface


/// \todo add LocalBinaryIntersectionIntegrandSum
/// \todo add operator+ to LocalBinaryIntersectionIntegrandInterface
/// \sa LocalQuaternaryIntersectionIntegrandSum
// template <class I, size_t r, size_t rC, class R, class F>
// class LocalBinaryIntersectionIntegrandSum
//  : public XT::Common::ParametricInterface
//  , public XT::Grid::IntersectionBoundObject<Intersection>
//{
//  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

//  using ThisType =
//      LocalBinaryIntersectionIntegrandInterface<Intersection, range_dim, range_dim_cols, RangeField, Field>;

// public:
//  using typename XT::Grid::IntersectionBoundObject<Intersection>::IntersectionType;
//  using ElementType = XT::Grid::extract_inside_element_t<Intersection>;

//  using I = Intersection;
//  using E = ElementType;
//  using D = typename ElementType::Geometry::ctype;
//  static const constexpr size_t d = E::dimension;
//  using F = Field;

//  using RF = RangeField;
//  static const constexpr size_t r = range_dim;
//  static const constexpr size_t rC = range_dim_cols;

//  using DomainType = FieldVector<D, d - 1>;
//  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, RF>;

//  LocalBinaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {})
//    : XT::Common::ParametricInterface(param_type)
//  {}

//  virtual ~LocalBinaryIntersectionIntegrandInterface() = default;

//  virtual std::unique_ptr<ThisType> copy() const = 0;

// protected:
//  virtual void post_bind(const IntersectionType& intrsctn) = 0;

// public:
//  /**
//   * Returns the polynomial order of the integrand, given the bases.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual int order(const LocalBasisType& inside_basis,
//                    const LocalBasisType& outside_basis,
//                    const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * Computes the evaluation of this integrand at the given point for each combination of functions from the two
//   bases.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual void evaluate(const LocalBasisType& inside_basis,
//                        const LocalBasisType& outside_basis,
//                        const DomainType& point_in_reference_element,
//                        DynamicMatrix<F>& result,
//                        const XT::Common::Parameter& param = {}) const = 0;

//  /**
//   * This method is provided for convenience and should not be used within library code.
//   *
//   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
//   **/
//  virtual DynamicMatrix<F> evaluate(const LocalBasisType& inside_basis,
//                                    const LocalBasisType& outside_basis,
//                                    const DomainType& point_in_reference_element,
//                                    const XT::Common::Parameter& param = {}) const
//  {
//    DynamicMatrix<F> result(inside_basis.size(param), outside_basis.size(param), 0);
//    evaluate(inside_basis, inside_basis, point_in_reference_element, result, param);
//    return result;
//  }
//}; // class LocalBinaryIntersectionIntegrandSum


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
    , left_(left.copy().release())
    , right_(right.copy().release())
  {}

  LocalQuaternaryIntersectionIntegrandSum(const ThisType& other)
    : BaseType(other)
    , left_(other.left_.access().copy().release())
    , right_(other.right_.access().copy().release())
  {}

  LocalQuaternaryIntersectionIntegrandSum(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
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
