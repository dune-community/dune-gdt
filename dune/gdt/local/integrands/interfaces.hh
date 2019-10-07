// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2014 - 2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/bound-object.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


// forwards (required for operator+, include are below)
template <class E, size_t r, size_t rC, class R, class F>
class LocalUnaryElementIntegrandSum;

template <class E, size_t t_r, size_t t_RC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryElementIntegrandSum;

template <class I, size_t r, size_t rC, class R, class F>
class LocalBinaryIntersectionIntegrandSum;

template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalQuaternaryIntersectionIntegrandSum;


/**
 * Interface for integrands in integrals over grid elements, which depend on one argument only (usually the test basis
 * in an integral-based functional).
 *
 * \note Regarding SMP: the integrand is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Element,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double,
          class Field = double>
class LocalUnaryElementIntegrandInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::ElementBoundObject<Element>
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalUnaryElementIntegrandInterface;

public:
  using E = Element;
  using D = typename Element::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using typename XT::Grid::ElementBoundObject<Element>::ElementType;
  using DomainType = FieldVector<D, d>;
  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

  LocalUnaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalUnaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Returns the polynomial order of the integrand, given the basis.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each function in the basis.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalBasisType& basis,
                        const DomainType& point_in_reference_element,
                        DynamicVector<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicVector<F> evaluate(const LocalBasisType& basis,
                                    const DomainType& point_in_reference_element,
                                    const XT::Common::Parameter& param = {}) const
  {
    DynamicVector<F> result(basis.size(param));
    evaluate(basis, point_in_reference_element, result, param);
    return result;
  }
}; // class LocalUnaryElementIntegrandInterface


/**
 * Interface for integrands in integrals over grid elements, which depend on two arguments (usually the test and ansatz
 * bases in an integral-based operator, aka two-form).
 *
 * \note Regarding SMP: the integrand is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Element,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class TestRangeField = double,
          class Field = double,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols,
          class AnsatzRangeField = TestRangeField>
class LocalBinaryElementIntegrandInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::ElementBoundObject<Element>
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalBinaryElementIntegrandInterface;

public:
  using E = Element;
  using D = typename Element::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using TR = TestRangeField;
  static const constexpr size_t t_r = test_range_dim;
  static const constexpr size_t t_rC = test_range_dim_cols;

  using AR = AnsatzRangeField;
  static const constexpr size_t a_r = ansatz_range_dim;
  static const constexpr size_t a_rC = ansatz_range_dim_cols;

  using typename XT::Grid::ElementBoundObject<Element>::ElementType;
  using DomainType = FieldVector<D, d>;
  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::ElementFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalBinaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalBinaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Returns the polynomial order of the integrand, given the bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalTestBasisType& test_basis,
                    const LocalAnsatzBasisType& ansatz_basis,
                    const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each combination of functions from the two bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalTestBasisType& test_basis,
                        const LocalAnsatzBasisType& ansatz_basis,
                        const DomainType& point_in_reference_element,
                        DynamicMatrix<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicMatrix<F> evaluate(const LocalTestBasisType& test_basis,
                                    const LocalAnsatzBasisType& ansatz_basis,
                                    const DomainType& point_in_reference_element,
                                    const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result(test_basis.size(param), ansatz_basis.size(param), 0);
    evaluate(test_basis, ansatz_basis, point_in_reference_element, result, param);
    return result;
  }
}; // class LocalBinaryElementIntegrandInterface


/**
 * Interface for integrands in integrals over grid intersections, which depend on two arguments (usually the bases on
 * the inside and outside the in an integral-based functional).
 *
 * \note Regarding SMP: the integrand is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Intersection,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double,
          class Field = double>
class LocalBinaryIntersectionIntegrandInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::IntersectionBoundObject<Intersection>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalBinaryIntersectionIntegrandInterface;

public:
  using typename XT::Grid::IntersectionBoundObject<Intersection>::IntersectionType;
  using ElementType = XT::Grid::extract_inside_element_t<Intersection>;

  using I = Intersection;
  using E = ElementType;
  using D = typename ElementType::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using RF = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using DomainType = FieldVector<D, d - 1>;
  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, RF>;

  LocalBinaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalBinaryIntersectionIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Returns the polynomial order of the integrand, given the bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalBasisType& inside_basis,
                    const LocalBasisType& outside_basis,
                    const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each combination of functions from the two bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalBasisType& inside_basis,
                        const LocalBasisType& outside_basis,
                        const DomainType& point_in_reference_element,
                        DynamicMatrix<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicMatrix<F> evaluate(const LocalBasisType& inside_basis,
                                    const LocalBasisType& outside_basis,
                                    const DomainType& point_in_reference_element,
                                    const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result(inside_basis.size(param), outside_basis.size(param), 0);
    evaluate(inside_basis, inside_basis, point_in_reference_element, result, param);
    return result;
  }
}; // class LocalBinaryIntersectionIntegrandInterface


template <class Intersection,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class TestRangeField = double,
          class Field = double,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols,
          class AnsatzRangeField = TestRangeField>
class LocalQuaternaryIntersectionIntegrandInterface
  : public XT::Common::ParametricInterface
  , public XT::Grid::IntersectionBoundObject<Intersection>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalQuaternaryIntersectionIntegrandInterface;

public:
  using typename XT::Grid::IntersectionBoundObject<Intersection>::IntersectionType;
  using ElementType = XT::Grid::extract_inside_element_t<Intersection>;

  using I = Intersection;
  using E = ElementType;
  using D = typename ElementType::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using TR = TestRangeField;
  static const constexpr size_t t_r = test_range_dim;
  static const constexpr size_t t_rC = test_range_dim_cols;

  using AR = AnsatzRangeField;
  static const constexpr size_t a_r = ansatz_range_dim;
  static const constexpr size_t a_rC = ansatz_range_dim_cols;

  using DomainType = FieldVector<D, d - 1>;
  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::ElementFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalQuaternaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalQuaternaryIntersectionIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  LocalQuaternaryIntersectionIntegrandSum<I, t_r, t_rC, TR, F, a_r, a_rC, AR> operator+(const ThisType& other) const
  {
    return LocalQuaternaryIntersectionIntegrandSum<I, t_r, t_rC, TR, F, a_r, a_rC, AR>(*this, other);
  }

  /**
   * Returns the polynomial order of the integrand, given the bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalTestBasisType& test_basis_inside,
                    const LocalAnsatzBasisType& ansatz_basis_inside,
                    const LocalTestBasisType& test_basis_outside,
                    const LocalAnsatzBasisType& ansatz_basis_outside,
                    const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each combination of functions from the bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalTestBasisType& test_basis_inside,
                        const LocalAnsatzBasisType& ansatz_basis_inside,
                        const LocalTestBasisType& test_basis_outside,
                        const LocalAnsatzBasisType& ansatz_basis_outside,
                        const DomainType& point_in_reference_intersection,
                        DynamicMatrix<F>& result_in_in,
                        DynamicMatrix<F>& result_in_out,
                        DynamicMatrix<F>& result_out_in,
                        DynamicMatrix<F>& result_out_out,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual std::array<DynamicMatrix<F>, 4> evaluate(const LocalTestBasisType& test_basis_inside,
                                                   const LocalAnsatzBasisType& ansatz_basis_inside,
                                                   const LocalTestBasisType& test_basis_outside,
                                                   const LocalAnsatzBasisType& ansatz_basis_outside,
                                                   const DomainType& point_in_reference_intersection,
                                                   const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result_in_in(test_basis_inside.size(param), ansatz_basis_inside.size(param), 0);
    DynamicMatrix<F> result_in_out(test_basis_inside.size(param), ansatz_basis_outside.size(param), 0);
    DynamicMatrix<F> result_out_in(test_basis_outside.size(param), ansatz_basis_inside.size(param), 0);
    DynamicMatrix<F> result_out_out(test_basis_outside.size(param), ansatz_basis_outside.size(param), 0);
    this->evaluate(test_basis_inside,
                   ansatz_basis_inside,
                   test_basis_outside,
                   ansatz_basis_outside,
                   point_in_reference_intersection,
                   result_in_in,
                   result_in_out,
                   result_out_in,
                   result_out_out,
                   param);
    return {result_in_in, result_in_out, result_out_in, result_out_out};
  } // ... apply(...)

protected:
  void ensure_size_and_clear_results(const LocalTestBasisType& test_basis_inside,
                                     const LocalAnsatzBasisType& ansatz_basis_inside,
                                     const LocalTestBasisType& test_basis_outside,
                                     const LocalAnsatzBasisType& ansatz_basis_outside,
                                     DynamicMatrix<F>& result_in_in,
                                     DynamicMatrix<F>& result_in_out,
                                     DynamicMatrix<F>& result_out_in,
                                     DynamicMatrix<F>& result_out_out,
                                     const XT::Common::Parameter& param) const
  {
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out);
    ensure_size_and_clear(result_out_in, rows_out, cols_in);
    ensure_size_and_clear(result_out_out, rows_out, cols_out);
  } // ... ensure_size_and_clear_results(...)
}; // class LocalQuaternaryIntersectionIntegrandInterface


} // namespace GDT
} // namespace Dune

#include "combined.hh"

#endif // DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH
