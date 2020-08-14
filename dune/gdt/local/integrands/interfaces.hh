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
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/bound-object.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>
#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


// forwards (required for operator+), includes are below
template <class E, size_t r, size_t rC, class R, class F>
class LocalUnaryElementIntegrandSum;

template <class E, size_t t_r, size_t t_RC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryElementIntegrandSum;

template <class I, size_t r, size_t rC, class R, class F>
class LocalUnaryIntersectionIntegrandSum;

template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryIntersectionIntegrandSum;

template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalUnaryAndBinaryIntersectionIntegrandSum;

template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalQuaternaryIntersectionIntegrandSum;


// forwards (required for with_ansatz), includes are below
template <class E, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryToUnaryElementIntegrand;

template <class I, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
class LocalBinaryToUnaryIntersectionIntegrand;


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
  , public XT::Common::WithLogger<
        LocalUnaryElementIntegrandInterface<Element, range_dim, range_dim_cols, RangeField, Field>>
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalUnaryElementIntegrandInterface;
  using Logger = XT::Common::WithLogger<
      LocalUnaryElementIntegrandInterface<Element, range_dim, range_dim_cols, RangeField, Field>>;

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
  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

  LocalUnaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                      const std::string& logging_prefix = "",
                                      const std::string& logging_id_ = "",
                                      const bool logging_disabled = true)
    : XT::Common::ParametricInterface(param_type)
    , Logger(logging_prefix.empty() ? "gdt" : logging_prefix,
             logging_id_.empty() ? "LocalUnaryElementIntegrand" : logging_id_,
             logging_disabled)
  {}

  virtual ~LocalUnaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy_as_unary_element_integrand() const = 0;

  LocalUnaryElementIntegrandSum<E, r, rC, R, F> operator+(const ThisType& other) const
  {
    return LocalUnaryElementIntegrandSum<E, r, rC, R, F>(*this, other);
  }

  /**
   * Returns the polynomial order of the integrand, given the basis.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each function in the basis.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalTestBasisType& basis,
                        const DomainType& point_in_reference_element,
                        DynamicVector<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicVector<F> evaluate(const LocalTestBasisType& basis,
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
  , public XT::Common::WithLogger<LocalBinaryElementIntegrandInterface<Element,
                                                                       test_range_dim,
                                                                       test_range_dim_cols,
                                                                       TestRangeField,
                                                                       Field,
                                                                       ansatz_range_dim,
                                                                       ansatz_range_dim_cols,
                                                                       AnsatzRangeField>>
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalBinaryElementIntegrandInterface;
  using Logger = XT::Common::WithLogger<LocalBinaryElementIntegrandInterface<Element,
                                                                             test_range_dim,
                                                                             test_range_dim_cols,
                                                                             TestRangeField,
                                                                             Field,
                                                                             ansatz_range_dim,
                                                                             ansatz_range_dim_cols,
                                                                             AnsatzRangeField>>;

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

  LocalBinaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                       const std::string& logging_prefix = "",
                                       const bool logging_disabled = true)
    : XT::Common::ParametricInterface(param_type)
    , XT::Grid::ElementBoundObject<Element>()
    , Logger(logging_prefix.empty() ? "LocalBinaryElementIntegrand" : logging_prefix, logging_disabled)
  {
    LOG_(debug) << "LocalBinaryElementIntegrandInterface(this=" << this << ", param_type=" << param_type << ")"
                << std::endl;
  }

  virtual ~LocalBinaryElementIntegrandInterface() = default;

  template <class... Args>
  LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TR, F, a_r, a_rC, AR> with_ansatz(Args&&... args) const
  {
    return LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TR, F, a_r, a_rC, AR>(*this, std::forward<Args>(args)...);
  }

  virtual std::unique_ptr<ThisType> copy_as_binary_element_integrand() const = 0;

  LocalBinaryElementIntegrandSum<E, t_r, t_rC, TR, F, a_r, a_rC, AR> operator+(const ThisType& other) const
  {
    return LocalBinaryElementIntegrandSum<E, t_r, t_rC, TR, F, a_r, a_rC, AR>(*this, other);
  }

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

protected:
  void ensure_size_and_clear_results(const LocalTestBasisType& test_basis,
                                     const LocalAnsatzBasisType& ansatz_basis,
                                     DynamicMatrix<F>& result,
                                     const XT::Common::Parameter& param) const
  {
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
  } // ... ensure_size_and_clear_results(...)
}; // class LocalBinaryElementIntegrandInterface


template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F_ = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalUnaryAndBinaryElementIntegrandInterface
  : public LocalUnaryElementIntegrandInterface<I, t_r, t_rC, TR, F_>
  , public LocalBinaryElementIntegrandInterface<I, t_r, t_rC, TR, F_, a_r, a_rC, AR>
{

  using ThisType = LocalUnaryAndBinaryElementIntegrandInterface;

protected:
  using UnaryBaseType = LocalUnaryElementIntegrandInterface<I, t_r, t_rC, TR, F_>;
  using BinaryBaseType = LocalBinaryElementIntegrandInterface<I, t_r, t_rC, TR, F_, a_r, a_rC, AR>;

public:
  /// \name Members and typedefs required for disambiuation.
  /// \{

  using typename UnaryBaseType::D;
  using typename UnaryBaseType::DomainType;
  using typename UnaryBaseType::E;
  using typename UnaryBaseType::ElementType;
  using typename UnaryBaseType::F;
  using UnaryBaseType::d;

  /// \}

  LocalUnaryAndBinaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                               const std::string& logging_prefix = "",
                                               const std::string& logging_id_ = "",
                                               const bool logging_disabled = true)
    : UnaryBaseType(param_type, logging_prefix, logging_id_, logging_disabled)
    , BinaryBaseType(param_type, logging_prefix, logging_id_, logging_disabled)
  {}

  virtual ~LocalUnaryAndBinaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy_as_unary_and_binary_element_integrand() const = 0;

  /// \name Methods required for disambiguation.
  /// \{

  using UnaryBaseType::is_parametric;
  using UnaryBaseType::parameter_type;
  using UnaryBaseType::parse_parameter;

  /// \}
}; // class LocalUnaryAndBinaryElementIntegrandInterface


/**
 * Interface for integrands in integrals over grid intersections, which depend on one argument (usually the test basis
 * on the inside of the intersection in an integral-based functional).
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
class LocalUnaryIntersectionIntegrandInterface
  : virtual public XT::Common::ParametricInterface
  , virtual public XT::Grid::IntersectionBoundObject<Intersection>
  , public XT::Common::WithLogger<
        LocalUnaryIntersectionIntegrandInterface<Intersection, range_dim, range_dim_cols, RangeField, Field>>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalUnaryIntersectionIntegrandInterface;
  using BaseParametricType = XT::Common::ParametricInterface;
  using BaseIntersectionBoundObjectType = XT::Grid::IntersectionBoundObject<Intersection>;
  using Logger = XT::Common::WithLogger<
      LocalUnaryIntersectionIntegrandInterface<Intersection, range_dim, range_dim_cols, RangeField, Field>>;

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
  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, RF>;

  LocalUnaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                           const std::string& logging_prefix = "",
                                           const std::string& logging_id_ = "",
                                           const bool logging_disabled = true)
    : BaseParametricType(param_type)
    , BaseIntersectionBoundObjectType()
    , Logger(logging_prefix.empty() ? "gdt" : logging_prefix,
             logging_id_.empty() ? "LocalUnaryIntersectionIntegrand" : logging_id_,
             logging_disabled)
  {}

  LocalUnaryIntersectionIntegrandInterface(const ThisType&) = default;

  virtual ~LocalUnaryIntersectionIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy_as_unary_intersection_integrand() const = 0;

  LocalUnaryIntersectionIntegrandSum<I, r, rC, RF, F> operator+(const ThisType& other) const
  {
    return LocalUnaryIntersectionIntegrandSum<I, r, rC, RF, F>(*this, other);
  }

  /**
   * Flag to document which element the basis is expected to be bound to.
   */
  virtual bool inside() const
  {
    return true;
  }

  /**
   * Returns the polynomial order of the integrand, given the bases.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each basis function.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const LocalTestBasisType& basis,
                        const DomainType& point_in_reference_intersection,
                        DynamicVector<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicVector<F> evaluate(const LocalTestBasisType& basis,
                                    const DomainType& point_in_reference_intersection,
                                    const XT::Common::Parameter& param = {}) const
  {
    DynamicVector<F> result(basis.size(param), 0);
    evaluate(basis, point_in_reference_intersection, result, param);
    return result;
  }

protected:
  void ensure_size_and_clear_results(const LocalTestBasisType& basis,
                                     DynamicVector<F>& result,
                                     const XT::Common::Parameter& param) const
  {
    const size_t size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    result *= 0;
  }
}; // class LocalUnaryIntersectionIntegrandInterface


/**
 * Interface for integrands in integrals over grid intersections, which depend on two arguments (usually the test and
 * ansatz bases on the inside of the intersection in an integral-based functional).
 *
 * \note Regarding SMP: the integrand is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Intersection,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class TestRangeField = double,
          class Field = double,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols,
          class AnsatzRangeField = TestRangeField>
class LocalBinaryIntersectionIntegrandInterface
  : virtual public XT::Common::ParametricInterface
  , virtual public XT::Grid::IntersectionBoundObject<Intersection>
  , public XT::Common::WithLogger<LocalBinaryIntersectionIntegrandInterface<Intersection,
                                                                            test_range_dim,
                                                                            test_range_dim_cols,
                                                                            TestRangeField,
                                                                            Field,
                                                                            ansatz_range_dim,
                                                                            ansatz_range_dim_cols,
                                                                            AnsatzRangeField>>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalBinaryIntersectionIntegrandInterface;
  using BaseParametricType = XT::Common::ParametricInterface;
  using BaseIntersectionBoundObjectType = XT::Grid::IntersectionBoundObject<Intersection>;
  using Logger = XT::Common::WithLogger<LocalBinaryIntersectionIntegrandInterface<Intersection,
                                                                                  test_range_dim,
                                                                                  test_range_dim_cols,
                                                                                  TestRangeField,
                                                                                  Field,
                                                                                  ansatz_range_dim,
                                                                                  ansatz_range_dim_cols,
                                                                                  AnsatzRangeField>>;

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

  LocalBinaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                            const std::string& logging_prefix = "",
                                            const std::string& logging_id_ = "",
                                            const bool logging_disabled = true)
    : BaseParametricType(param_type)
    , BaseIntersectionBoundObjectType()
    , Logger(logging_prefix.empty() ? "gdt" : logging_prefix,
             logging_id_.empty() ? "LocalBinaryIntersectionIntegrand" : logging_id_,
             logging_disabled)
  {}

  LocalBinaryIntersectionIntegrandInterface(const ThisType& other) = default;

  virtual ~LocalBinaryIntersectionIntegrandInterface() = default;

  template <class... Args>
  LocalBinaryToUnaryIntersectionIntegrand<I, t_r, t_rC, TR, F, a_r, a_rC, AR> with_ansatz(Args&&... args) const
  {
    return LocalBinaryToUnaryIntersectionIntegrand<I, t_r, t_rC, TR, F, a_r, a_rC, AR>(*this,
                                                                                       std::forward<Args>(args)...);
  }

  virtual std::unique_ptr<ThisType> copy_as_binary_intersection_integrand() const = 0;

  LocalBinaryIntersectionIntegrandSum<I, t_r, t_rC, TR, F, a_r, a_rC, AR> operator+(const ThisType& other) const
  {
    return LocalBinaryIntersectionIntegrandSum<I, t_r, t_rC, TR, F, a_r, a_rC, AR>(*this, other);
  }

  /**
   * Flag to document which element the bases are expected to be bound to.
   */
  virtual bool inside() const
  {
    return true;
  }

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
                        const DomainType& point_in_reference_intersection,
                        DynamicMatrix<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual DynamicMatrix<F> evaluate(const LocalTestBasisType& test_basis,
                                    const LocalAnsatzBasisType& ansatz_basis,
                                    const DomainType& point_in_reference_intersection,
                                    const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result(test_basis.size(param), ansatz_basis.size(param), 0);
    evaluate(test_basis, ansatz_basis, point_in_reference_intersection, result, param);
    return result;
  }

protected:
  void ensure_size_and_clear_results(const LocalTestBasisType& test_basis,
                                     const LocalAnsatzBasisType& ansatz_basis,
                                     DynamicMatrix<F>& result,
                                     const XT::Common::Parameter& param) const
  {
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
  } // ... ensure_size_and_clear_results(...)
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
  , public XT::Common::WithLogger<LocalQuaternaryIntersectionIntegrandInterface<Intersection,
                                                                                test_range_dim,
                                                                                test_range_dim_cols,
                                                                                TestRangeField,
                                                                                Field,
                                                                                ansatz_range_dim,
                                                                                ansatz_range_dim_cols,
                                                                                AnsatzRangeField>>
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalQuaternaryIntersectionIntegrandInterface;
  using Logger = XT::Common::WithLogger<LocalQuaternaryIntersectionIntegrandInterface<Intersection,
                                                                                      test_range_dim,
                                                                                      test_range_dim_cols,
                                                                                      TestRangeField,
                                                                                      Field,
                                                                                      ansatz_range_dim,
                                                                                      ansatz_range_dim_cols,
                                                                                      AnsatzRangeField>>;

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

  LocalQuaternaryIntersectionIntegrandInterface(const XT::Common::ParameterType& param_type = {},
                                                const std::string& logging_prefix = "",
                                                const std::string& logging_id_ = "",
                                                const bool logging_disabled = true)
    : XT::Common::ParametricInterface(param_type)
    , XT::Grid::IntersectionBoundObject<Intersection>()
    , Logger(logging_prefix.empty() ? "gdt" : logging_prefix,
             logging_id_.empty() ? "LocalQuaternaryIntersectionIntegrand" : logging_id_,
             logging_disabled)
  {}

  virtual ~LocalQuaternaryIntersectionIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy_as_quaternary_intersection_integrand() const = 0;

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
#include "conversion.hh"

#endif // DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH
