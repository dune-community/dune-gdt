// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/interfaces/local-functions.hh>

namespace Dune {
namespace GDT {


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
class LocalUnaryElementIntegrandInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalUnaryElementIntegrandInterface<Element, range_dim, range_dim_cols, RangeField, Field>;

public:
  using E = Element;
  using D = typename Element::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using ElementType = E;
  using DomainType = FieldVector<D, d>;
  using LocalBasisType = XT::Functions::LocalFunctionSetInterface<E, r, rC, R>;

  LocalUnaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalUnaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * This method needs to be called on each grid element before calling order() or evalaute(). It is supposed to be a
   * (nearly-)no-op if already bound to the element.
   */
  virtual ThisType& bind(const ElementType& element) = 0;

  /**
   * Returns the polynomial order of the integrand, given the basis.
   *
   * \note Undefined behaviour if not bound!
   */
  virtual int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each function in the basis.
   *
   * \note Undefined behaviour if not bound!
   */
  virtual void evaluate(const LocalBasisType& basis,
                        const DomainType& point_in_reference_element,
                        DynamicVector<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
    * This method is provided for convenience and should not be used within library code.
   *
   * \note Undefined behaviour if not bound!
    */
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
class LocalBinaryElementIntegrandInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalBinaryElementIntegrandInterface<Element,
                                                        test_range_dim,
                                                        test_range_dim_cols,
                                                        TestRangeField,
                                                        Field,
                                                        ansatz_range_dim,
                                                        ansatz_range_dim_cols,
                                                        AnsatzRangeField>;

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

  using ElementType = E;
  using DomainType = FieldVector<D, d>;
  using LocalTestBasisType = XT::Functions::LocalFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::LocalFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalBinaryElementIntegrandInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalBinaryElementIntegrandInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * This method needs to be called on each grid element before calling order() or evalaute(). It is supposed to be a
   * (nearly-)no-op if already bound to the element.
   */
  virtual ThisType& bind(const ElementType& element) = 0;

  /**
   * Returns the polynomial order of the integrand, given the bases.
   *
   * \note Undefined behaviour if not bound!
   */
  virtual int order(const LocalTestBasisType& test_basis,
                    const LocalAnsatzBasisType& ansatz_basis,
                    const XT::Common::Parameter& param = {}) const = 0;

  /**
   * Computes the evaluation of this integrand at the given point for each combination of functions from the two bases.
   *
   * \note Undefined behaviour if not bound!
   */
  virtual void evaluate(const LocalTestBasisType& test_basis,
                        const LocalAnsatzBasisType& ansatz_basis,
                        const DomainType& point_in_reference_element,
                        DynamicMatrix<F>& result,
                        const XT::Common::Parameter& param = {}) const = 0;

  /**
    * This method is provided for convenience and should not be used within library code.
   *
   * \note Undefined behaviour if not bound!
    */
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


#if 0
/**
 *  \brief  Interface for local evaluations that depend on an intersection.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 */
template <class Traits, size_t numArguments>
class LocalFaceIntegrandInterface
{
  static_assert(AlwaysFalse<Traits>::value, "There is no interface for this numArguments!");
};


/**
 *  \brief  Interface for unary codim 1 evaluations.
 */
template <class Traits>
class LocalFaceIntegrandInterface<Traits, 1> : public XT::CRTPInterface<LocalFaceIntegrandInterface<Traits, 1>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

protected:
  typedef EntityType E;
  typedef DomainFieldType D;
  static const size_t d = dimDomain;

public:
  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template <class R, size_t r, size_t rC>
  size_t order(const LocalfunctionTupleType& localFunctionsTuple,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctionsTuple, testBase));
    return this->as_imp().order(localFunctionsTuple, testBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r                   dimRange of the testBase
   *  \tparam rC                  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void evaluate(const LocalfunctionTupleType& localFunctionsTuple,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctionsTuple, testBase, intersection, localPoint, ret));
  }
}; // class LocalFaceIntegrandInterface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 1 evaluations.
 */
template <class Traits>
class LocalFaceIntegrandInterface<Traits, 2> : public XT::CRTPInterface<LocalFaceIntegrandInterface<Traits, 2>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

protected:
  typedef EntityType E;
  typedef DomainFieldType D;
  static const size_t d = dimDomain;

public:
  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase,ansatzBase}
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(const LocalfunctionTupleType& localFunctionsTuple,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctionsTuple, testBase, ansatzBase));
    return this->as_imp().order(localFunctionsTuple, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFunctionsTuple,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
                const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(
        this->as_imp().evaluate(localFunctionsTuple, testBase, ansatzBase, intersection, localPoint, ret));
  }
}; // class LocalFaceIntegrandInterface< Traits, 2 >


/**
 *  \brief  Interface for quaternary codim 1 evaluations.
 */
template <class Traits>
class LocalFaceIntegrandInterface<Traits, 4> : public XT::CRTPInterface<LocalFaceIntegrandInterface<Traits, 4>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

protected:
  typedef EntityType E;
  typedef DomainFieldType D;
  static const size_t d = dimDomain;

public:
  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().localFunctions(entity));
    return this->as_imp().localFunctions(entity);
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(const LocalfunctionTupleType localFunctionsEntity,
               const LocalfunctionTupleType localFunctionsNeighbor,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBaseEntity,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBaseEntity,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBaseNeighbor,
               const XT::Functions::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBaseNeighbor) const
  {
    CHECK_CRTP(this->as_imp().order(localFunctionsEntity,
                                    localFunctionsNeighbor,
                                    testBaseEntity,
                                    ansatzBaseEntity,
                                    testBaseNeighbor,
                                    ansatzBaseNeighbor));
    return this->as_imp().order(localFunctionsEntity,
                                localFunctionsNeighbor,
                                testBaseEntity,
                                ansatzBaseEntity,
                                testBaseNeighbor,
                                ansatzBaseNeighbor);
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFunctionsEntity,
                const LocalfunctionTupleType& localFunctionsNeighbor,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>&
                    testBaseEntity,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>&
                    ansatzBaseEntity,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>&
                    testBaseNeighbor,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>&
                    ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet,
                Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet,
                Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().evaluate(localFunctionsEntity,
                                                localFunctionsNeighbor,
                                                testBaseEntity,
                                                ansatzBaseEntity,
                                                testBaseNeighbor,
                                                ansatzBaseNeighbor,
                                                intersection,
                                                localPoint,
                                                entityEntityRet,
                                                neighborNeighborRet,
                                                entityNeighborRet,
                                                neighborEntityRet));
  }
}; // class LocalFaceIntegrandInterface< Traits, 4 >
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_INTERFACES_HH
