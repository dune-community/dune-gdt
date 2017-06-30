// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class LocalizableFunctionImp, class Traits>
class LocalProductIntegrand;

template <class LocalizableFunctionImp, class Traits>
class LocalFVProductIntegrand;


namespace internal {


/**
 *  \brief Traits for the Product evaluation.
 */
template <class LocalizableFunctionImp>
class LocalProductIntegrandTraits
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be a localizable function.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef LocalProductIntegrand<LocalizableFunctionType, LocalProductIntegrandTraits> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
}; // class LocalProductIntegrandTraits

template <class LocalizableFunctionImp>
class LocalFVProductIntegrandTraits : public LocalProductIntegrandTraits<LocalizableFunctionImp>
{
public:
  typedef LocalFVProductIntegrand<LocalizableFunctionImp, LocalFVProductIntegrandTraits> derived_type;
}; // class LocalFVProductIntegrandTraits


} // namespace internal


/**
 *  \brief  Computes a product evaluation.
 */
template <class LocalizableFunctionImp, class TraitsImp = internal::LocalProductIntegrandTraits<LocalizableFunctionImp>>
class LocalProductIntegrand : public LocalVolumeIntegrandInterface<TraitsImp, 1>,
                              public LocalVolumeIntegrandInterface<TraitsImp, 2>,
                              public LocalFaceIntegrandInterface<TraitsImp, 1>,
                              public LocalFaceIntegrandInterface<TraitsImp, 2>
{
  using LocalVolumeIntegrandInterface<TraitsImp, 1>::as_imp;

public:
  typedef TraitsImp Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  LocalProductIntegrand(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {
  }

  /// \name Required by all variants of LocalVolumeIntegrandInterface
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /// \}
  /// \name Required by LocalVolumeIntegrandInterface< ..., 1 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT>
  size_t order(const LocalfunctionTupleType& localFuncs,
               const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>&
                   testBase) const
  {
    return this->as_imp().order(*std::get<0>(localFuncs), testBase);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class R, size_t rT, size_t rCT>
  void
  evaluate(const LocalfunctionTupleType& localFuncs,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
           const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    this->as_imp().evaluate(*std::get<0>(localFuncs), testBase, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalFaceIntegrandInterface< ..., 1 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void
  evaluate(const LocalfunctionTupleType& localFuncs,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    this->as_imp().evaluate(*std::get<0>(localFuncs), testBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(
      const LocalfunctionTupleType& localFuncs,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
      const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
      Dune::DynamicMatrix<R>& ret) const
  {
    this->as_imp().evaluate(*std::get<0>(localFuncs), testBase, ansatzBase, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalFaceIntegrandInterface< ..., 2 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(
      const LocalfunctionTupleType& localFunctions_in,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
      Dune::DynamicMatrix<R>& ret) const
  {
    this->as_imp().evaluate(*std::get<0>(localFunctions_in), testBase, ansatzBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 > and LocalFaceIntegrandInterface< ..., 2 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase)
      const
  {
    return this->as_imp().order(*std::get<0>(localFuncs), testBase, ansatzBase);
  }

  /// \}
  /// \name Actual implementations of order
  /// \{

  /**
   * \note   for `LocalVolumeIntegrandInterface< ..., 1 >`
   * \return localFunction.order() + testBase.order()
   */
  template <class R, size_t rL, size_t rCL, size_t rT, size_t rCT>
  size_t
  order(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase)
      const
  {
    return localFunction.order() + testBase.order();
  }

  /**
   * \note   for
   *         - `LocalVolumeIntegrandInterface< ..., 2 >`
   *         - `LocalFaceIntegrandInterface< ..., 2 >`
   * \return localFunction.order() + testBase.order() + ansatzBase.order()
   */
  template <class R, size_t rL, size_t rCL, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase)
      const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  }

  /// \}
  /// \name Actual implementations of evaluate
  /// \{

  /**
   * \note for `LocalVolumeIntegrandInterface< ..., 1 >`
   */
  template <class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test base
    const size_t size = testBase.size();
    const auto testValues = testBase.evaluate(localPoint);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  for `LocalVolumeIntegrandInterface< ..., 2 >`
   */
  template <class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
           const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate bases
    const auto rows = testBase.size();
    const auto cols = ansatzBase.size();
    const auto testValues = testBase.evaluate(localPoint);
    const auto ansatzValues = ansatzBase.evaluate(localPoint);
    // compute product
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        retRow[jj] = functionValue * (testValues[ii] * ansatzValues[jj]);
      }
    }
  } // ... evaluate(...)

  /**
   * \note for `LocalFaceIntegrandInterface< ..., 1 >`
   */
  template <class IntersectionType, class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate test base
    const size_t size = testBase.size();
    const auto testValues = testBase.evaluate(localPointEntity);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  for `LocalFaceIntegrandInterface< ..., 2 >`
   */
  template <class IntersectionType, class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0.0;
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate bases
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    const auto testValues = testBase.evaluate(localPointEntity);
    const auto ansatzValues = ansatzBase.evaluate(localPointEntity);
    // compute product
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        retRow[jj] = functionValue * (testValues[ii] * ansatzValues[jj]);
      }
    }
  } // ... evaluate(...)

  /// \}

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class LocalProductIntegrand


/**
 *  \brief  Computes a product evaluation for a finite volume basis
 */
template <class LocalizableFunctionImp,
          class TraitsImp = internal::LocalFVProductIntegrandTraits<LocalizableFunctionImp>>
class LocalFVProductIntegrand : public LocalProductIntegrand<LocalizableFunctionImp, TraitsImp>
{
  typedef LocalProductIntegrand<LocalizableFunctionImp, TraitsImp> BaseType;

public:
  using typename BaseType::Traits;
  using typename BaseType::LocalizableFunctionType;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  LocalFVProductIntegrand(const LocalizableFunctionType& inducingFunction)
    : BaseType(inducingFunction)
  {
  }

  using BaseType::evaluate;
  using BaseType::order;

  /// \}
  /// \name Actual implementations of evaluate
  /// \{

  /**
   * \note for `LocalVolumeIntegrandInterface< ..., 1 >`
   */
  template <class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test base
    const size_t size = testBase.size();
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue[ii];
    }
  } // ... evaluate(...)

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  for `LocalVolumeIntegrandInterface< ..., 2 >`
   */
  template <class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
           const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate bases
    const auto rows = testBase.size();
    const auto cols = ansatzBase.size();
    // compute product
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii)
      ret[ii][ii] = functionValue;
  } // ... evaluate(...)

  /**
   * \note for `LocalFaceIntegrandInterface< ..., 1 >`
   */
  template <class IntersectionType, class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    // evaluate local function
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate test base
    const size_t size = testBase.size();
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii)
      ret[ii] = functionValue[ii];
  } // ... evaluate(...)

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  for `LocalFaceIntegrandInterface< ..., 2 >`
   */
  template <class IntersectionType, class R, size_t r>
  void
  evaluate(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0.0;
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate bases
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    // compute product
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii)
      ret[ii][ii] = functionValue;
  } // ... evaluate(...)

  /// \}
}; // class LocalFVProductIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
