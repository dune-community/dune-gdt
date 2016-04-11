// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class LocalizableFunctionImp>
class LocalProductIntegrand;


namespace internal {


/**
 *  \brief Traits for the Product evaluation.
 */
template <class LocalizableFunctionImp>
class LocalProductIntegrandTraits
{
  static_assert(Stuff::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be a localizable function.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef LocalProductIntegrand<LocalizableFunctionType> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
}; // class LocalProductIntegrandTraits


} // namespace internal


/**
 *  \brief  Computes a product evaluation.
 */
template <class LocalizableFunctionImp>
class LocalProductIntegrand
    : public LocalVolumeIntegrandInterface<internal::LocalProductIntegrandTraits<LocalizableFunctionImp>, 1>,
      public LocalVolumeIntegrandInterface<internal::LocalProductIntegrandTraits<LocalizableFunctionImp>, 2>,
      public LocalFaceIntegrandInterface<internal::LocalProductIntegrandTraits<LocalizableFunctionImp>, 1>,
      public LocalFaceIntegrandInterface<internal::LocalProductIntegrandTraits<LocalizableFunctionImp>, 2>
{
public:
  typedef internal::LocalProductIntegrandTraits<LocalizableFunctionImp> Traits;
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
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    return order(*std::get<0>(localFuncs), testBase);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class R, size_t rT, size_t rCT>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    evaluate(*std::get<0>(localFuncs), testBase, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalFaceIntegrandInterface< ..., 1 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    evaluate(*std::get<0>(localFuncs), testBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    evaluate(*std::get<0>(localFuncs), testBase, ansatzBase, localPoint, ret);
  }

  /// \}
  /// \name Required by LocalFaceIntegrandInterface< ..., 2 >
  /// \{

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFunctions_in,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    evaluate(*std::get<0>(localFunctions_in), testBase, ansatzBase, intersection, localPoint, ret);
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
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return order(*std::get<0>(localFuncs), testBase, ansatzBase);
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
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
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
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
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
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test base
    const size_t size     = testBase.size();
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
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate bases
    const auto rows         = testBase.size();
    const auto cols         = ansatzBase.size();
    const auto testValues   = testBase.evaluate(localPoint);
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
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    ret *= 0.0;
    // evaluate local function
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto functionValue    = localFunction.evaluate(localPointEntity);
    // evaluate test base
    const size_t size     = testBase.size();
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
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    ret *= 0.0;
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate bases
    const size_t rows       = testBase.size();
    const size_t cols       = ansatzBase.size();
    const auto testValues   = testBase.evaluate(localPointEntity);
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

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
