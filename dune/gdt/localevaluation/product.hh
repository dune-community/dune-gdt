// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_PRODUCT_HH
#define DUNE_GDT_EVALUATION_PRODUCT_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forward, to be used in the traits
template <class LocalizableFunctionImp>
class Product;


/**
 *  \brief Traits for the Product evaluation.
 */
template <class LocalizableFunctionImp>
class ProductTraits
{
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef Product<LocalizableFunctionType> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};


/**
 *  \brief  Computes a product evaluation.
 */
template <class LocalizableFunctionImp>
class Product : public LocalEvaluation::Codim0Interface<ProductTraits<LocalizableFunctionImp>, 1>,
                public LocalEvaluation::Codim0Interface<ProductTraits<LocalizableFunctionImp>, 2>,
                public LocalEvaluation::Codim1Interface<ProductTraits<LocalizableFunctionImp>, 1>,
                public LocalEvaluation::Codim1Interface<ProductTraits<LocalizableFunctionImp>, 2>
{
public:
  typedef ProductTraits<LocalizableFunctionImp> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  Product(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * These are the methods required to fulfill the interfaces
   * \{
   */

  /**
   * \brief extracts the local function and calls the correct order() method
   * \note  required by `LocalEvaluation::Codim0Interface< ..., 1 >`
   */
  template <class R, int rT, int rCT>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return order(*localFunction, testBase);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   * \note  required by `LocalEvaluation::Codim0Interface< ..., 1 >`
   */
  template <class R, int rT, int rCT>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, localPoint, ret);
  }

  /**
   * \brief extracts the local function and calls the correct order() method
   * \note  required by
   *        - `LocalEvaluation::Codim0Interface< ..., 2 >`
   *        - `LocalEvaluation::Codim1Interface< ..., 2 >`
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return order(*localFunction, testBase, ansatzBase);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   * \note  required by `LocalEvaluation::Codim0Interface< ..., 2 >`
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   * \note  required by `LocalEvaluation::Codim1Interface< ..., 1 >`
   */
  template <class IntersectionType, class R, int r, int rC>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, intersection, localPoint, ret);
  }

  /**
   * \brief extracts the local function and calls the correct evaluate() method
   * \note  required by `LocalEvaluation::Codim1Interface< ..., 2 >`
   */
  template <class IntersectionType, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFunctions);
    evaluate(*localFunction, testBase, ansatzBase, intersection, localPoint, ret);
  }

  /**
   * \}
   * end of: These are the methods required to fulfill the interfaces
   */

  /**
   * These are the methods that do the actual computations
   * \{
   */

  /**
   * \note implementation for `LocalEvaluation::Codim0Interface< ..., 1 >`
   */
  template <class R, int rL, int rCL, int rT, int rCT>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    return localFunction.order() + testBase.order();
  }

  /**
   * \note implementation for `LocalEvaluation::Codim0Interface< ..., 1 >`
   */
  template <class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    // checks
    typedef Dune::FieldVector<R, 1> RangeType;
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test base
    const size_t size = testBase.size();
    std::vector<RangeType> testValues(size, RangeType(0));
    testBase.evaluate(localPoint, testValues);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

  /**
   * \note implementation for
   *       - `LocalEvaluation::Codim0Interface< ..., 2 >`
   *       - `LocalEvaluation::Codim1Interface< ..., 2 >`
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  }

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  implementation for `LocalEvaluation::Codim0Interface< ..., 2 >`
   */
  template <class R, int r>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate bases
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>::RangeType RangeType;
    std::vector<RangeType> testValues(rows, RangeType(0));
    std::vector<RangeType> ansatzValues(cols, RangeType(0));
    testBase.evaluate(localPoint, testValues);
    ansatzBase.evaluate(localPoint, ansatzValues);
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
   * \note implementation for `LocalEvaluation::Codim1Interface< ..., 1 >`
   */
  template <class IntersectionType, class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    // checks
    typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
    typedef Dune::FieldVector<R, 1> RangeType;
    // evaluate local function
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const RangeType functionValue     = localFunction.evaluate(localPointEntity);
    // evaluate test base
    const size_t size = testBase.size();
    std::vector<RangeType> testValues(size, RangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

  /**
   * \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   * \note  implementation for `LocalEvaluation::Codim1Interface< ..., 2 >`
   */
  template <class IntersectionType, class R, int r>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // evaluate bases
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>::RangeType RangeType;
    std::vector<RangeType> testValues(rows, RangeType(0));
    std::vector<RangeType> ansatzValues(cols, RangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    ansatzBase.evaluate(localPointEntity, ansatzValues);
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
   * \}
   * end of: These are the methods that do the actual computations
   */

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class Product

} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_PRODUCT_HH
