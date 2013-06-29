// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_PRODUCT_HH
#define DUNE_GDT_EVALUATION_PRODUCT_HH

#include <tuple>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/gdt/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forward, to be used in the traits
template< class LocalizableFunctionImp >
class Product;


/**
 *  \brief Traits for the Product evaluation.
 */
template< class LocalizableFunctionImp >
class ProductTraits
{
public:
  typedef Product< LocalizableFunctionImp > derived_type;
  typedef LocalizableFunctionImp            LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableFunctionImp >::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
};


/**
 *  \brief  Computes a product evaluation.
 */
template< class LocalizableFunctionImp >
class Product
  : public LocalEvaluation::Codim0Interface< ProductTraits< LocalizableFunctionImp >, 1 >
  , public LocalEvaluation::Codim1Interface< ProductTraits< LocalizableFunctionImp >, 1 >
{
public:
  typedef ProductTraits< LocalizableFunctionImp >   Traits;
  typedef typename Traits::LocalizableFunctionType  LocalizableFunctionType;

  Product(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {}

  template< class EntityType >
  std::tuple< typename LocalizableFunctionType::template LocalFunction< EntityType >::Type > localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.localFunction(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class L, class T, class D, int d, class R, int rT, int rCT >
  int order(const std::tuple< L >& localFunctions,
            const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase) const
  {
    const auto& localFunction = std::get< 0 >(localFunctions);
    return order(localFunction, testBase);
  } // int order(...)

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template< class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  int order(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rL, rCL >& localFunction,
            const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase) const
  {
    if (localFunction.order() < 0)
      return -1;
    else
      return localFunction.order() + testBase.order();
  } // int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class L, class T, class D, int d, class R, int rT, int rCT >
   void evaluate(const std::tuple< L >& localFunctions,
                 const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
                 const Dune::FieldVector< D, d >& localPoint,
                 Dune::DynamicVector< R >& ret) const
  {
    const auto& localFunction = std::get< 0 >(localFunctions);
    evaluate(localFunction, testBase, localPoint, ret);
  }

  template< class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rL, rCL >& /*localFunction*/,
                const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& /*testBase*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicVector< R >& /*ret*/) const
  {
    dune_static_assert((Dune::AlwaysFalse< R >::value),
                       "ERROR: not implemented for this combination of dimensions!");
  }

  /**
   *  \brief computes a scalar product evaluation.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template< class L, class T, class D, int d, class R >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, 1, 1 >& localFunction,
                const BaseFunctionSetInterface< T, D, d, R, 1, 1 >& testBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    // checks
    typedef Dune::FieldVector< R, 1 > RangeType;
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test base
    const size_t size = testBase.size();
    std::vector< RangeType > testValues(size, RangeType(0));
    testBase.evaluate(localPoint, testValues);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

  template< class L, class T, class IntersectionType, class D, int d, class R, int r, int rC >
  void evaluate(const std::tuple< L >& localFunctions,
                const BaseFunctionSetInterface< T, D, d, R, r, rC >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    const auto& localFunction = std::get< 0 >(localFunctions);
    evaluate(localFunction, testBase, intersection, localPoint, ret);
  }

  template< class L, class T, class IntersectionType, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rL, rCL >& /*localFunction*/,
                const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& /*testBase*/,
                const IntersectionType& /*intersection*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicVector< R >& /*ret*/) const
  {
    dune_static_assert((Dune::AlwaysFalse< R >::value),
                       "ERROR: not implemented for this combination of dimensions!");
  }

  template< class L, class T, class IntersectionType, class D, int d, class R >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, 1, 1 >& localFunction,
                const BaseFunctionSetInterface< T, D, d, R, 1, 1 >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    // checks
    typedef Dune::FieldVector< D, d > DomainType;
    typedef Dune::FieldVector< R, 1 > RangeType;
    // evaluate local function
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const RangeType functionValue = localFunction.evaluate(localPointEntity);
    // evaluate test base
    const size_t size = testBase.size();
    std::vector< RangeType > testValues(size, RangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class Product

} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_PRODUCT_HH
