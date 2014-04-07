// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
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
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
};


/**
 *  \brief  Computes a product evaluation.
 */
template< class LocalizableFunctionImp >
class Product
  : public LocalEvaluation::Codim0Interface< ProductTraits< LocalizableFunctionImp >, 1 >
  , public LocalEvaluation::Codim0Interface< ProductTraits< LocalizableFunctionImp >, 2 >
  , public LocalEvaluation::Codim1Interface< ProductTraits< LocalizableFunctionImp >, 1 >
{
public:
  typedef ProductTraits< LocalizableFunctionImp >   Traits;
  typedef typename Traits::LocalizableFunctionType  LocalizableFunctionType;

  Product(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalfunctionType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class E, class D, int d, class R, int rT, int rCT >
  size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    return order(*localFunction, testBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& localFunction,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase) const
  {
    return localFunction.order() + testBase.order();
  } // int order(...)

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    return order(*localFunction, testBase, ansatzBase);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& localFunction,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT >
   void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                 const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                 const Dune::FieldVector< D, d >& localPoint,
                 Dune::DynamicVector< R >& ret) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicVector< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  /**
   *  \brief computes a scalar product evaluation.
   */
  template< class E, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
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

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatzBase*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  /**
   *  \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   */
  template< class E, class D, int d, class R, int r >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate bases
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >::RangeType RangeType;
    std::vector< RangeType > testValues(rows, RangeType(0));
    std::vector< RangeType > ansatzValues(cols, RangeType(0));
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

  template< class E, class IntersectionType, class D, int d, class R, int r, int rC >
  void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, intersection, localPoint, ret);
  }

  template< class E, class IntersectionType, class D, int d, class R, int rL, int rCL, int rT, int rCT >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                const IntersectionType& /*intersection*/,
                const Dune::FieldVector< D, d - 1>& /*localPoint*/,
                Dune::DynamicVector< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  template< class E, class IntersectionType, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
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
