// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_PRODUCT_HH
#define DUNE_GDT_EVALUATION_PRODUCT_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

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
                public LocalEvaluation::Codim1Interface<ProductTraits<LocalizableFunctionImp>, 1>
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

  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunction*/,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/,
           const Dune::FieldVector<DomainFieldType, dimDomain>& /*localPoint*/, Dune::DynamicMatrix<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  } // ... evaluate(...)

  template <class IntersectionType, class R, int rL, int rCL, int rT, int rCT>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunction*/,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
           const IntersectionType& /*intersection*/,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& /*localPoint*/,
           Dune::DynamicVector<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  } // ... evaluate(...)

  /**
   * \defgroup codim0_1 ´´Required by LocalEvaluation::Codim0Interface< ProductTraits< LocalizableFunctionImp >, 1 >!``
   * \{
   */

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, int rT, int rCT>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return redirect_order(*localFunction, testBase);
  } // ... order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class R, int rT, int rCT>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    redirect_evaluate(*localFunction, testBase, localPoint, ret);
  } // ... evaluate(...)

  /**
   * \}
   */

  /**
   * \defgroup codim0_2 ´´Required by LocalEvaluation::Codim0Interface< ProductTraits< LocalizableFunctionImp >, 2 >!``
   * \{
   */

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return redirect_order(*localFunction, testBase, ansatzBase);
  } // ... order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    redirect_evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  } // ... evaluate(...)

  /**
   * \}
   */

  /**
   * \defgroup codim1_1 ´´Required by LocalEvaluation::Codim1Interface< ProductTraits< LocalizableFunctionImp >, 1 >!``
   * \{
   */

  template <class IntersectionType, class R, int r, int rC>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    redirect_evaluate(*localFunction, testBase, intersection, localPoint, ret);
  } // ... evaluate(...)

  /**
   * \}
   */

private:
  template <class R, int rL, int rCL, int rT, int rCT>
  void redirect_evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunction*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
      const Dune::FieldVector<DomainFieldType, dimDomain>& /*localPoint*/, Dune::DynamicVector<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  } // ... redirect_evaluate(...)

  /**
   * \addtogroup codim1_1
   * \addtogroup codim2_1
   * \return localFunction.order() + testBase.order()
   */
  template <class R, int rL, int rCL, int rT, int rCT>
  size_t redirect_order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    return localFunction.order() + testBase.order();
  } // size_t redirect_order(...)

  /**
   *  \addtogroup codim1_1
   *  \brief computes a scalar product evaluation.
   */
  template <class R>
  void
  redirect_evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
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
  } // ... redirect_evaluate(...)

  /**
   *  \addtogroup codim1_2
   *  \{
   *  \return localFunction.order() + testBase.order() + ansatzBase.order()
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  size_t redirect_order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  } // ... redirect_order(...)

  /**
   *  \brief Computes a product evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   */
  template <class R, int r>
  void
  redirect_evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
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
  } // ... redirect_evaluate(...)

  /**
   * \}
   */

  /**
   * \addtogroup codim2_1
   */
  template <class IntersectionType, class R>
  void
  redirect_evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                    const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                    const IntersectionType& intersection,
                    const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
                    Dune::DynamicVector<R>& ret) const
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
  } // ... redirect_evaluate(...)

  const LocalizableFunctionType& inducingFunction_;
}; // class Product

} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_PRODUCT_HH
