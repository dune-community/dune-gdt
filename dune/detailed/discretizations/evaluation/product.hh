#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Evaluation {


// forward, to be used in the traits
class Product;


/**
 *  \brief Traits for the Product evaluation.
 */
class ProductTraits
{
public:
  typedef Product derived_type;
};


/**
 *  \brief  Computes a product evaluation.
 */
class Product : public UnaryEvaluationInterface<ProductTraits>
{
public:
  typedef ProductTraits Traits;

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + testBase.order()
   */
  template <class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  int order(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
            const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase) const
  {
    if (localFunction.order() < 0)
      return -1;
    else
      return localFunction.order() + testBase.order();
  } // int order(...)

  template <class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& /*localFunction*/,
                       const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& /*testBase*/,
                       const Dune::FieldVector<D, d>& /*localPoint*/, Dune::DynamicVector<R>& /*ret*/)
  {
    dune_static_assert((Dune::AlwaysFalse<R>::value), "ERROR: not implemented for this combination of dimensions!");
  }

  /**
   *  \brief  Computes a product evaluation.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template <class L, class T, class D, int d, class R>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, 1, 1>& localFunction,
                       const BaseFunctionSetInterface<T, D, d, R, 1, 1>& testBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret)
  {
    // checks
    typedef Dune::FieldVector<R, 1> RangeType;
    // evaluate local function
    RangeType functionValue = localFunction.evaluate(localPoint);
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
}; // class Product

} // namespace Evaluation
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH
