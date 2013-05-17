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
//template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
class Product;


//template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
class /*Evaluation*/ProductTraits
{
public:
  typedef /*Evaluation*/Product/*< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >*/ derived_type;
//  typedef DomainFieldImp                                  DomainFieldType;
//  static const unsigned int                               dimDomain = domainDim;
//  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
//  typedef RangeFieldImp                                 RangeFieldType;
//  static const unsigned int                             dimRangeRows = rangeDimRows;
//  static const unsigned int                             dimRangeCols = rangeDimCols;
};


///**
// *  \brief Product evaluation for matrix valued inducing functions (not implemented).
// *
// *        See specialization for rangeDimCols = 1 for scalar and vector valued inducing functions.
// */
//template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols >
//class EvaluationProduct
//{
//public:
//  EvaluationProduct() = delete;
//};


//template< class DomainFieldImp, int domainDim, class RangeFieldImp >
class /*Evaluation*/Product/*< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >*/
  : public UnaryEvaluationInterface< /*Evaluation*/ProductTraits/*< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >*/ >
{
public:
  typedef /*Evaluation*/ProductTraits/*< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >*/ Traits;

//  typedef typename Traits::DomainFieldType  DomainFieldType;
//  static const unsigned int                 dimDomain = Traits::dimDomain;
//  typedef typename Traits::DomainType       DomainType;

//  typedef typename Traits::RangeFieldType RangeFieldType;
//  static const unsigned int               dimRangeRows = Traits::dimRangeRows;
//  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  template< class L, class T, class A, class D, int d, class R, int r, int rC >
  static void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, r, rC >& /*localFunction*/,
                       const BaseFunctionSetInterface< T, D, d, R, r, rC >& /*testBase*/,
                       const Dune::FieldVector< D, d >& /*localPoint*/,
                       Dune::DynamicMatrix< R >& /*ret*/)
  {
    dune_static_assert((Dune::AlwaysFalse< R >::value),
                       "ERROR: not implemented for this combination of dimensions d, r and rC!");
  }

  /**
   *  \brief  Computes a product evaluation.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   */
  template< class L, class T, class D, int d, class R >
  static void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, 1, 1 >& localFunction,
                       const BaseFunctionSetInterface< T, D, d, R, 1, 1 >& testBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicVector< R >& ret)
  {
    // checks
    typedef Dune::FieldVector< R, 1 > RangeType;
    // evaluate local function
    RangeType functionValue = localFunction.evaluate(localPoint);
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
}; // class Product

} // namespace Evaluation
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH
