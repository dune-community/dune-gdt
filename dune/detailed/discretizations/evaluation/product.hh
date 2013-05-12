#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH

#include <dune/common/dynvector.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class EvaluationProduct;


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class EvaluationProductTraits
{
public:
  typedef EvaluationProduct<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> derived_type;
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
};


/**
 *  \brief Product evaluation for matrix valued inducing functions (not implemented).
 *
 *        See specialization for rangeDimCols = 1 for scalar and vector valued inducing functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class EvaluationProduct
{
public:
  EvaluationProduct() = delete;
};


template <class DomainFieldImp, int domainDim, class RangeFieldImp>
class EvaluationProduct<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public UnaryEvaluationInterface<EvaluationProductTraits<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>>
{
public:
  typedef EvaluationProductTraits<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;
  typedef typename Traits::DomainType DomainType;

  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = Traits::dimRangeRows;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  /**
   *  \brief  Computes a product evaluation.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   */
  template <class L, class T>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, DomainFieldType, dimDomain, RangeFieldType,
                                                                 dimRangeRows, dimRangeCols>& localFunction,
                       const BaseFunctionSetInterface<T>& testBase, const DomainType& localPoint,
                       Dune::DynamicVector<RangeFieldType>& ret)
  {
    // checks
    typedef Dune::FieldVector<RangeFieldType, dimRangeRows> RangeType;
    typedef typename BaseFunctionSetInterface<T>::RangeType TestRangeType;
    // evaluate local function
    RangeType functionValue(0);
    localFunction.evaluate(localPoint, functionValue);
    // evaluate test base
    const size_t size = testBase.size();
    std::vector<TestRangeType> testValues(size, TestRangeType(0));
    testBase.evaluate(localPoint, testValues);
    // compute product
    assert(ret.size() >= size);
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] = functionValue * testValues[ii];
    }
  } // ... evaluate(...)
}; // class EvaluationProduct

} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_PRODUCT_HH
