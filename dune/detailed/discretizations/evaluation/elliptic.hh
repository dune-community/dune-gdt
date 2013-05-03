#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH

#include <memory>

#include <dune/common/dynmatrix.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class EvaluationElliptic;


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class EvaluationEllipticTraits
{
public:
  typedef EvaluationElliptic<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> derived_type;
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
};


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class EvaluationElliptic
{
public:
  EvaluationElliptic() = delete;
};


template <class DomainFieldImp, int domainDim, class RangeFieldImp>
class EvaluationElliptic<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public BinaryEvaluationInterface<EvaluationEllipticTraits<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>>
{
public:
  typedef EvaluationEllipticTraits<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;
  typedef typename Traits::DomainType DomainType;

  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = Traits::dimRangeRows;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  /**
   *  \brief  Computes an elliptic evaluation.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   */
  template <class L, class T, class A>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, DomainFieldType, dimDomain, RangeFieldType,
                                                                 dimRangeRows, dimRangeCols>& localFunction,
                       const BaseFunctionSetInterface<T>& testBase, const BaseFunctionSetInterface<A>& ansatzBase,
                       const DomainType& localPoint, Dune::DynamicMatrix<RangeFieldType>& ret)
  {
    // check
    typedef typename BaseFunctionSetInterface<A>::JacobianRangeType AnsatzJacobianRangeType;
    typedef typename BaseFunctionSetInterface<A>::JacobianRangeType TestJacobianRangeType;
    typedef Dune::FieldVector<RangeFieldType, dimRangeRows> RangeType;
    // evaluate local function
    RangeType functionValue(0);
    localFunction.evaluate(localPoint, functionValue);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector<TestJacobianRangeType> testGradients(rows, TestJacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector<AnsatzJacobianRangeType> ansatzGradients(cols, AnsatzJacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        const RangeFieldType gradientProduct = ansatzGradients[jj][0] * testGradients[jj][0];
        ret[ii][jj]                          = functionValue * gradientProduct;
      }
    }
  } // ... evaluate(...)
}; // end class EvaluationElliptic


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif //
