#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Evaluation {


// forward, to be used in the traits
// template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
class Elliptic;


// template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
class EllipticTraits
{
public:
  typedef Elliptic /*< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >*/ derived_type;
  //  typedef DomainFieldImp                                  DomainFieldType;
  //  static const unsigned int                               dimDomain = domainDim;
  //  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  //  typedef RangeFieldImp                                 RangeFieldType;
  //  static const unsigned int                             dimRangeRows = rangeDimRows;
  //  static const unsigned int                             dimRangeCols = rangeDimCols;
};


class Elliptic : public BinaryEvaluationInterface<EllipticTraits>
{
public:
  typedef /*Evaluation*/ EllipticTraits /*< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >*/ Traits;

  //  typedef typename Traits::DomainFieldType  DomainFieldType;
  //  static const unsigned int                 dimDomain = Traits::dimDomain;
  //  typedef typename Traits::DomainType       DomainType;

  //  typedef typename Traits::RangeFieldType RangeFieldType;
  //  static const unsigned int               dimRangeRows = Traits::dimRangeRows;
  //  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  template <class L, class T, class A, class D, int d, class R, int r, int rC>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, r, rC>& /*localFunction*/,
                       const BaseFunctionSetInterface<T, D, d, R, r, rC>& /*testBase*/,
                       const BaseFunctionSetInterface<A, D, d, R, r, rC>& /*ansatzBase*/,
                       const Dune::FieldVector<D, d>& /*localPoint*/, Dune::DynamicMatrix<R>& /*ret*/)
  {
    dune_static_assert((Dune::AlwaysFalse<R>::value),
                       "ERROR: not implemented for this combination of dimensions d, r and rC!");
  }


  /**
   *  \brief  Computes an elliptic evaluationf for a scalar local function and scalar basefunctionsets.
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   */
  template <class L, class T, class A, class D, int d, class RangeFieldType>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, RangeFieldType, 1, 1>& localFunction,
                       const BaseFunctionSetInterface<T, D, d, RangeFieldType, 1, 1>& testBase,
                       const BaseFunctionSetInterface<A, D, d, RangeFieldType, 1, 1>& ansatzBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<RangeFieldType>& ret)
  {
    typedef typename BaseFunctionSetInterface<A, D, d, RangeFieldType, 1, 1>::RangeType RangeType;
    typedef typename BaseFunctionSetInterface<A, D, d, RangeFieldType, 1, 1>::JacobianRangeType JacobianRangeType;
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const RangeFieldType gradientProduct = ansatzGradients[jj][0] * testGradients[ii][0];
        retRow[jj]                           = functionValue * gradientProduct;
      }
    }
  } // ... evaluate(...)
}; // class Elliptic


} // namespace Evaluation
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif //
