#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH

#include <vector>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/matrix.hh>

#include "../basefunctionset/interface.hh"
#include "../evaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template< class BinaryEvaluationImp, class LocalizableFunctionImp >
class LocalOperatorCodim0Integral;


template< class BinaryEvaluationImp, class LocalizableFunctionImp >
class LocalOperatorCodim0IntegralTraits
{
public:
  typedef LocalOperatorCodim0Integral< BinaryEvaluationImp, LocalizableFunctionImp > derived_type;
  typedef BinaryEvaluationInterface< typename BinaryEvaluationImp::Traits > BinaryEvaluationType;
  typedef LocalizableFunctionImp                                            LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableFunctionImp >::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
}; // class LocalOperatorCodim0IntegralTraits


template< class BinaryEvaluationImp, class LocalizableFunctionImp >
class LocalOperatorCodim0Integral
  : public LocalOperatorCodim0Interface< LocalOperatorCodim0IntegralTraits< BinaryEvaluationImp, LocalizableFunctionImp > >
{
public:
  typedef LocalOperatorCodim0IntegralTraits< BinaryEvaluationImp, LocalizableFunctionImp > Traits;
  typedef typename Traits::BinaryEvaluationType     BinaryEvaluationType;
  typedef typename Traits::LocalizableFunctionType  LocalizableFunctionType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  LocalOperatorCodim0Integral(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {}

  const LocalizableFunctionType& inducingFunction() const
  {
    return inducingFunction_;
  }

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  /**
   *  \brief      Applies the local operator.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \tparam A   Traits of the BaseFunctionSetInterface implementation, representing the type of the ansatzBase
   *  \attention  ret is assumed to be zero!
   */
  template< class T, class A, class RangeFieldType >
  void apply(const BaseFunctionSetInterface< T >& testBase,
             const BaseFunctionSetInterface< A >& ansatzBase,
             Dune::DynamicMatrix< RangeFieldType >& ret,
             std::vector< Dune::DynamicMatrix< RangeFieldType > >& tmpLocalMatrices) const
  {
    // checks
    typedef BaseFunctionSetInterface< A > AnsatzBaseType;
    typedef BaseFunctionSetInterface< T > TestBaseType;
    static const unsigned int dimDomain = AnsatzBaseType::dimDomain;
    dune_static_assert((dimDomain == TestBaseType::dimDomain), "ERROR: BaseFunctionSets do not match!");
    typedef typename AnsatzBaseType::DomainFieldType DomainFieldType;
    typedef typename AnsatzBaseType::DomainType DomainType;
    // local inducing function
    const auto& entity = ansatzBase.entity();
    const auto localFunction = inducingFunction_.localFunction(entity);
    // quadrature
    typedef Dune::QuadratureRules< DomainFieldType, dimDomain > VolumeQuadratureRules;
    typedef Dune::QuadratureRule< DomainFieldType, dimDomain > VolumeQuadratureType;
    assert(localFunction.order() >= 0 && "Not implemented for negative integration orders!");
    const size_t quadratureOrder = std::max(int(localFunction.order()), 0)
                                   + std::max(int(ansatzBase.order()) - 1, 0)
                                   + std::max(int(testBase.order()) - 1, 0);
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), 2*quadratureOrder + 1);
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const DomainType x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight = quadPointIt->weight();
      // clear tmp matrix
      Dune::Stuff::Common::clear(tmpLocalMatrices[0]);
      // evaluate the local operation
      BinaryEvaluationType::evaluate(localFunction, ansatzBase, testBase, x, tmpLocalMatrices[0]);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow = ret[ii];
        const auto& tmpRow = tmpLocalMatrices[0][ii];
        for (size_t jj = 0; jj < cols; ++jj)
          retRow[jj] += tmpRow[jj] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class LocalOperatorCodim0Integral


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH
