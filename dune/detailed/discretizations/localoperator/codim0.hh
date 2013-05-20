#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH

#include <vector>
#include <utility>

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
namespace LocalOperator {


// forward, to be used in the traits
template <class BinaryEvaluationImp, class LocalizableFunctionImp>
class Codim0Integral;


template <class BinaryEvaluationImp, class LocalizableFunctionImp>
class Codim0IntegralTraits
{
public:
  typedef Codim0Integral<BinaryEvaluationImp, LocalizableFunctionImp> derived_type;
  typedef LocalEvaluation::BinaryInterface<typename BinaryEvaluationImp::Traits> BinaryEvaluationType;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf<Dune::Stuff::LocalizableFunction, LocalizableFunctionImp>::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
}; // class LocalOperatorCodim0IntegralTraits


template <class BinaryEvaluationImp, class LocalizableFunctionImp>
class Codim0Integral
    : public LocalOperator::Codim0Interface<Codim0IntegralTraits<BinaryEvaluationImp, LocalizableFunctionImp>>
{
public:
  typedef Codim0IntegralTraits<BinaryEvaluationImp, LocalizableFunctionImp> Traits;
  typedef typename Traits::BinaryEvaluationType BinaryEvaluationType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  Codim0Integral(const LocalizableFunctionType& function)
    : function_(function)
    , evaluation_()
  {
  }

  Codim0Integral(const LocalizableFunctionType& function, const BinaryEvaluationImp evaluation)
    : function_(function)
    , evaluation_(evaluation)
  {
  }

  template <class... Args>
  Codim0Integral(const LocalizableFunctionType& function, Args&&... args)
    : function_(function)
    , evaluation_(std::forward<Args>(args)...)
  {
  }

  const LocalizableFunctionType& inducingFunction() const
  {
    return function_;
  }

  const BinaryEvaluationType& inducingEvaluation() const
  {
    return evaluation_;
  }

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template <class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void apply(const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
             const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase, Dune::DynamicMatrix<R>& ret,
             std::vector<Dune::DynamicMatrix<R>>& tmpLocalMatrices) const
  {
    // local inducing function
    const auto& entity       = ansatzBase.entity();
    const auto localFunction = function_.localFunction(entity);
    // quadrature
    typedef Dune::QuadratureRules<D, d> VolumeQuadratureRules;
    typedef Dune::QuadratureRule<D, d> VolumeQuadratureType;
    const int quadratureOrder = evaluation_.order(localFunction, ansatzBase, testBase);
    assert(quadratureOrder >= 0 && "Not implemented for negative integration orders!");
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), 2 * quadratureOrder + 1);
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector<D, d> x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight  = quadPointIt->weight();
      // clear tmp matrix
      Dune::Stuff::Common::clear(tmpLocalMatrices[0]);
      // evaluate the local operation
      evaluation_.evaluate(localFunction, ansatzBase, testBase, x, tmpLocalMatrices[0]);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow       = ret[ii];
        const auto& tmpRow = tmpLocalMatrices[0][ii];
        for (size_t jj = 0; jj < cols; ++jj)
          retRow[jj] += tmpRow[jj] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const LocalizableFunctionType& function_;
  const BinaryEvaluationImp evaluation_;
}; // class Codim0Integral


} // namespace LocalOperator
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTEGRAL_HH
