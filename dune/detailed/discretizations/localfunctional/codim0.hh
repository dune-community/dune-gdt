#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH

#include <vector>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/vector.hh>

#include "../basefunctionset/interface.hh"
#include "../evaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template <class UnaryEvaluationImp, class LocalizableFunctionImp>
class LocalFunctionalCodim0Integral;


template <class UnaryEvaluationImp, class LocalizableFunctionImp>
class LocalFunctionalCodim0IntegralTraits
{
public:
  typedef LocalFunctionalCodim0Integral<UnaryEvaluationImp, LocalizableFunctionImp> derived_type;
  typedef UnaryEvaluationInterface<typename UnaryEvaluationImp::Traits> UnaryEvaluationType;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf<Dune::Stuff::LocalizableFunction, LocalizableFunctionImp>::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
}; // class LocalFunctionalCodim0IntegralTraits


template <class UnaryEvaluationImp, class LocalizableFunctionImp>
class LocalFunctionalCodim0Integral
    : public LocalFunctionalCodim0Interface<LocalFunctionalCodim0IntegralTraits<UnaryEvaluationImp,
                                                                                LocalizableFunctionImp>>
{
public:
  typedef LocalFunctionalCodim0IntegralTraits<UnaryEvaluationImp, LocalizableFunctionImp> Traits;
  typedef typename Traits::UnaryEvaluationType UnaryEvaluationType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  LocalFunctionalCodim0Integral(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {
  }

  const LocalizableFunctionType& inducingFunction() const
  {
    return inducingFunction_;
  }

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  /**
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \attention  ret is assumed to be zero!
   */
  template <class T, class D, int d, class R, int r, int rC>
  void apply(const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase, Dune::DynamicVector<R>& ret,
             std::vector<Dune::DynamicVector<R>>& tmpLocalVectors) const
  {
    // local inducing function
    const auto& entity       = testBase.entity();
    const auto localFunction = inducingFunction_.localFunction(entity);
    // quadrature
    typedef Dune::QuadratureRules<D, d> VolumeQuadratureRules;
    typedef Dune::QuadratureRule<D, d> VolumeQuadratureType;
    assert(localFunction.order() >= 0 && "Not implemented for negative integration orders!");
    const size_t quadratureOrder                 = std::max(int(localFunction.order()), 0) + std::max(int(testBase.order()), 0);
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), 2 * quadratureOrder + 1);
    // check vector and tmp storage
    const size_t size = testBase.size();
    assert(ret.size() >= size);
    assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector<D, d> x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight  = quadPointIt->weight();
      // clear tmp vector
      Dune::Stuff::Common::clear(tmpLocalVectors[0]);
      // evaluate the local operation
      UnaryEvaluationType::evaluate(localFunction, testBase, x, tmpLocalVectors[0]);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii) {
        ret[ii] += tmpLocalVectors[0][ii] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class LocalFunctionalCodim0Integral


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH
