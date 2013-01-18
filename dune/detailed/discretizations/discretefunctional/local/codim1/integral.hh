#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONAL_LOCAL_CODIM1_INTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONAL_LOCAL_CODIM1_INTEGRAL_HH

#include <dune/stuff/common/header/disable_warnings.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/header/reenable_warnings.hh>

#include <dune/stuff/common/vector.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace DiscreteFunctional {
namespace Local {
namespace Codim1 {
namespace Integral {

template <class LocalEvaluationImp>
class Boundary
{
public:
  typedef LocalEvaluationImp LocalEvaluationType;

  typedef Boundary<LocalEvaluationType> ThisType;

  typedef typename LocalEvaluationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  Boundary(const LocalEvaluationType& _localEvaluation)
    : localEvaluation_(_localEvaluation)
  {
  }

  const LocalEvaluationType& localEvaluation() const
  {
    return localEvaluation_;
  }

  unsigned int numTmpObjectsRequired() const
  {
    return 1;
  }

  template <class LocalTestBaseFunctionSetType, class IntersectionType, class LocalVectorType>
  void applyLocal(const LocalTestBaseFunctionSetType localTestBaseFunctionSet, const IntersectionType& intersection,
                  LocalVectorType& localVector, std::vector<LocalVectorType>& tmpLocalVectors) const
  {
    // sanity checks
    const unsigned int size = localTestBaseFunctionSet.size();
    assert(localVector.size() >= size);
    if (tmpLocalVectors.size() < numTmpObjectsRequired())
      tmpLocalVectors.resize(
          numTmpObjectsRequired(),
          LocalVectorType(localTestBaseFunctionSet.baseFunctionSet().space().map().maxLocalSize(), RangeFieldType(0)));

    // quadrature
    const unsigned int quadratureOrder = localEvaluation_.order() + localTestBaseFunctionSet.order();
    typedef Dune::QuadratureRules<DomainFieldType, IntersectionType::mydimension> FaceQuadratureRules;
    typedef Dune::QuadratureRule<DomainFieldType, IntersectionType::mydimension> FaceQuadratureType;
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule(intersection.type(), 2 * quadratureOrder + 1);

    // loop over all quadrature points
    for (typename FaceQuadratureType::const_iterator quadPoint = faceQuadrature.begin();
         quadPoint != faceQuadrature.end();
         ++quadPoint) {
      // coordinates
      const typename IntersectionType::LocalCoordinate xIntersection = quadPoint->position();
      const DomainType xEntity                                       = intersection.geometryInInside().global(xIntersection);
      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement(xIntersection);
      const double quadratureWeight  = quadPoint->weight();
      // clear target matrix
      Dune::Stuff::Common::clear(tmpLocalVectors[0]);
      // evaluate the local operation
      localEvaluation_.evaluateLocal(localTestBaseFunctionSet, xEntity, tmpLocalVectors[0]);
      // compute integral
      for (unsigned int i = 0; i < size; ++i)
        localVector[i] += tmpLocalVectors[0][i] * integrationFactor * quadratureWeight;
    } // loop over all quadrature points
  } // void apply(...)

private:
  Boundary(const ThisType&);
  ThisType& operator=(const ThisType&);

  const LocalEvaluationType& localEvaluation_;
}; // class Boundary

} // namespace Integral
} // namespace Codim1
} // namespace Local
} // namespace DiscreteFunctional
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONAL_LOCAL_CODIM0_INTEGRAL_HH
