
#ifndef DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
#define DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune fem-tools includes
#include "../../../tools/function/functiontools.hh" // should be removed in the end!

namespace Dune {

namespace Functionals {

/**
  * \brief      This class represents an L2 functional.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
template <class InducingFunctionImp>
class L2Functional
{
public:
  typedef InducingFunctionImp InducingFunctionType;

  typedef typename InducingFunctionType::RangeFieldType RangeFieldType;

  L2Functional(const InducingFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {
  }

  template <class DiscreteFunctionType>
  RangeFieldType operator()(const DiscreteFunctionType& discreteFunction) const
  {
    RangeFieldType ret = 0.0;

    // some types we will need
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;
    typedef typename EntityIteratorType::Entity EntityType;
    typedef typename EntityType::Geometry EntityGeometryType;
    typedef typename DiscreteFunctionType::RangeType RangeType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    // some things we will need
    const DiscreteFunctionSpaceType& discreteFunctionSpace = discreteFunction.space();

    // do gridwalk
    const EntityIteratorType BehindLastEntityIterator = discreteFunctionSpace.end();
    for (EntityIteratorType entityIterator = discreteFunctionSpace.begin(); entityIterator != BehindLastEntityIterator;
         ++entityIterator) {
      // entity and geometry
      const EntityType& entity                 = *entityIterator;
      const EntityGeometryType& entityGeometry = entity.geometry();

      // quadrature
      typedef CachingQuadrature<GridPartType, 0> EntityQuadratureType;
      typedef typename EntityQuadratureType::CoordinateType EntityCoordinateType;
      const int quadratureOrder = (2 * discreteFunctionSpace.order() + 1);
      const EntityQuadratureType entityQuadrature(entity, quadratureOrder);
      const unsigned int numberOfQuadraturePoints = entityQuadrature.nop();

      // local function and basefunction set
      const LocalFunctionType& localFunction    = discreteFunction.localFunction(entity);
      const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace.baseFunctionSet(entity);
      const unsigned numberOfLocalDoFs          = baseFunctionSet.numBaseFunctions();

      // do loop over all local DoFs
      for (unsigned int localDoF = 0; localDoF < numberOfLocalDoFs; ++localDoF) {
        // value of the local DoF
        const RangeFieldType localDoFValue = localFunction[localDoF];

        // value of the L2 functional, applied to the local basefunction, associated with the local DoF
        RangeFieldType localFunctionalValue = 0.0;
        // do walk over quadrature points
        for (unsigned int quadraturePoint = 0; quadraturePoint < numberOfQuadraturePoints; ++quadraturePoint) {
          // coordinates
          const EntityCoordinateType xReferenceElement = entityQuadrature.point(quadraturePoint);
          const EntityCoordinateType xWorld            = entityGeometry.global(xReferenceElement);

          // integration factors
          const double integrationFactor = entityGeometry.integrationElement(xReferenceElement);
          const double quadratureWeight  = entityQuadrature.weight(quadraturePoint);

          // evaluate function and basefunction
          RangeType functionValue = 0.0;
          inducingFunction_.evaluate(xWorld, functionValue);
          RangeType baseFunctionValue = 0.0;
          baseFunctionSet.evaluate(localDoF, xReferenceElement, baseFunctionValue);

          // compute integral
          localFunctionalValue += integrationFactor * quadratureWeight * functionValue * baseFunctionValue;

        } // done walk over quadrature points

        // add local product
        ret += localDoFValue * localFunctionalValue;

      } // done loop over all local DoFs

    } // done gridwalk

    return ret;
  }

private:
  const InducingFunctionType& inducingFunction_;

}; // end class L2Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
