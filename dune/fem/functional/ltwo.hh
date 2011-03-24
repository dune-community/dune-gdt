
#ifndef DUNE_FEM_FUNCTIONALS_FUNCTIONAL_LTWO_HH
#define DUNE_FEM_FUNCTIONALS_FUNCTIONAL_LTWO_HH

// dune common includes
#include <dune/common/fvector.hh>

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune fem-functionals includes
#include <dune/fem/common/localvector.hh>

namespace Dune
{

namespace Functionals
{

namespace Functional
{

/**
  * \brief      This class represents an L2 functional.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
template< class DiscreteFunctionSpaceImp, class InducingFunctionImp >
class L2
{
public:

  typedef InducingFunctionImp
    InducingFunctionType;

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef typename InducingFunctionType::RangeFieldType
    RangeFieldType;

  typedef Dune::Functionals::Common::LocalVector< RangeFieldType >
    LocalDoFVectorType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
    EntityIteratorType;

  typedef typename EntityIteratorType::Entity
    EntityType;

  L2( const DiscreteFunctionSpaceType& discreteFunctionSpace, const InducingFunctionType& inducingFunction )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      inducingFunction_( inducingFunction )
  {
  }

  ~L2()
  {
  }

  /**
    * \brief      This operator represents the application of the functional to a discrete function.
    *
    * \todo       Doc me, please!
    **/
  template< class DiscreteFunctionType >
  const RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
  {
    RangeFieldType ret = 0.0;

    // some types we will need
    typedef typename DiscreteFunctionType::RangeType
      RangeType;
    typedef typename DiscreteFunctionType::LocalFunctionType
      LocalFunctionType;

    // do gridwalk
    const EntityIteratorType BehindLastEntity = discreteFunctionSpace_.end();
    for ( EntityIteratorType entityIterator = discreteFunctionSpace_.begin();
          entityIterator != BehindLastEntity;
          ++entityIterator )
    {
      // entity and geometry
      const EntityType& entity = *entityIterator;

      // local function and basefunction set
      const LocalFunctionType& localFunction = discreteFunction.localFunction( entity );
      const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet( entity );

      // local DoF and functional vector
      const LocalDoFVectorType localDoFVector( localFunction );
      const LocalDoFVectorType localFunctionalVector = applyLocal( entity, baseFunctionSet );

      // compute product
      ret += localDoFVector * localFunctionalVector;

    } // done gridwalk

    return ret;

  } // end operator()

  /**
    * \brief      This function respresents the local application of the functional
    *             to a local basefunctionset.
    *
    * \todo       Doc me, please!
    **/
  const LocalDoFVectorType applyLocal(  const EntityType& entity,
                                        const BaseFunctionSetType& baseFunctionSet ) const
  {
    // some types we will need
    typedef typename EntityType::Geometry
      EntityGeometryType;
    typedef typename BaseFunctionSetType::RangeType
      RangeType;

    // init return vector
    const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();
    LocalDoFVectorType ret( numberOfLocalDoFs );

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    typedef CachingQuadrature< GridPartType, 0 >
      EntityQuadratureType;
    typedef typename EntityQuadratureType::CoordinateType
      EntityCoordinateType;
    const unsigned quadratureOrder = (2 * discreteFunctionSpace_.order() + 1);
    const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
    const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

    // do loop over all local DoFs
    for(  unsigned int localDoF = 0;
          localDoF < numberOfLocalDoFs;
          ++localDoF )
    {
      // value of the L2 functional, applied to the local basefunction, associated with the local DoF
      RangeFieldType localFunctionalValue = 0.0;

      // do walk over quadrature points
      for(  unsigned int quadraturePoint = 0;
            quadraturePoint < numberOfQuadraturePoints;
            ++quadraturePoint )
      {
        // coordinates
        const EntityCoordinateType xReferenceElement = entityQuadrature.point( quadraturePoint );
        const EntityCoordinateType xWorld = entityGeometry.global( xReferenceElement );

        // integration factors
        const double integrationFactor = entityGeometry.integrationElement( xReferenceElement );
        const double quadratureWeight = entityQuadrature.weight( quadraturePoint );

        // evaluate function and basefunction
        RangeType functionValue = 0.0;
        inducingFunction_.evaluate( xWorld, functionValue );
        RangeType baseFunctionValue = 0.0;
        baseFunctionSet.evaluate( localDoF, xReferenceElement, baseFunctionValue );

        // compute integral
        localFunctionalValue += integrationFactor * quadratureWeight * functionValue * baseFunctionValue;

      } // done walk over quadrature points

      // set local vector
      ret[localDoF] = localFunctionalValue;

    } // done loop over all local DoFs

    return ret;

  } // end applyLocal()

  const DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const InducingFunctionType& inducingFunction_;

}; // end class L2

} // end namespace Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_FUNCTIONAL_LTWO_HH
