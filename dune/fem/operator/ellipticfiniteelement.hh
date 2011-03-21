#ifndef DUNE_FEM_FUNCTIONALS_ELLIPTICFINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_ELLIPTICFINITEELEMENT_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/dofvector/localmatrix.hh>

namespace Dune
{

namespace Fem
{

namespace Functional
{

namespace Operator
{

template< class DiscreteFunctionSpaceImp, class InducingFunctionImp >
class EllipticFiniteElement
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef InducingFunctionImp
    InducingFunctionType;

  typedef typename InducingFunctionType::RangeFieldType
    RangeFieldType;

  typedef Dune::Fem::Functionals::LocalMatrix< RangeFieldType >
    LocalMatrixType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
    EntityIteratorType;

  typedef typename EntityIteratorType::Entity
    EntityType;

  EllipticFiniteElement(  const DiscreteFunctionSpaceType& discreteFunctionSpace,
                          const InducingFunctionType& inducingFunction )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      inducingFunction_( inducingFunction )
  {
  }

  ~EllipticFiniteElement()
  {
  }

  LocalMatrixType applyLocal( const EntityType& entity )
  {

    // some types we will need
    typedef typename EntityType::Geometry
      EntityGeometryType;

    typedef typename BaseFunctionSetType::RangeType
      RangeType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType
      JacobianRangeType;

    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet( entity );
    const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();

    // init return matrix
    LocalMatrixType ret( numberOfLocalDoFs, numberOfLocalDoFs );

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    typedef CachingQuadrature< GridPartType, 0 >
      EntityQuadratureType;
    typedef typename EntityQuadratureType::CoordinateType
      EntityCoordinateType;
    const unsigned quadratureOrder = discreteFunctionSpace_.order();
    const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
    const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

    // do loop over all local DoFs (i)
    for(  unsigned int i = 0;
          i < numberOfLocalDoFs;
          ++i )
    {

      // do loop over all local DoFs (i)
      for(  unsigned int j = 0;
            j < numberOfLocalDoFs;
            ++j )
      {

        // value of the operator, applied to the local basefunctions
        RangeFieldType operator_i_j = 0.0;

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

          // jacobian inverse transposed
          const typename EntityGeometryType::Jacobian& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed( xReferenceElement );

          // evaluate function and basefunction
          RangeType functionValue = 0.0;
          inducingFunction_.evaluate( xWorld, functionValue );

          JacobianRangeType gradientBaseFunction_i_untransposed = 0.0;
          baseFunctionSet.jacobian( i, xReferenceElement, gradientBaseFunction_i_untransposed );
          JacobianRangeType gradientBaseFunction_i = 0.0;
          jacobianInverseTransposed.mv( gradientBaseFunction_i_untransposed, gradientBaseFunction_i );

          JacobianRangeType gradientBaseFunction_j_untransposed = 0.0;
          baseFunctionSet.jacobian( i, xReferenceElement, gradientBaseFunction_j_untransposed );
          JacobianRangeType gradientBaseFunction_j = 0.0;
          jacobianInverseTransposed.mv( gradientBaseFunction_j_untransposed, gradientBaseFunction_j );

          // compute integral
          operator_i_j += integrationFactor * quadratureWeight * functionValue * gradientBaseFunction_i * gradientBaseFunction_j;

        } // done walk over quadrature points

        // set local matrix
        ret.set( i, j, operator_i_j );

      } // done loop over all local DoFs (i)

    } // done loop over all local DoFs (i)

  }

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const InducingFunctionType& inducingFunction_;

};

} // end of namespace Operator

} // end of namespace Functional

} // end of namespace Fem

} // end of namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_ELLIPTICFINITEELEMENT_HH
