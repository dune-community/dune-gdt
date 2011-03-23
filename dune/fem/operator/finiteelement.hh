#ifndef DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH

// dune fem includes
#include <dune/fem/function/common/function.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/common/localmatrix.hh>
#include <dune/fem/common/localvector.hh>
#include <dune/fem/common/localbasefunction.hh>
#include <dune/fem/functional/ltwo.hh>

namespace Dune
{

namespace Functionals
{

namespace Operator
{

template< class DiscreteFunctionSpaceImp, class LocalOperationImp >
class FiniteElement
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef LocalOperationImp
    LocalOperationType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef Dune::Functionals::Common::LocalMatrix< RangeFieldType >
    LocalMatrixType;

  typedef Dune::Functionals::Common::LocalVector< RangeFieldType >
    LocalVectorType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  typedef typename DiscreteFunctionSpaceType::IteratorType
    EntityIteratorType;

  typedef typename EntityIteratorType::Entity
    EntityType;

  typedef typename EntityType::Geometry
    EntityGeometryType;

  typedef CachingQuadrature< GridPartType, 0 >
    EntityQuadratureType;

  typedef Dune::Functionals::Common::LocalBaseFunctionProvider< DiscreteFunctionSpaceType >
    LocalBaseFunctionProviderType;

  typedef typename LocalBaseFunctionProviderType::LocalBaseFunctionType
    LocalBaseFunctionType;

  FiniteElement(  const DiscreteFunctionSpaceType& discreteFunctionSpace,
                  const LocalOperationType& localOperation )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      localOperation_( localOperation ),
      localBaseFunctionProvider_( discreteFunctionSpace )
  {
  }

  ~FiniteElement()
  {
  }
  
  LocalMatrixType applyLocal( const EntityType& entity )
  {

    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet( entity );
    const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();

    // init return matrix
    LocalMatrixType ret( numberOfLocalDoFs, numberOfLocalDoFs );

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
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
          const DomainType x = entityQuadrature.point( quadraturePoint );

          // integration factors
          const double integrationFactor = entityGeometry.integrationElement( x );
          const double quadratureWeight = entityQuadrature.weight( quadraturePoint );

          // get local basefucntions
          const LocalBaseFunctionType phi_i = localBaseFunctionProvider_.provide( entity, i );
          const LocalBaseFunctionType phi_j = localBaseFunctionProvider_.provide( entity, j );

          // apply local operation
          const double localOperationEvaluated = localOperation_.operate( phi_i, phi_j, x );

          // compute integral
          operator_i_j += integrationFactor * quadratureWeight * localOperationEvaluated;

        } // done walk over quadrature points

        // set local matrix
        ret.set( i, j, operator_i_j );

      } // done loop over all local DoFs (i)

    } // done loop over all local DoFs (i)

    return ret;

  }

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const LocalOperationType& localOperation_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class FiniteElement

} // end of namespace Operator

} // end of namespace Functionals

} // end of namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
