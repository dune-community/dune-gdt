
#ifndef DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH

// dune common includes
#include <dune/common/fvector.hh>

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune fem-functionals includes
#include <dune/fem/common/localbasefunction.hh>
#include <dune/fem/common/localvector.hh>
#include <dune/fem/localoperation/interface.hh>

namespace Dune
{

namespace Functionals
{

namespace Functional
{

/**
  * \brief      This class represents a finite element functional.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
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

  /**
    * \brief  Copy constructor
    **/
  FiniteElement( const FiniteElement& other )
    : discreteFunctionSpace_( other.space() ),
      localOperation_( other.localOperation() ),
      localBaseFunctionProvider_( discreteFunctionSpace_ )
  {
  }

  /**
    * \brief  Destructor
    *         Does nothing.
    **/
  ~FiniteElement()
  {
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

  const LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  const LocalVectorType applyLocal( const EntityType& entity ) const
  {
    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet( entity );
    const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();

    // init return vector
    LocalVectorType ret( numberOfLocalDoFs );

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    const unsigned quadratureOrder = discreteFunctionSpace_.order();
    const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
    const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

    // do loop over all local DoFs
    for(  unsigned int i = 0;
          i < numberOfLocalDoFs;
          ++i )
    {
      // value of the functional, applied to the local basefunction
      RangeFieldType localFunctionalValue = 0.0;

      // get local basefunctions
      const LocalBaseFunctionType phi_i = localBaseFunctionProvider_.provide( entity, i );

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

        // apply local operation
        const double localOperationEvaluated = localOperation_.operate( phi_i, x );

        // compute integral
        localFunctionalValue += integrationFactor * quadratureWeight * localOperationEvaluated;

      } // done walk over quadrature points

      // set local vector
      ret[i] = localFunctionalValue;

    } // done loop over all local DoFs

    return ret;

  } // end applyLocal()

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const LocalOperationType localOperation_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class FiniteElement

template< class DiscreteFunctionSpaceImp, class LocalOperationImp >
class FiniteElementLOP
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef Dune::Functionals::LocalOperation::Interface< FunctionSpaceType, LocalOperationImp >
    LocalOperationType;

  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef Dune::Functionals::Common::LocalVector< RangeFieldType >
    LocalVectorType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef Dune::Functionals::Common::LocalBaseFunctionProvider< DiscreteFunctionSpaceType >
    LocalBaseFunctionProviderType;

  typedef typename LocalBaseFunctionProviderType::LocalBaseFunctionType
    LocalBaseFunctionType;

  FiniteElementLOP( const DiscreteFunctionSpaceType& discreteFunctionSpace,
                    const LocalOperationType& localOperation )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      localOperation_( localOperation ),
      localBaseFunctionProvider_( discreteFunctionSpace )
  {
  }

  /**
    * \brief  Destructor
    *         Does nothing.
    **/
  ~FiniteElementLOP()
  {
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

  const LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  const LocalVectorType applyLocal( const EntityType& entity ) const
  {
    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet( entity );
    const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();

    // init return vector
    LocalVectorType ret( numberOfLocalDoFs );

    // do loop over all local DoFs
    for(  unsigned int i = 0;
          i < numberOfLocalDoFs;
          ++i )
    {
      // get local basefunctions
      const LocalBaseFunctionType localBaseFunction_i = localBaseFunctionProvider_.provide( entity, i );

      // set local vector
      ret[i] = localOperation_.operateLocal( localBaseFunction_i );

    } // done loop over all local DoFs

    return ret;

  } // end applyLocal()

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const LocalOperationType localOperation_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class FiniteElementLOP

} // end namespace Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH
