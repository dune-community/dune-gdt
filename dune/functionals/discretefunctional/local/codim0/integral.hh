/**
  \file   integration.hh
  **/

#ifndef DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH
#define DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH

// dune-common includes
#include <dune/common/dynmatrix.hh>

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-functionals includes
#include <dune/functionals/basefunctionset/local/wrapper.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctional
{

namespace Local
{

namespace Codim0
{

template< class LocalEvaluationImp >
class Integral
{
public:

  typedef LocalEvaluationImp
    LocalEvaluationType;

  typedef Integral< LocalEvaluationType >
    ThisType;

  typedef typename LocalEvaluationType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  Integral( const LocalEvaluationType localEvaluation )
    : localEvaluation_( localEvaluation )
  {
  }

//private:
  //! copy constructor
  Integral( const ThisType& other )
    : localEvaluation_( other.localEvaluation() )
  {
  }

//public:
  LocalEvaluationType localEvaluation() const
  {
    return localEvaluation_;
  }

  template< class LocalTestBaseFunctionSetType, class LocalVectorType >
  void applyLocal( const LocalTestBaseFunctionSetType localTestBaseFunctionSet,
                   LocalVectorType& ret ) const
  {
    // clear target
    for( unsigned int i = 0; i < ret.size(); ++i )
    {
      ret[i] = 0.0;
    }

    // some types
    typedef typename LocalTestBaseFunctionSetType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    typedef Dune::CachingQuadrature< GridPartType, 0 >
      VolumeQuadratureType;

    // some stuff
    const unsigned int size = localTestBaseFunctionSet.size();
    const unsigned int quadratureOrder = localEvaluation_.order() + localTestBaseFunctionSet.order();
    const VolumeQuadratureType volumeQuadrature( localTestBaseFunctionSet.entity(), quadratureOrder );
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    // do loop over all quadrature points
    LocalVectorType tmpVector( size );
    for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
    {
      // local coordinate
      const DomainType x = volumeQuadrature.point( q );

      // integration factors
      const double integrationFactor = localTestBaseFunctionSet.entity().geometry().integrationElement( x );
      const double quadratureWeight = volumeQuadrature.weight( q );

      // evaluate the local evaluation
      localEvaluation_.evaluate( localTestBaseFunctionSet, x, tmpVector );

      // compute integral
      for( unsigned int i = 0; i < size; ++i )
      {
        ret[i] += tmpVector[i] * integrationFactor * quadratureWeight;
      }
    } // done loop over all quadrature points
  } // end method applyLocal

private:

  //! assignment operator
  ThisType& operator=( const ThisType& other );

  const LocalEvaluationType localEvaluation_;

}; // end class Integral

template< class InducingOperatorImp, class InducingDiscreteFunctionImp >
class IntegralInduced
{
public:

  typedef InducingOperatorImp
    InducingOperatorType;

  typedef InducingDiscreteFunctionImp
    InducingDiscreteFunctionType;

  typedef IntegralInduced< InducingOperatorType, InducingDiscreteFunctionType >
    ThisType;

  typedef typename InducingDiscreteFunctionType::RangeFieldType
    RangeFieldType;

  typedef typename InducingDiscreteFunctionType::DomainType
    DomainType;

  IntegralInduced(  const InducingOperatorType& inducingOperator,
                    const InducingDiscreteFunctionType& inducingDiscreteFunction )
    : inducingOperator_( inducingOperator ),
      inducingDiscreteFunction_( inducingDiscreteFunction )
  {
  }

  //! copy constructor
  IntegralInduced( const ThisType& other )
    : inducingOperator_( other.inducingOperator() ),
      inducingDiscreteFunction_( other.inducingDiscreteFunction() )
  {
  }

  InducingOperatorType inducingOperator() const
  {
    return inducingOperator_;
  }

  InducingDiscreteFunctionType inducingDiscreteFunction() const
  {
    return inducingDiscreteFunction_;
  }

  template< class LocalTestBaseFunctionSetType, class LocalVectorType >
  void applyLocal( const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                   LocalVectorType& ret ) const
  {
    // clear target
    for( unsigned int i = 0; i < ret.size(); ++i )
    {
      ret[i] = 0.0;
    }

    // some types
    typedef typename LocalTestBaseFunctionSetType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    typedef Dune::CachingQuadrature< GridPartType, 0 >
      VolumeQuadratureType;

    typedef typename LocalTestBaseFunctionSetType::EntityType
      EntityType;

    typedef typename InducingDiscreteFunctionType::ConstLocalFunctionType
      InducingLocalFunctionType;

    typedef typename Dune::Functionals::BaseFunctionSet::Local::Wrapper< InducingLocalFunctionType >
      InducingBaseFunctionSetType;

    typedef Dune::DynamicMatrix< RangeFieldType >
      LocalMatrixType;

    const EntityType& entity = localTestBaseFunctionSet.entity();

    // wrap inducing local function
    const InducingLocalFunctionType inducingLocalFunction = inducingDiscreteFunction_.localFunction( entity );
    const InducingBaseFunctionSetType inducingBaseFunctionSet( inducingLocalFunction );

    // some stuff
    const unsigned int size = localTestBaseFunctionSet.size();
    const unsigned int quadratureOrder =
      inducingOperator_.localEvaluation().order() + inducingLocalFunction.order() + localTestBaseFunctionSet.order();
    const VolumeQuadratureType volumeQuadrature( entity, quadratureOrder );
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    // do loop over all quadrature points
    LocalMatrixType tmpMatrix( 1, size );
    for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
    {
      // local coordinate
      const DomainType x = volumeQuadrature.point( q );

      // integration factors
      const double integrationFactor = entity.geometry().integrationElement( x );
      const double quadratureWeight = volumeQuadrature.weight( q );

      // evaluate the local evaluation
      inducingOperator_.localEvaluation().evaluate( inducingBaseFunctionSet, localTestBaseFunctionSet, x, tmpMatrix );

      // compute integral
      for( unsigned int i = 0; i < size; ++i )
      {
        ret[i] += tmpMatrix[0][i] * integrationFactor * quadratureWeight;
      }
    } // done loop over all quadrature points

  } // end method applyLocal

private:
  //! assignment operator
  ThisType& operator=( const ThisType& );

  const InducingOperatorType inducingOperator_;
  const InducingDiscreteFunctionType inducingDiscreteFunction_;

}; // end class IntegralInduced

} // end namespace Codim0

} // end namespace Local

} // end namespace DiscreteFunctional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH
