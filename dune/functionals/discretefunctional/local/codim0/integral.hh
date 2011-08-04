/**
  \file   integration.hh
  **/

#ifndef DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH
#define DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-functionals includes
#include <dune/functionals/common/localvector.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctional {

namespace Local {

namespace Codim0 {

template <class LocalEvaluationImp>
class Integral
{
public:
  typedef LocalEvaluationImp LocalEvaluationType;

  typedef typename LocalEvaluationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  Integral(const LocalEvaluationType& localEvaluation)
    : localEvaluation_(localEvaluation)
  {
  }

  template <class LocalTestBaseFunctionSetType, class LocalVectorType>
  void applyLocal(const LocalTestBaseFunctionSetType localTestBaseFunctionSet, LocalVectorType& localVector) const
  {
    // some types
    typedef typename LocalTestBaseFunctionSetType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef Dune::CachingQuadrature<GridPartType, 0> VolumeQuadratureType;

    //    typedef typename LocalTestBaseFunctionSetType::LocalBaseFunctionType
    //      LocalTestBaseFunctionType;

    // some stuff
    //    const unsigned numberOfLocalTestDoFs = localTestBaseFunctionSet.numBaseFunctions();
    const unsigned int quadratureOrder = 1 + localTestBaseFunctionSet.order();
    const VolumeQuadratureType volumeQuadrature(localTestBaseFunctionSet.entity(), quadratureOrder);
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    //    // do loop over all local test DoFs
    //    for( unsigned int j = 0; j < numberOfLocalTestDoFs; ++j )
    //    {
    //      const LocalTestBaseFunctionType localTestBaseFunction_j = localTestBaseFunctionSet.baseFunction( j );

    // do loop over all quadrature points
    LocalVectorType tmpVector(localTestBaseFunctionSet.size());
    for (unsigned int q = 0; q < numberOfQuadraturePoints; ++q) {
      // local coordinate
      const DomainType x = volumeQuadrature.point(q);

      // integration factors
      const double integrationFactor = localTestBaseFunctionSet.entity().geometry().integrationElement(x);
      const double quadratureWeight  = volumeQuadrature.weight(q);

      // evaluate the local evaluation
      localEvaluation_.evaluate(localTestBaseFunctionSet, x, tmpVector);

      // compute integral
      tmpVector *= integrationFactor * quadratureWeight;
      localVector += tmpVector;
    } // done loop over all quadrature points

    //    // set local vector (the = is important, since we dont assume a clean vector)
    //    localVector[j] = functional_j;

    //    } // done loop over all local test DoFs

  } // end method applyLocal

private:
  const LocalEvaluationType& localEvaluation_;

}; // end class Integral

// template< class InducingOperatorImp, class InducingDiscreteFunctionImp >
// class IntegralInduced
//{
// public:

//  typedef InducingOperatorImp
//    InducingOperatorType;

//  typedef InducingDiscreteFunctionImp
//    InducingDiscreteFunctionType;

//  typedef IntegralInduced< InducingOperatorType, InducingDiscreteFunctionType >
//    ThisType;

//  typedef typename InducingDiscreteFunctionType::RangeFieldType
//    RangeFieldType;

//  typedef typename InducingDiscreteFunctionType::DomainType
//    DomainType;

//  typedef Dune::Functionals::Common::LocalVector< RangeFieldType >
//    LocalVectorType;

//  IntegralInduced(  const InducingOperatorType& inducingOperator,
//                    const InducingDiscreteFunctionType& inducingDiscreteFunction )
//    : inducingOperator_( inducingOperator ),
//      inducingDiscreteFunction_( inducingDiscreteFunction )
//  {
//  }

//  IntegralInduced( const ThisType& other )
//    : inducingOperator_( other.inducingOperator() ),
//      inducingDiscreteFunction_( other.inducingDiscreteFunction() )
//  {
//  }

//  const InducingOperatorType& inducingOperator() const
//  {
//    return inducingOperator_;
//  }

//  const InducingDiscreteFunctionType& inducingDiscreteFunction() const
//  {
//    return inducingDiscreteFunction_;
//  }

//  template< class LocalTestBaseFunctionSetType >
//  void applyLocal( const LocalTestBaseFunctionSetType localTestBaseFunctionSet,
//                   LocalVectorType& localVector ) const
//  {
//    // some types
//    typedef typename LocalTestBaseFunctionSetType::DiscreteFunctionSpaceType
//      DiscreteFunctionSpaceType;

//    typedef typename DiscreteFunctionSpaceType::GridPartType
//      GridPartType;

//    typedef Dune::CachingQuadrature< GridPartType, 0 >
//      VolumeQuadratureType;

//    typedef typename LocalTestBaseFunctionSetType::LocalBaseFunctionType
//      LocalTestBaseFunctionType;

//    typedef typename LocalTestBaseFunctionType::EntityType
//      EntityType;

//    typedef typename InducingDiscreteFunctionType::ConstLocalFunctionType
//      InducingLocalFunctionType;

//    const EntityType& entity = localTestBaseFunctionSet.entity();

//    const InducingLocalFunctionType inducingLocalFunction = inducingDiscreteFunction_.localFunction( entity );

//    // some stuff
//    const unsigned numberOfLocalTestDoFs = localTestBaseFunctionSet.numBaseFunctions();
//    const unsigned int quadratureOrder = 1 + localTestBaseFunctionSet.order();
//    const VolumeQuadratureType volumeQuadrature( entity, quadratureOrder );
//    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

//    // do loop over all local test DoFs
//    for( unsigned int j = 0; j < numberOfLocalTestDoFs; ++j )
//    {
//      const LocalTestBaseFunctionType localTestBaseFunction_j = localTestBaseFunctionSet.baseFunction( j );

//      // do loop over all quadrature points
//      RangeFieldType functional_j( 0.0 );
//      for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
//      {
//        // local coordinate
//        const DomainType x = volumeQuadrature.point( q );

//        // integration factors
//        const double integrationFactor = entity.geometry().integrationElement( x );
//        const double quadratureWeight = volumeQuadrature.weight( q );

//        // evaluate the local evaluation
//        const RangeFieldType localEvaluationEvalauted = inducingOperator_.localEvaluation().evaluate(
//        inducingLocalFunction, localTestBaseFunction_j, x );

//        // compute integral
//        functional_j += integrationFactor * quadratureWeight * localEvaluationEvalauted;
//      } // done loop over all quadrature points

//      // set local vector (the = is important, since we dont assume a clean vector)
//      localVector[j] = functional_j;

//    } // done loop over all local test DoFs

//  } // end method applyLocal

// private:

//  const InducingOperatorType& inducingOperator_;
//  const InducingDiscreteFunctionType& inducingDiscreteFunction_;

//}; // end class IntegralInduced

} // end namespace Codim0

} // end namespace Local

} // end namespace DiscreteFunctional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_LOCAL_INTEGRATION_HH
