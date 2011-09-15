/**
  \file   integral.hh
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH

// dune-fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune
{

namespace DetailedDiscretizations
{

namespace DiscreteOperator
{

namespace Local
{

namespace Codim1
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

//  template< class InducingDiscreteFunctionType >
//  class LocalFunctional
//  {
//  public:
//    typedef Dune::Functionals::DiscreteFunctional::Local::Codim0::IntegralInduced<  ThisType,
//                                                                                    InducingDiscreteFunctionType >
//      Type;
//  };

  Integral( const LocalEvaluationType localEvaluation )
    : localEvaluation_( localEvaluation )
  {
  }

  //! copy constructor
  Integral( const ThisType& other )
    : localEvaluation_( other.localEvaluation() )
  {
  }

  const LocalEvaluationType& localEvaluation() const
  {
    return localEvaluation_;
  }

//  template< class InducingDiscreteFunctionType >
//  typename LocalFunctional< InducingDiscreteFunctionType >::Type
//    localFunctional( const InducingDiscreteFunctionType& inducingDiscreteFunction ) const
//  {
//    typedef Dune::Functionals::DiscreteFunctional::Local::Codim0::IntegralInduced<  ThisType,
//                                                                                    InducingDiscreteFunctionType >
//      LocalFunctionalType;

//    return LocalFunctionalType( *this, inducingDiscreteFunction );
//  } // end method localFunctional

  template< class IntersectionType,
           class LocalAnsatzBaseFunctionSetType, /*Entity*/
           class LocalTestBaseFunctionSetType, /*Neighour*/
           class LocalMatrixType >
  void applyLocal( const IntersectionType& intersection,
                   const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                   const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                   LocalMatrixType& localMatrix ) const
  {
    // clear target matrix
    for( unsigned int i = 0; i < localMatrix.rows(); ++i )
    {
      for( unsigned int j = 0; j < localMatrix.cols(); ++j )
      {
        localMatrix[i][j] = 0.0;
      }
    }

    // some types
    typedef typename LocalAnsatzBaseFunctionSetType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    typedef Dune::CachingQuadrature< GridPartType, 1 >
      FaceQuadratureType;

    // some stuff
    const GridPartType& gridPart = localAnsatzBaseFunctionSet.space().gridPart();
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    const unsigned int cols = localTestBaseFunctionSet.size();
    const unsigned int quadratureOrder = localEvaluation_.order() + localAnsatzBaseFunctionSet.order() + localTestBaseFunctionSet.order();
    const FaceQuadratureType faceQuadrature(  gridPart,
                                              intersection,
                                              quadratureOrder,
                                              FaceQuadratureType::INSIDE );
    const unsigned int numberOfQuadraturePoints = faceQuadrature.nop();

    // some tmp storage
    LocalMatrixType tmpMatrix( rows, cols );

    // do loop over all quadrature points
    for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
    {
      // local coordinates
      const DomainType x = faceQuadrature.point( q );
      const DomainType xInEntity = intersection.geometryInInside().global( x );
      const DomainType xInNeighbour = intersection.geometryInOutside().global( x );
      const DomainType unitOuterNormal = intersection.unitOuterNormal();

      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement( x );
      const double quadratureWeight = faceQuadrature.weight( q );

      // evaluate the local operation
      localEvaluation_.evaluate(  localAnsatzBaseFunctionSet,
                                  localTestBaseFunctionSet,
                                  xInEntity,
                                  xInNeighbour,
                                  unitOuterNormal,
                                  tmpMatrix );

      // compute integral
      for( unsigned int i = 0; i < rows; ++i )
      {
        for( unsigned int j = 0; j < cols; ++j )
        {
          localMatrix[i][j] += tmpMatrix[i][j] * integrationFactor * quadratureWeight;
        }
      }
    } // done loop over all quadrature points

  } // end method applyLocal

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const LocalEvaluationType localEvaluation_;

}; // end class Codim0Integration

} // end namespace Codim1

} // end namespace Local

} // end namespace DiscreteOperator

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH
