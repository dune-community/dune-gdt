/**
  \file   boundaryintegral.hh
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH

// dune grid includes
#include <dune/grid/common/quadraturerules.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/matrix.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace DiscreteOperator
{

namespace Local
{

namespace Codim1
{

/**
  \brief  Local operator for inner intersections, i.e. those who have an inner codim 0 entity (Entity or En) and an
          outer codim 0 neighbouring entity (Neighbour or Ne).
  **/
template< class LocalEvaluationImp >
class BoundaryIntegral
{
public:

  typedef LocalEvaluationImp
    LocalEvaluationType;

  typedef BoundaryIntegral< LocalEvaluationType >
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

  BoundaryIntegral( const LocalEvaluationType localEvaluation )
    : localEvaluation_( localEvaluation )
  {
  }

  //! copy constructor
  BoundaryIntegral( const ThisType& other )
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

  unsigned int numTmpObjectsRequired() const
  {
    return 1;
  }

  template< class LocalAnsatzBaseFunctionSetType,
            class LocalTestBaseFunctionSetType,
            class IntersectionType,
            class LocalMatrixType >
  void applyLocal(  const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                    const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                    const IntersectionType& intersection,
                    LocalMatrixType& localMatrix,
                    std::vector< LocalMatrixType >& tmpLocalMatrices ) const
  {
    // some types
    typedef typename LocalAnsatzBaseFunctionSetType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridViewType
      GridViewType;

    typedef Dune::QuadratureRules< double, IntersectionType::mydimension >
      FaceQuadratureRules;

    typedef Dune::QuadratureRule< double, IntersectionType::mydimension >
      FaceQuadratureType;

    typedef typename IntersectionType::LocalCoordinate
      LocalCoordinateType;

    // some stuff
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    const unsigned int cols = localTestBaseFunctionSet.size();
    const unsigned int quadratureOrder =
      localEvaluation_.order() + localAnsatzBaseFunctionSet.order() + localTestBaseFunctionSet.order();
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule( intersection.type(), 2*quadratureOrder+1 );

    // make sure, that the target matrices are big enough
    assert( localMatrix.rows() >= rows );
    assert( localMatrix.cols() >= cols );

    // clear target matrices
    Dune::HelperTools::Common::Matrix::clear( localMatrix );

    // check tmp local matrices
    if( tmpLocalMatrices.size() < 1 )
    {
      tmpLocalMatrices.resize(  1,
                                LocalMatrixType(  localAnsatzBaseFunctionSet.baseFunctionSet().space().map().maxLocalSize(),
                                                  localTestBaseFunctionSet.baseFunctionSet().space().map().maxLocalSize(),
                                                  RangeFieldType( 0.0 ) ) );
    }

    const typename FaceQuadratureType::const_iterator quadratureEnd = faceQuadrature.end();
    for (typename FaceQuadratureType::const_iterator quadPoint = faceQuadrature.begin(); quadPoint != quadratureEnd; ++quadPoint )
    {
      // local coordinates
      const LocalCoordinateType x = quadPoint->position();

      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement( x );
      const double quadratureWeight = quadPoint->weight();

      // evaluate the local operation
      localEvaluation_.evaluateLocal( localAnsatzBaseFunctionSet,
                                      localTestBaseFunctionSet,
                                      intersection,
                                      x,
                                      tmpLocalMatrices[0] );

      // compute integral (see below)
      addToIntegral( tmpLocalMatrices[0], integrationFactor, quadratureWeight, rows, cols, localMatrix );

    } // done loop over all quadrature points
  } // end method applyLocal

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  template< class LocalMatrixType,
            class RangeFieldType >
  void addToIntegral( const LocalMatrixType& localMatrix,
                      const RangeFieldType& integrationFactor,
                      const RangeFieldType& quadratureWeight,
                      const unsigned int rows,
                      const unsigned int cols,
                      LocalMatrixType& ret ) const
  {
    for( unsigned int i = 0; i < rows; ++i )
    {
      for( unsigned int j = 0; j < cols; ++j )
      {
        ret[i][j] += localMatrix[i][j] * integrationFactor * quadratureWeight;
      }
    }
  } // end method addToIntegral

  const LocalEvaluationType localEvaluation_;

}; // end class BoundaryIntegral

} // end namespace Codim1

} // end namespace Local

} // end namespace DiscreteOperator

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH
