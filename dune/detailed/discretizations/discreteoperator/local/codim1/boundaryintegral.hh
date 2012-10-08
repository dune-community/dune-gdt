/**
  \file   boundaryintegral.hh
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH

// dune-common
#include <dune/common/densematrix.hh>

// dune-geometry
#include <dune/geometry/quadraturerules.hh>

// dune-stuff
#include <dune/stuff/common/matrix.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace DiscreteOperator {

namespace Local {

namespace Codim1 {

/**
  \brief  Local operator for inner intersections, i.e. those who have an inner codim 0 entity (Entity or En) and an
          outer codim 0 neighbouring entity (Neighbour or Ne).
  **/
template< class LocalEvaluationImp >
class BoundaryIntegral
{
public:
  typedef LocalEvaluationImp LocalEvaluationType;

  typedef BoundaryIntegral< LocalEvaluationType > ThisType;

  typedef typename LocalEvaluationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

//  template< class InducingDiscreteFunctionType >
//  class LocalFunctional
//  {
//  public:
//    typedef Dune::Functionals::DiscreteFunctional::Local::Codim0::IntegralInduced<  ThisType,
//                                                                                    InducingDiscreteFunctionType >
//      Type;
//  };

  BoundaryIntegral(const LocalEvaluationType& localEvaluation)
    : localEvaluation_( localEvaluation )
  {}

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
  void applyLocal(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                  const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                  const IntersectionType& intersection,
                  LocalMatrixType& localMatrix,
                  std::vector< LocalMatrixType >& tmpLocalMatrices) const
  {
    // some stuff
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    const unsigned int cols = localTestBaseFunctionSet.size();

    // make sure, that the target matrix is big enough
    assert(localMatrix.rows() >= rows);
    assert(localMatrix.cols() >= cols);

    // clear target matrix
    Dune::Stuff::Common::Matrix::clear(localMatrix);

    // check tmp local matrices
    if (tmpLocalMatrices.size() < numTmpObjectsRequired())
    {
      tmpLocalMatrices.resize(numTmpObjectsRequired(),
                              LocalMatrixType(
                                localAnsatzBaseFunctionSet.baseFunctionSet().space().map().maxLocalSize(),
                                localTestBaseFunctionSet.baseFunctionSet().space().map().maxLocalSize(),
                                RangeFieldType(0.0)));
    } // check tmp local matrices

    // quadrature
    const unsigned int quadratureOrder = localEvaluation_.order()
      + localAnsatzBaseFunctionSet.order()
      + localTestBaseFunctionSet.order();
    typedef Dune::QuadratureRules< DomainFieldType, IntersectionType::mydimension > FaceQuadratureRules;
    typedef Dune::QuadratureRule< DomainFieldType, IntersectionType::mydimension > FaceQuadratureType;
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule(intersection.type(), 2*quadratureOrder + 1);

    // loop over all quadrature points
    for (typename FaceQuadratureType::const_iterator quadPoint = faceQuadrature.begin();
         quadPoint != faceQuadrature.end();
         ++quadPoint) {
      // local coordinates
      const  typename IntersectionType::LocalCoordinate x = quadPoint->position();

      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement( x );
      const double quadratureWeight = quadPoint->weight();

      // evaluate the local operation
      localEvaluation_.evaluateLocal(localAnsatzBaseFunctionSet,
                                     localTestBaseFunctionSet,
                                     intersection,
                                     x,
                                     tmpLocalMatrices[0]);

      // compute integral (see below)
      addToIntegral(tmpLocalMatrices[0],
                    integrationFactor,
                    quadratureWeight,
                    rows,
                    cols,
                    localMatrix);

    } // loop over all quadrature points
  } // end method applyLocal

private:
  BoundaryIntegral(const ThisType&);
  ThisType& operator=( const ThisType& );

  template< class MatrixImp, class RangeFieldType >
  void addToIntegral(const Dune::DenseMatrix< MatrixImp >& localMatrix,
                     const RangeFieldType& integrationFactor,
                     const RangeFieldType& quadratureWeight,
                     const unsigned int rows,
                     const unsigned int cols,
                     Dune::DenseMatrix< MatrixImp >& ret) const
  {
    // loop over all rows
    for(unsigned int i = 0; i < rows; ++i) {
      // get row
      const typename Dune::DenseMatrix< MatrixImp >::const_row_reference localMatrixRow = localMatrix[i];
      typename Dune::DenseMatrix< MatrixImp >::row_reference retRow = ret[i];
      // loop over all cols
      for(unsigned int j = 0; j < cols; ++j) {
        retRow[j] += localMatrixRow[j] * integrationFactor * quadratureWeight;
      } // loop over all cols
    } // loop over all rows
  } // void addToIntegral(...) const

  const LocalEvaluationType& localEvaluation_;
}; // end class BoundaryIntegral

} // end namespace Codim1

} // end namespace Local

} // end namespace DiscreteOperator

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_BOUNDARYINTEGRAL_HH
