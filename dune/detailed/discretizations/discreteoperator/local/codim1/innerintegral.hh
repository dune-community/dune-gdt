/**
  \file   innerintegral.hh
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INNERINTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INNERINTEGRAL_HH

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
          outer codim 0 Neighboring entity (Neighbor or Ne).
  **/
template <class LocalEvaluationImp>
class InnerIntegral
{
public:
  typedef LocalEvaluationImp LocalEvaluationType;

  typedef InnerIntegral<LocalEvaluationType> ThisType;

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

  InnerIntegral(const LocalEvaluationType& localEvaluation)
    : localEvaluation_(localEvaluation)
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
    return 4;
  }

  /**
    \todo Rename Entity -> En, Neighbor -> Ne
    **/
  template <class LocalAnsatzBaseFunctionSetEntityType, class LocalTestBaseFunctionSetEntityType,
            class LocalAnsatzBaseFunctionSetNeighborType, class LocalTestBaseFunctionSetNeighborType,
            class IntersectionType, class LocalMatrixType>
  void applyLocal(const LocalAnsatzBaseFunctionSetEntityType& localAnsatzBaseFunctionSetEntity,
                  const LocalTestBaseFunctionSetEntityType& localTestBaseFunctionSetEntity,
                  const LocalAnsatzBaseFunctionSetNeighborType& localAnsatzBaseFunctionSetNeighbor,
                  const LocalTestBaseFunctionSetNeighborType& localTestBaseFunctionSetNeighbor,
                  const IntersectionType& intersection, LocalMatrixType& localMatrixEnEn,
                  LocalMatrixType& localMatrixNeNe, LocalMatrixType& localMatrixEnNe, LocalMatrixType& localMatrixNeEn,
                  std::vector<LocalMatrixType>& tmpLocalMatrices) const
  {
    // preparations
    const unsigned int rowsEn = localAnsatzBaseFunctionSetEntity.size();
    const unsigned int rowsNe = localAnsatzBaseFunctionSetNeighbor.size();
    const unsigned int colsEn = localTestBaseFunctionSetEntity.size();
    const unsigned int colsNe = localTestBaseFunctionSetNeighbor.size();

    // make sure, that the target matrices are big enough
    assert(localMatrixEnEn.rows() >= rowsEn);
    assert(localMatrixEnEn.cols() >= colsEn);
    assert(localMatrixEnNe.rows() >= rowsEn);
    assert(localMatrixEnNe.cols() >= colsNe);
    assert(localMatrixNeEn.rows() >= rowsNe);
    assert(localMatrixNeEn.cols() >= colsEn);
    assert(localMatrixNeNe.rows() >= rowsNe);
    assert(localMatrixNeNe.cols() >= colsNe);

    // clear target matrices
    Dune::Stuff::Common::Matrix::clear(localMatrixEnEn);
    Dune::Stuff::Common::Matrix::clear(localMatrixEnNe);
    Dune::Stuff::Common::Matrix::clear(localMatrixNeEn);
    Dune::Stuff::Common::Matrix::clear(localMatrixNeNe);

    // check tmp local matrices
    if (tmpLocalMatrices.size() < numTmpObjectsRequired()) {
      tmpLocalMatrices.resize(
          numTmpObjectsRequired(),
          LocalMatrixType(std::max(localAnsatzBaseFunctionSetEntity.baseFunctionSet().space().map().maxLocalSize(),
                                   localAnsatzBaseFunctionSetNeighbor.baseFunctionSet().space().map().maxLocalSize()),
                          std::max(localTestBaseFunctionSetEntity.baseFunctionSet().space().map().maxLocalSize(),
                                   localTestBaseFunctionSetNeighbor.baseFunctionSet().space().map().maxLocalSize()),
                          RangeFieldType(0.0)));
    } // check tmp local matrices

    // quadrature
    const unsigned int quadratureOrder =
        localEvaluation_.order()
        + std::max(localAnsatzBaseFunctionSetEntity.order(), localAnsatzBaseFunctionSetNeighbor.order())
        + std::max(localTestBaseFunctionSetEntity.order(), localTestBaseFunctionSetNeighbor.order());
    typedef Dune::QuadratureRules<DomainFieldType, IntersectionType::mydimension> FaceQuadratureRules;
    typedef Dune::QuadratureRule<DomainFieldType, IntersectionType::mydimension> FaceQuadratureType;
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule(intersection.type(), 2 * quadratureOrder + 1);

    // do loop over all quadrature points
    for (typename FaceQuadratureType::const_iterator quadPoint = faceQuadrature.begin();
         quadPoint != faceQuadrature.end();
         ++quadPoint) {
      // local coordinates
      typedef typename IntersectionType::LocalCoordinate LocalCoordinateType;
      const LocalCoordinateType x = quadPoint->position();

      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement(x);
      const double quadratureWeight  = quadPoint->weight();

      // evaluate the local operation
      localEvaluation_.evaluateLocal(localAnsatzBaseFunctionSetEntity,
                                     localTestBaseFunctionSetEntity,
                                     localAnsatzBaseFunctionSetNeighbor,
                                     localTestBaseFunctionSetNeighbor,
                                     intersection,
                                     x,
                                     tmpLocalMatrices[0], /*EnEn*/
                                     tmpLocalMatrices[1], /*NeNe*/
                                     tmpLocalMatrices[2], /*EnNe*/
                                     tmpLocalMatrices[3]); /*NeEn*/

      // compute integral (see below)
      addToIntegral(tmpLocalMatrices[0], integrationFactor, quadratureWeight, rowsEn, colsEn, localMatrixEnEn);
      addToIntegral(tmpLocalMatrices[1], integrationFactor, quadratureWeight, rowsNe, colsNe, localMatrixNeNe);
      addToIntegral(tmpLocalMatrices[2], integrationFactor, quadratureWeight, rowsEn, colsNe, localMatrixEnNe);
      addToIntegral(tmpLocalMatrices[3], integrationFactor, quadratureWeight, rowsNe, colsEn, localMatrixNeEn);
    } // done loop over all quadrature points
  } // void applyLocal

private:
  InnerIntegral(const ThisType& other);
  ThisType& operator=(const ThisType&);

  template <class MatrixImp, class RangeFieldType>
  void addToIntegral(const Dune::DenseMatrix<MatrixImp>& localMatrix, const RangeFieldType& integrationFactor,
                     const RangeFieldType& quadratureWeight, const unsigned int rows, const unsigned int cols,
                     Dune::DenseMatrix<MatrixImp>& ret) const
  {
    // loop over all rows
    for (unsigned int i = 0; i < rows; ++i) {
      // get row
      const typename Dune::DenseMatrix<MatrixImp>::const_row_reference localMatrixRow = localMatrix[i];
      typename Dune::DenseMatrix<MatrixImp>::row_reference retRow                     = ret[i];
      // loop over all cols
      for (unsigned int j = 0; j < cols; ++j) {
        retRow[j] += localMatrixRow[j] * integrationFactor * quadratureWeight;
      } // loop over all cols
    } // loop over all rows
  } // end method addToIntegral

  const LocalEvaluationType& localEvaluation_;
}; // end class InnerIntegral

} // end namespace Codim1

} // end namespace Local

} // end namespace DiscreteOperator

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INNERINTEGRAL_HH
