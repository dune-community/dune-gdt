/**
  \file   integral.hh
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH

// dune-fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/matrix.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace DiscreteOperator {

namespace Local {

namespace Codim1 {

template <class LocalEvaluationImp>
class Integral
{
public:
  typedef LocalEvaluationImp LocalEvaluationType;

  typedef Integral<LocalEvaluationType> ThisType;

  typedef typename LocalEvaluationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  //  template< class InducingDiscreteFunctionType >
  //  class LocalFunctional
  //  {
  //  public:
  //    typedef Dune::Functionals::DiscreteFunctional::Local::Codim0::IntegralInduced<  ThisType,
  //                                                                                    InducingDiscreteFunctionType >
  //      Type;
  //  };

  Integral(const LocalEvaluationType localEvaluation)
    : localEvaluation_(localEvaluation)
  {
  }

  //! copy constructor
  Integral(const ThisType& other)
    : localEvaluation_(other.localEvaluation())
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
    \todo Rename Entity -> En, Neighbour -> Ne
    **/
  template <class LocalAnsatzBaseFunctionSetEntityType, class LocalAnsatzBaseFunctionSetNeighbourType,
            class LocalTestBaseFunctionSetEntityType, class LocalTestBaseFunctionSetNeighbourType,
            class IntersectionType, class LocalMatrixType>
  void applyLocal(const LocalAnsatzBaseFunctionSetEntityType& localAnsatzBaseFunctionSetEntity,
                  const LocalAnsatzBaseFunctionSetNeighbourType& localAnsatzBaseFunctionSetNeighbour,
                  const LocalTestBaseFunctionSetEntityType& localTestBaseFunctionSetEntity,
                  const LocalTestBaseFunctionSetNeighbourType& localTestBaseFunctionSetNeighbour,
                  const IntersectionType& intersection, LocalMatrixType& localMatrixEnEn,
                  LocalMatrixType& localMatrixEnNe, LocalMatrixType& localMatrixNeEn, LocalMatrixType& localMatrixNeNe,
                  std::vector<LocalMatrixType>& tmpLocalMatrices) const
  {
    // some types
    typedef typename LocalAnsatzBaseFunctionSetEntityType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef Dune::CachingQuadrature<GridPartType, 1> FaceQuadratureType;

    // some stuff
    const GridPartType& gridPart = localAnsatzBaseFunctionSetEntity.space().gridPart();
    const unsigned int rowsEn    = localAnsatzBaseFunctionSetEntity.size();
    const unsigned int rowsNe    = localAnsatzBaseFunctionSetNeighbour.size();
    const unsigned int colsEn    = localTestBaseFunctionSetEntity.size();
    const unsigned int colsNe = localTestBaseFunctionSetNeighbour.size();
    const unsigned int quadratureOrder =
        localEvaluation_.order()
        + std::max(localAnsatzBaseFunctionSetEntity.order(), localAnsatzBaseFunctionSetNeighbour.order())
        + std::max(localTestBaseFunctionSetEntity.order(), localTestBaseFunctionSetNeighbour.order());
    const FaceQuadratureType faceQuadrature(gridPart, intersection, quadratureOrder, FaceQuadratureType::INSIDE);
    const unsigned int numberOfQuadraturePoints = faceQuadrature.nop();

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
    Dune::HelperTools::Common::Matrix::clear(localMatrixEnEn);
    Dune::HelperTools::Common::Matrix::clear(localMatrixEnNe);
    Dune::HelperTools::Common::Matrix::clear(localMatrixNeEn);
    Dune::HelperTools::Common::Matrix::clear(localMatrixNeNe);

    // check tmp local matrices
    if (tmpLocalMatrices.size() < 3) {
      tmpLocalMatrices.resize(
          3,
          LocalMatrixType(localAnsatzBaseFunctionSetEntity.baseFunctionSet().space().map().maxLocalSize(),
                          localTestBaseFunctionSetEntity.baseFunctionSet().space().map().maxLocalSize(),
                          RangeFieldType(0.0)));
    }

    // do loop over all quadrature points
    for (unsigned int q = 0; q < numberOfQuadraturePoints; ++q) {
      // local coordinates
      const DomainType x = faceQuadrature.point(q);

      // integration factors
      const double integrationFactor = intersection.geometry().integrationElement(x);
      const double quadratureWeight  = faceQuadrature.weight(q);

      // evaluate the local operation
      localEvaluation_.evaluate(localAnsatzBaseFunctionSetEntity,
                                localAnsatzBaseFunctionSetNeighbour,
                                localTestBaseFunctionSetEntity,
                                localTestBaseFunctionSetNeighbour,
                                intersection,
                                x,
                                tmpLocalMatrices[0], /*EnEn*/
                                tmpLocalMatrices[1], /*EnNe*/
                                tmpLocalMatrices[2], /*NeEn*/
                                tmpLocalMatrices[3]); /*NeNe*/

      // compute integral
      addToIntegral(localMatrixEnEn, integrationFactor, quadratureWeight, rowsEn, colsEn, tmpLocalMatrices[0]);
      addToIntegral(localMatrixEnNe, integrationFactor, quadratureWeight, rowsEn, colsNe, tmpLocalMatrices[1]);
      addToIntegral(localMatrixNeEn, integrationFactor, quadratureWeight, rowsNe, colsEn, tmpLocalMatrices[2]);
      addToIntegral(localMatrixNeNe, integrationFactor, quadratureWeight, rowsNe, colsNe, tmpLocalMatrices[3]);

    } // done loop over all quadrature points

  } // end method applyLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class LocalMatrixType, class RangeFieldType>
  void addToIntegral(const LocalMatrixType& localMatrix, const RangeFieldType& integrationFactor,
                     const RangeFieldType& quadratureWeight, const unsigned int rows, const unsigned int cols,
                     LocalMatrixType& ret) const
  {
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
        ret[i][j] += localMatrix[i][j] * integrationFactor * quadratureWeight;
      }
    }
  } // end method addToIntegral

  const LocalEvaluationType localEvaluation_;

}; // end class Codim0Integration

} // end namespace Codim1

} // end namespace Local

} // end namespace DiscreteOperator

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // end DUNE_DETAILED_DISCRETIZATIONS_DISCRETEOPERATOR_LOCAL_CODIM1_INTEGRAL_HH
