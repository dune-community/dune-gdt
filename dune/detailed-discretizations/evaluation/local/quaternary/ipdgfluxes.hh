#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH

namespace Dune
{

namespace DetailedDiscretizations
{

namespace Evaluation
{

namespace Local
{

namespace Quaternary
{

template< class FunctionSpaceImp >
class IPDGFluxesInner
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef IPDGFluxesInner< FunctionSpaceType >
    ThisType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef Dune::FemTools::Function::Runtime< FunctionSpaceType >
    InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  IPDGFluxesInner( const std::string expression = "[1.0;1.0;1.0]", const int order = 1 )
    : inducingFunction_( expression ),
      order_( std::max( 0, order ) )
  {
  }

  //! copy constructor
  IPDGFluxesInner( const ThisType& other )
    : inducingFunction_( other.inducingFunction() ),
      order_( other.order() )
  {
  }

  //! returns the inducing function
  InducingFunctionType inducingFunction() const
  {
    return inducingFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    \attention  Assumes, that the return matrices are empty, because we do multiple += !
    **/
  template< class LocalAnsatzBaseFunctionSetEntityType,
            class LocalAnsatzBaseFunctionSetNeighbourType,
            class LocalTestBaseFunctionSetEntityType,
            class LocalTestBaseFunctionSetNeighbourType,
            class IntersectionType,
            class LocalMatrixType >
  void evaluateLocal( const LocalAnsatzBaseFunctionSetEntityType& localAnsatzBaseFunctionSetEntity,
                      const LocalAnsatzBaseFunctionSetNeighbourType& localAnsatzBaseFunctionSetNeighbour,
                      const LocalTestBaseFunctionSetEntityType& localTestBaseFunctionSetEntity,
                      const LocalTestBaseFunctionSetNeighbourType& localTestBaseFunctionSetNeighbour,
                      const IntersectionType& intersection,
                      const DomainType& localPointEntity,
                      const DomainType& localPointNeighbour,
                      LocalMatrixType& entityEntityRet,
                      LocalMatrixType& entityNeighbourRet,
                      LocalMatrixType& neighbourEntityRet,
                      LocalMatrixType& neighbourNeighbourRet ) const
  {
    // some stuff
    const DomainType globalPoint = localAnsatzBaseFunctionSet.entity().geometry().global( localPointEntity );
    const DomainType unitOuterNormal = intersection.unitOuterNormal();

    // evaluate ansatz entity basefunctionset
    const unsigned int rowsEntity = localAnsatzBaseFunctionSetEntity.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetEntityEvaluations( rowsEntity, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetEntityGradients( rowsEntity, JacobianRangeType( 0.0 ) );
    localAnsatzBaseFunctionSetEntity.evaluate( localPointEntity, localAnsatzBaseFunctionSetEntityEvaluations );
    localAnsatzBaseFunctionSetEntity.jacobian( localPointEntity, localAnsatzBaseFunctionSetEntityGradients );

    // evaluate ansatz neighbour basefunctionset
    const unsigned int rowsNeighbour = localAnsatzBaseFunctionSetNeighbour.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetNeighbourEvaluations( rowsNeighbour, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetNeighbourGradients( rowsNeighbour, JacobianRangeType( 0.0 ) );
    localAnsatzBaseFunctionSetNeighbour.evaluate( localPointEntity, localAnsatzBaseFunctionSetNeighbourEvaluations );
    localAnsatzBaseFunctionSetNeighbour.jacobian( localPointEntity, localAnsatzBaseFunctionSetNeighbourGradients );

    // evaluate test entity basefunctionset
    const unsigned int colsEntity = localTestBaseFunctionSetEntity.size();
    std::vector< RangeType > localTestBaseFunctionSetEntityEvaluations( colsEntity, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localTestBaseFunctionSetEntityGradients( colsEntity, JacobianRangeType( 0.0 ) );
    localTestBaseFunctionSetEntity.evaluate( localPointEntity, localTestBaseFunctionSetEntityEvaluations );
    localTestBaseFunctionSetEntity.jacobian( localPointEntity, localTestBaseFunctionSetEntityGradients );

    // evaluate test entity basefunctionset
    const unsigned int colsNeighbour = localTestBaseFunctionSetNeighbour.size();
    std::vector< RangeType > localTestBaseFunctionSetNeighbourEvaluations( colsNeighbour, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localTestBaseFunctionSetNeighbourGradients( colsNeighbour, JacobianRangeType( 0.0 ) );
    localTestBaseFunctionSetNeighbour.evaluate( localPointNeighbour, localTestBaseFunctionSetNeighbourEvaluations );
    localTestBaseFunctionSetNeighbour.jacobian( localPointNeighbour, localTestBaseFunctionSetNeighbourGradients );

    // evaluate inducing function
    RangeType functionValue( 0.0 );
    inducingFunction_.evaluate( globalPoint, functionValue );

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = 1.0 / std::pow( intersection.geometry().volume(), 1.0 );

    // entity entity combinations
    // do loop over all ansatz and test basefunctions
    assert( entityEntityRet.rows() == rowsEntity );
    assert( entityEntityRet.cols() == colsEntity );
    for( unsigned int i = 0; i < rowsEntity; ++i )
    {
      for( unsigned int j = 0; j < colsEntity; ++j )
      {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityEntityRet[i][j] += -0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evalautionTimesNormalTimesGradient = localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          entityEntityRet[i][j] += -0.5 * functionValue * evalautionTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evalautionTimesEvaluation = localAnsatzBaseFunctionSetEntityEvaluations[i] * localTestBaseFunctionSetEntityEvaluations[j];
          entityEntityRet[i][j] += penaltyParameter * evalautionTimesEvaluation;
        }
      }
    } // done entity entity combinations

    // do entity neighbour combinations
    assert( entityNeighbourRet.rows() == rowsEntity );
    assert( entityNeighbourRet.cols() == colsNeighbour );
    for( unsigned int i = 0; i < rowsEntity; ++i )
    {
      for( unsigned int j = 0; j < colsNeighbour; ++j )
      {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetNeighbourGradients[j][0] * unitOuterNormal;
          const RangeType evaluationTimesGradientTimesNormal = localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityNeighbourRet[i][j] += -0.5 * functionValue * evaluationTimesGradientTimesNormal;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetNeighbourEvaluations[j] * normalTimesGradient;
          entityNeighbourRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation = localTestBaseFunctionSetNeighbourEvaluations[j] * localAnsatzBaseFunctionSetEntityEvaluations[i];
          entityNeighbourRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done entity neighbour combinations

    // do neighbour entity combinations
    assert( neighbourEntityRet.rows() == rowsNeighbour );
    assert( neighbourEntityRet.cols() == colsEntity );
    for( unsigned int i = 0; i < rowsNeighbour; ++i )
    {
      for( unsigned int j = 0; j < colsEntity; ++j )
      {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetNeighbourEvaluations[i] * gradientTimesNormal;
          neighbourEntityRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetNeighbourGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          neighbourEntityRet[i][j] += -0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvalaution = localTestBaseFunctionSetEntityEvaluations[j] * localAnsatzBaseFunctionSetNeighbourEvaluations[i];
          neighbourEntityRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvalaution;
        }
      }
    } // done neighbour entity combinations

    // do neighbour neighbour combinations
    assert( neighbourNeighbourRet.rows() == rowsNeighbour );
    assert( neighbourNeighbourRet.cols() == colsNeighbour );
    for( unsigned int i = 0; i < rowsNeighbour; ++i )
    {
      for( unsigned int j = 0; j < colsNeighbour; ++j )
      {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetNeighbourGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvalaution = localAnsatzBaseFunctionSetNeighbourEvaluations[i] * gradientTimesNormal;
          neighbourNeighbourRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvalaution;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetNeighbourGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetNeighbourEvaluations[j] * normalTimesGradient;
          neighbourNeighbourRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation = localTestBaseFunctionSetNeighbourEvaluations[j] * localAnsatzBaseFunctionSetNeighbourEvaluations[i];
          neighbourNeighbourRet[i][j] += penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done neighbour neighbour combinations

  } // end method evaluateLocal


private:
  //! assignment operator
  ThisType& operator=( const ThisType& );

  const InducingFunctionType inducingFunction_;
  unsigned int order_;
}; // end class IPDGFluxesInner

} // end namespace Quaternary

} // end namespace Local

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
