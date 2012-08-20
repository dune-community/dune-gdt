#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/densematrix.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace Evaluation
{

namespace Local
{

namespace Quaternary
{

namespace IPDGfluxes
{

/**
  \todo       Implement penalty parameter
  \todo       Implement different constructor, for function and discretefunction
  **/
template< class FunctionSpaceImp,
          class InducingFunctionImp = Dune::Stuff::Function::Expression< typename FunctionSpaceImp::DomainFieldType, FunctionSpaceImp::DimDomain, typename FunctionSpaceImp::RangeFieldType, FunctionSpaceImp::DimRange > >
class Inner
{
public:
  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef InducingFunctionImp InducingFunctionType;

  typedef Inner< FunctionSpaceType, InducingFunctionType >
    ThisType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  Inner(const Dune::shared_ptr< const InducingFunctionType > inducingFunction,
        const unsigned int order,
        const RangeFieldType penaltyFactor)
    : inducingFunction_(inducingFunction)
    , order_(order)
    , penaltyFactor_(penaltyFactor)
  {}

  //! returns the inducing function
  const Dune::shared_ptr< const InducingFunctionType > inducingFunction() const
  {
    return inducingFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    \attention  Assumes, that the return matrices are empty, because we do multiple +=
    \todo       There are currently four nested loops, these can be simplyfied to two (see msdg-matlab)
    \todo       Rename Entity -> En, Neighbor -> Ne (otherwise stuff is too long)
    **/
  template< class LocalAnsatzBaseFunctionSetEntityType,
            class LocalTestBaseFunctionSetEntityType,
            class LocalAnsatzBaseFunctionSetNeighborType,
            class LocalTestBaseFunctionSetNeighborType,
            class IntersectionType,
            class LocalPointType,
            class LocalMatrixImp >
  void evaluateLocal( const LocalAnsatzBaseFunctionSetEntityType& localAnsatzBaseFunctionSetEntity,
                      const LocalTestBaseFunctionSetEntityType& localTestBaseFunctionSetEntity,
                      const LocalAnsatzBaseFunctionSetNeighborType& localAnsatzBaseFunctionSetNeighbor,
                      const LocalTestBaseFunctionSetNeighborType& localTestBaseFunctionSetNeighbor,
                      const IntersectionType& intersection,
                      const LocalPointType& localPoint,
                      Dune::DenseMatrix< LocalMatrixImp >& entityEntityRet,
                      Dune::DenseMatrix< LocalMatrixImp >& NeighborNeighborRet,
                      Dune::DenseMatrix< LocalMatrixImp >& entityNeighborRet,
                      Dune::DenseMatrix< LocalMatrixImp >& NeighborEntityRet) const
  {
    // some stuff
    const DomainType globalPoint = intersection.geometry().global(localPoint);
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType localPointNeighbor = intersection.geometryInOutside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    const RangeType zeroRange(0.0);
    const JacobianRangeType zeroJacobainRange(0.0);
    RangeFieldType gradientTimesNormal(0.0);
    RangeType gradientTimesNormalTimesEvaluation(0.0);

    // evaluate ansatz entity basefunctionset
    const unsigned int rowsEntity = localAnsatzBaseFunctionSetEntity.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetEntityEvaluations(rowsEntity, zeroRange);
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetEntityGradients(rowsEntity, zeroJacobainRange);
    localAnsatzBaseFunctionSetEntity.evaluate(localPointEntity, localAnsatzBaseFunctionSetEntityEvaluations);
    localAnsatzBaseFunctionSetEntity.jacobian(localPointEntity, localAnsatzBaseFunctionSetEntityGradients);

    // evaluate ansatz Neighbor basefunctionset
    const unsigned int rowsNeighbor = localAnsatzBaseFunctionSetNeighbor.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetNeighborEvaluations(rowsNeighbor, zeroRange);
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetNeighborGradients(rowsNeighbor, zeroJacobainRange);
    localAnsatzBaseFunctionSetNeighbor.evaluate(localPointNeighbor, localAnsatzBaseFunctionSetNeighborEvaluations);
    localAnsatzBaseFunctionSetNeighbor.jacobian(localPointNeighbor, localAnsatzBaseFunctionSetNeighborGradients);

    // evaluate test entity basefunctionset
    const unsigned int colsEntity = localTestBaseFunctionSetEntity.size();
    std::vector< RangeType > localTestBaseFunctionSetEntityEvaluations(colsEntity, zeroRange);
    std::vector< JacobianRangeType > localTestBaseFunctionSetEntityGradients(colsEntity, zeroJacobainRange);
    localTestBaseFunctionSetEntity.evaluate(localPointEntity, localTestBaseFunctionSetEntityEvaluations);
    localTestBaseFunctionSetEntity.jacobian(localPointEntity, localTestBaseFunctionSetEntityGradients);

    // evaluate test Neighbor basefunctionset
    const unsigned int colsNeighbor = localTestBaseFunctionSetNeighbor.size();
    std::vector< RangeType > localTestBaseFunctionSetNeighborEvaluations(colsNeighbor, zeroRange);
    std::vector< JacobianRangeType > localTestBaseFunctionSetNeighborGradients(colsNeighbor, zeroJacobainRange);
    localTestBaseFunctionSetNeighbor.evaluate(localPointNeighbor, localTestBaseFunctionSetNeighborEvaluations);
    localTestBaseFunctionSetNeighbor.jacobian(localPointNeighbor, localTestBaseFunctionSetNeighborGradients);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_->evaluate(globalPoint, functionValue);

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = penaltyFactor_ / std::pow(intersection.geometry().volume(), 1.0);

    // entity entity combinations
    assert(entityEntityRet.rows() == rowsEntity);
    assert(entityEntityRet.cols() == colsEntity);
    // loop over all ansatz functions
    for(unsigned int i = 0; i < rowsEntity; ++i ) {
      // get row of ret matrix
      typename Dune::DenseMatrix< LocalMatrixImp >::row_reference entityEntityRetRow = entityEntityRet[i];
      // loop over all test function
      for( unsigned int j = 0; j < colsEntity; ++j ) { {
          gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityEntityRetRow[j] += -0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        } {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evalautionTimesNormalTimesGradient = localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          entityEntityRetRow[j] += -0.5 * functionValue * evalautionTimesNormalTimesGradient;
        } {
          const RangeFieldType evalautionTimesEvaluation = localAnsatzBaseFunctionSetEntityEvaluations[i] * localTestBaseFunctionSetEntityEvaluations[j];
          entityEntityRetRow[j] += penaltyParameter * evalautionTimesEvaluation;
      } } // loop over all test function
    } // loop over all ansatz functions

    // do entity neighbor combinations
    assert( entityNeighborRet.rows() == rowsEntity );
    assert( entityNeighborRet.cols() == colsNeighbor );
    for( unsigned int i = 0; i < rowsEntity; ++i )
    {
      for( unsigned int j = 0; j < colsNeighbor; ++j )
      {
        {
          gradientTimesNormal = localTestBaseFunctionSetNeighborGradients[j][0] * unitOuterNormal;
          const RangeType evaluationTimesGradientTimesNormal = localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityNeighborRet[i][j] += -0.5 * functionValue * evaluationTimesGradientTimesNormal;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetNeighborEvaluations[j] * normalTimesGradient;
          entityNeighborRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation = localTestBaseFunctionSetNeighborEvaluations[j] * localAnsatzBaseFunctionSetEntityEvaluations[i];
          entityNeighborRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done entity Neighbor combinations

    // do Neighbor entity combinations
    assert( NeighborEntityRet.rows() == rowsNeighbor );
    assert( NeighborEntityRet.cols() == colsEntity );
    for( unsigned int i = 0; i < rowsNeighbor; ++i )
    {
      for( unsigned int j = 0; j < colsEntity; ++j )
      {
        {
          gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetNeighborEvaluations[i] * gradientTimesNormal;
          NeighborEntityRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetNeighborGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          NeighborEntityRet[i][j] += -0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvalaution = localTestBaseFunctionSetEntityEvaluations[j] * localAnsatzBaseFunctionSetNeighborEvaluations[i];
          NeighborEntityRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvalaution;
        }
      }
    } // done Neighbor entity combinations

    // do Neighbor Neighbor combinations
    assert( NeighborNeighborRet.rows() == rowsNeighbor );
    assert( NeighborNeighborRet.cols() == colsNeighbor );
    for( unsigned int i = 0; i < rowsNeighbor; ++i )
    {
      for( unsigned int j = 0; j < colsNeighbor; ++j )
      {
        {
          gradientTimesNormal = localTestBaseFunctionSetNeighborGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvalaution = localAnsatzBaseFunctionSetNeighborEvaluations[i] * gradientTimesNormal;
          NeighborNeighborRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvalaution;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetNeighborGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient = localTestBaseFunctionSetNeighborEvaluations[j] * normalTimesGradient;
          NeighborNeighborRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation = localTestBaseFunctionSetNeighborEvaluations[j] * localAnsatzBaseFunctionSetNeighborEvaluations[i];
          NeighborNeighborRet[i][j] += penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done Neighbor Neighbor combinations

  } // end method evaluateLocal

private:
  Inner(const ThisType&);
  ThisType& operator=(const ThisType&);

  const Dune::shared_ptr< const InducingFunctionType > inducingFunction_;
  const unsigned int order_;
  const RangeFieldType penaltyFactor_;
}; // end class Inner

#if 0
/**
  \todo       Implement penalty parameter
  \todo       Implement different constructor, for function and discretefunction
  **/
template< class FunctionSpaceImp >
class Dirichlet
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef Dirichlet< FunctionSpaceType >
    ThisType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef Dune::HelperTools::Function::Runtime< FunctionSpaceType >
    InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Dirichlet( const std::string expression = "[1.0;1.0;1.0]", const int order = 1 )
    : inducingFunction_( expression ),
      order_( std::max( 0, order ) )
  {
  }

  //! copy constructor
  Dirichlet( const ThisType& other )
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
    \attention  Assumes, that the return vectors are empty, because we do multiple +=
    **/
  template< class LocalAnsatzBaseFunctionSetType,
            class LocalTestBaseFunctionSetType,
            class IntersectionType,
            class LocalPointType,
            class LocalMatrixType >
  void evaluateLocal( const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                      const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                      const IntersectionType& intersection,
                      const LocalPointType& localPoint,
                      LocalMatrixType& ret ) const
  {
    // some stuff
    const DomainType globalPoint = intersection.geometry().global( localPoint );
    const DomainType localPointEntity = intersection.geometryInInside().global( localPoint );
    const DomainType unitOuterNormal = intersection.unitOuterNormal( localPoint );

    // evaluate ansatz basefunctionset
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetEvaluations( rows, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetGradients( rows, JacobianRangeType( 0.0 ) );
    localAnsatzBaseFunctionSet.evaluate( localPointEntity, localAnsatzBaseFunctionSetEvaluations );
    localAnsatzBaseFunctionSet.jacobian( localPointEntity, localAnsatzBaseFunctionSetGradients );

    // evaluate test basefunctionset
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector< RangeType > localTestBaseFunctionSetEvaluations( cols, RangeType( 0.0 ) );
    std::vector< JacobianRangeType > localTestBaseFunctionSetGradients( cols, JacobianRangeType( 0.0 ) );
    localTestBaseFunctionSet.evaluate( localPointEntity, localTestBaseFunctionSetEvaluations );
    localTestBaseFunctionSet.jacobian( localPointEntity, localTestBaseFunctionSetGradients );

    // evaluate inducing function
    RangeType functionValue( 0.0 );
    inducingFunction_.evaluate( globalPoint, functionValue );

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = 20.0 / std::pow( intersection.geometry().volume(), 1.0 );

    // do loop over all ansatz and test basefunctions
    assert( ret.rows() == rows );
    assert( ret.cols() == cols );
    for( unsigned int i = 0; i < rows; ++i )
    {
      for( unsigned int j = 0; j < cols; ++j )
      {
        {
          gradientTimesNormal = localTestBaseFunctionSetGradients[j][0] * unitOuterNormal;
          gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetEvaluations[i] * gradientTimesNormal;
          ret[i][j] += -1.0 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetGradients[i][0];
          const RangeType evalautionTimesNormalTimesGradient = localTestBaseFunctionSetEvaluations[j] * normalTimesGradient;
          ret[i][j] += -1.0 * functionValue * evalautionTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evalautionTimesEvaluation = localAnsatzBaseFunctionSetEvaluations[i] * localTestBaseFunctionSetEvaluations[j];
          ret[i][j] += penaltyParameter * evalautionTimesEvaluation;
        }
      }
    } // done loop over all ansatz and test basefunctions

  } // end method evaluateLocal

private:
  //! assignment operator
  ThisType& operator=( const ThisType& );

  const InducingFunctionType inducingFunction_;
  const unsigned int order_;
}; // end class Dirichlet
#endif

} // end namespace IPDGfluxes

} // end namespace Quaternary

} // end namespace Local

} // end namespace Evaluation

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
