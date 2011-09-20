#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH

// dune-helper-tools includes
#include <dune/helper-tools/function/runtime.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace Evaluation {

namespace Local {

namespace Quaternary {

namespace IPDGFluxes {

/**
  \todo       Implement penalty parameter
  \todo       Implement different constructor, for function and discretefunction
  **/
template <class FunctionSpaceImp>
class Inner
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef Inner<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Inner(const std::string expression = "[1.0;1.0;1.0]", const int order = 1)
    : inducingFunction_(expression)
    , order_(std::max(0, order))
  {
  }

  //! copy constructor
  Inner(const ThisType& other)
    : inducingFunction_(other.inducingFunction())
    , order_(other.order())
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
    \attention  Assumes, that the return matrices are empty, because we do multiple +=
    \todo       There are currently four nested loops, these can be simplyfied to two (see msdg-matlab)
    \todo       Rename Entity -> En, Neighbour -> Ne (otherwise stuff is too long)
    **/
  template <class LocalAnsatzBaseFunctionSetEntityType, class LocalAnsatzBaseFunctionSetNeighbourType,
            class LocalTestBaseFunctionSetEntityType, class LocalTestBaseFunctionSetNeighbourType,
            class IntersectionType, class LocalMatrixType>
  void evaluateLocal(const LocalAnsatzBaseFunctionSetEntityType& localAnsatzBaseFunctionSetEntity,
                     const LocalAnsatzBaseFunctionSetNeighbourType& localAnsatzBaseFunctionSetNeighbour,
                     const LocalTestBaseFunctionSetEntityType& localTestBaseFunctionSetEntity,
                     const LocalTestBaseFunctionSetNeighbourType& localTestBaseFunctionSetNeighbour,
                     const IntersectionType& intersection, const DomainType& localPoint,
                     LocalMatrixType& entityEntityRet, LocalMatrixType& entityNeighbourRet,
                     LocalMatrixType& neighbourEntityRet, LocalMatrixType& neighbourNeighbourRet) const
  {
    // some stuff
    const DomainType globalPoint     = intersection.geometry().global(localPoint);
    const DomainType localPointEn    = intersection.geometryInInside().global(localPoint);
    const DomainType localPointNe    = intersection.geometryInOutside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);

    // evaluate ansatz entity basefunctionset
    const unsigned int rowsEntity = localAnsatzBaseFunctionSetEntity.size();
    std::vector<RangeType> localAnsatzBaseFunctionSetEntityEvaluations(rowsEntity, RangeType(0.0));
    std::vector<JacobianRangeType> localAnsatzBaseFunctionSetEntityGradients(rowsEntity, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSetEntity.evaluate(localPointEn, localAnsatzBaseFunctionSetEntityEvaluations);
    localAnsatzBaseFunctionSetEntity.jacobian(localPointEn, localAnsatzBaseFunctionSetEntityGradients);

    // evaluate ansatz neighbour basefunctionset
    const unsigned int rowsNeighbour = localAnsatzBaseFunctionSetNeighbour.size();
    std::vector<RangeType> localAnsatzBaseFunctionSetNeighbourEvaluations(rowsNeighbour, RangeType(0.0));
    std::vector<JacobianRangeType> localAnsatzBaseFunctionSetNeighbourGradients(rowsNeighbour, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSetNeighbour.evaluate(localPointNe, localAnsatzBaseFunctionSetNeighbourEvaluations);
    localAnsatzBaseFunctionSetNeighbour.jacobian(localPointNe, localAnsatzBaseFunctionSetNeighbourGradients);

    // evaluate test entity basefunctionset
    const unsigned int colsEntity = localTestBaseFunctionSetEntity.size();
    std::vector<RangeType> localTestBaseFunctionSetEntityEvaluations(colsEntity, RangeType(0.0));
    std::vector<JacobianRangeType> localTestBaseFunctionSetEntityGradients(colsEntity, JacobianRangeType(0.0));
    localTestBaseFunctionSetEntity.evaluate(localPointEn, localTestBaseFunctionSetEntityEvaluations);
    localTestBaseFunctionSetEntity.jacobian(localPointEn, localTestBaseFunctionSetEntityGradients);

    // evaluate test neighbour basefunctionset
    const unsigned int colsNeighbour = localTestBaseFunctionSetNeighbour.size();
    std::vector<RangeType> localTestBaseFunctionSetNeighbourEvaluations(colsNeighbour, RangeType(0.0));
    std::vector<JacobianRangeType> localTestBaseFunctionSetNeighbourGradients(colsNeighbour, JacobianRangeType(0.0));
    localTestBaseFunctionSetNeighbour.evaluate(localPointNe, localTestBaseFunctionSetNeighbourEvaluations);
    localTestBaseFunctionSetNeighbour.jacobian(localPointNe, localTestBaseFunctionSetNeighbourGradients);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_.evaluate(globalPoint, functionValue);

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = 1.0 / std::pow(intersection.geometry().volume(), 1.0);

    // entity entity combinations
    // do loop over all ansatz and test basefunctions
    assert(entityEntityRet.rows() == rowsEntity);
    assert(entityEntityRet.cols() == colsEntity);
    for (unsigned int i = 0; i < rowsEntity; ++i) {
      for (unsigned int j = 0; j < colsEntity; ++j) {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvaluation =
              localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityEntityRet[i][j] += -0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evalautionTimesNormalTimesGradient =
              localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          entityEntityRet[i][j] += -0.5 * functionValue * evalautionTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evalautionTimesEvaluation =
              localAnsatzBaseFunctionSetEntityEvaluations[i] * localTestBaseFunctionSetEntityEvaluations[j];
          entityEntityRet[i][j] += penaltyParameter * evalautionTimesEvaluation;
        }
      }
    } // done entity entity combinations

    // do entity neighbour combinations
    assert(entityNeighbourRet.rows() == rowsEntity);
    assert(entityNeighbourRet.cols() == colsNeighbour);
    for (unsigned int i = 0; i < rowsEntity; ++i) {
      for (unsigned int j = 0; j < colsNeighbour; ++j) {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetNeighbourGradients[j][0] * unitOuterNormal;
          const RangeType evaluationTimesGradientTimesNormal =
              localAnsatzBaseFunctionSetEntityEvaluations[i] * gradientTimesNormal;
          entityNeighbourRet[i][j] += -0.5 * functionValue * evaluationTimesGradientTimesNormal;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetEntityGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient =
              localTestBaseFunctionSetNeighbourEvaluations[j] * normalTimesGradient;
          entityNeighbourRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation =
              localTestBaseFunctionSetNeighbourEvaluations[j] * localAnsatzBaseFunctionSetEntityEvaluations[i];
          entityNeighbourRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done entity neighbour combinations

    // do neighbour entity combinations
    assert(neighbourEntityRet.rows() == rowsNeighbour);
    assert(neighbourEntityRet.cols() == colsEntity);
    for (unsigned int i = 0; i < rowsNeighbour; ++i) {
      for (unsigned int j = 0; j < colsEntity; ++j) {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetEntityGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvaluation =
              localAnsatzBaseFunctionSetNeighbourEvaluations[i] * gradientTimesNormal;
          neighbourEntityRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient =
              unitOuterNormal * localAnsatzBaseFunctionSetNeighbourGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient =
              localTestBaseFunctionSetEntityEvaluations[j] * normalTimesGradient;
          neighbourEntityRet[i][j] += -0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvalaution =
              localTestBaseFunctionSetEntityEvaluations[j] * localAnsatzBaseFunctionSetNeighbourEvaluations[i];
          neighbourEntityRet[i][j] += -1.0 * penaltyParameter * evaluationTimesEvalaution;
        }
      }
    } // done neighbour entity combinations

    // do neighbour neighbour combinations
    assert(neighbourNeighbourRet.rows() == rowsNeighbour);
    assert(neighbourNeighbourRet.cols() == colsNeighbour);
    for (unsigned int i = 0; i < rowsNeighbour; ++i) {
      for (unsigned int j = 0; j < colsNeighbour; ++j) {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetNeighbourGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvalaution =
              localAnsatzBaseFunctionSetNeighbourEvaluations[i] * gradientTimesNormal;
          neighbourNeighbourRet[i][j] += 0.5 * functionValue * gradientTimesNormalTimesEvalaution;
        }
        {
          const RangeFieldType normalTimesGradient =
              unitOuterNormal * localAnsatzBaseFunctionSetNeighbourGradients[i][0];
          const RangeType evaluationTimesNormalTimesGradient =
              localTestBaseFunctionSetNeighbourEvaluations[j] * normalTimesGradient;
          neighbourNeighbourRet[i][j] += 0.5 * functionValue * evaluationTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evaluationTimesEvaluation =
              localTestBaseFunctionSetNeighbourEvaluations[j] * localAnsatzBaseFunctionSetNeighbourEvaluations[i];
          neighbourNeighbourRet[i][j] += penaltyParameter * evaluationTimesEvaluation;
        }
      }
    } // done neighbour neighbour combinations

  } // end method evaluateLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const InducingFunctionType inducingFunction_;
  unsigned int order_;
}; // end class Inner

/**
  \todo       Implement penalty parameter
  \todo       Implement different constructor, for function and discretefunction
  **/
template <class FunctionSpaceImp>
class Dirichlet
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef Dirichlet<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Dirichlet(const std::string expression = "[1.0;1.0;1.0]", const int order = 1)
    : inducingFunction_(expression)
    , order_(std::max(0, order))
  {
  }

  //! copy constructor
  Dirichlet(const ThisType& other)
    : inducingFunction_(other.inducingFunction())
    , order_(other.order())
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
  template <class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType, class IntersectionType,
            class LocalMatrixType>
  void evaluateLocal(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                     const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const IntersectionType& intersection,
                     const DomainType& localPoint, LocalMatrixType& ret) const
  {
    // some stuff
    const DomainType globalPoint     = intersection.geometry().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);

    // evaluate ansatz basefunctionset
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector<RangeType> localAnsatzBaseFunctionSetEvaluations(rows, RangeType(0.0));
    std::vector<JacobianRangeType> localAnsatzBaseFunctionSetGradients(rows, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSet.evaluate(localPoint, localAnsatzBaseFunctionSetEvaluations);
    localAnsatzBaseFunctionSet.jacobian(localPoint, localAnsatzBaseFunctionSetGradients);

    // evaluate test basefunctionset
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector<RangeType> localTestBaseFunctionSetEvaluations(cols, RangeType(0.0));
    std::vector<JacobianRangeType> localTestBaseFunctionSetGradients(cols, JacobianRangeType(0.0));
    localTestBaseFunctionSet.evaluate(localPoint, localTestBaseFunctionSetEvaluations);
    localTestBaseFunctionSet.jacobian(localPoint, localTestBaseFunctionSetGradients);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_.evaluate(globalPoint, functionValue);

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = 1.0 / std::pow(intersection.geometry().volume(), 1.0);

    // do loop over all ansatz and test basefunctions
    assert(ret.rows() == rows);
    assert(ret.cols() == cols);
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
        {
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetGradients[j][0] * unitOuterNormal;
          const RangeType gradientTimesNormalTimesEvaluation =
              localAnsatzBaseFunctionSetEvaluations[i] * gradientTimesNormal;
          ret[i][j] += -1.0 * functionValue * gradientTimesNormalTimesEvaluation;
        }
        {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetGradients[i][0];
          const RangeType evalautionTimesNormalTimesGradient =
              localTestBaseFunctionSetEvaluations[j] * normalTimesGradient;
          ret[i][j] += -1.0 * functionValue * evalautionTimesNormalTimesGradient;
        }
        {
          const RangeFieldType evalautionTimesEvaluation =
              localAnsatzBaseFunctionSetEvaluations[i] * localTestBaseFunctionSetEvaluations[j];
          ret[i][j] += penaltyParameter * evalautionTimesEvaluation;
        }
      }
    } // done loop over all ansatz and test basefunctions

  } // end method evaluateLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const InducingFunctionType inducingFunction_;
  unsigned int order_;
}; // end class Dirichlet

} // end namespace IPDGFluxes

} // end namespace Quaternary

} // end namespace Local

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_QUATERNARY_IPDGFLUXES_HH
