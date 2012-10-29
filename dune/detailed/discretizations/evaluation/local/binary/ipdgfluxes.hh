#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_IPDGFLUXES_HH

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/densematrix.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>
//#include <dune/stuff/common/print.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Evaluation {

namespace Local {

namespace Binary {

namespace IPDGfluxes {

/**
  \todo       Implement penalty parameter
  \todo       Implement different constructor, for function and discretefunction
  **/
template< class FunctionSpaceImp,
          class InducingFunctionImp = Dune::Stuff::Function::Expression< typename FunctionSpaceImp::DomainFieldType, FunctionSpaceImp::DimDomain, typename FunctionSpaceImp::RangeFieldType, FunctionSpaceImp::DimRange > >
class Dirichlet
{
public:

  typedef FunctionSpaceImp FunctionSpaceType;

  typedef InducingFunctionImp InducingFunctionType;

  typedef Dirichlet< FunctionSpaceType, InducingFunctionType > ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  Dirichlet(const Dune::shared_ptr< const InducingFunctionType > inducingFunction,
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
    \attention  Assumes, that the return vectors are empty, because we do multiple +=
    **/
  template< class LocalAnsatzBaseFunctionSetType,
            class LocalTestBaseFunctionSetType,
            class IntersectionType,
            class LocalPointType,
            class LocalMatrixImp >
  void evaluateLocal(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                     const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                     const IntersectionType& intersection,
                     const LocalPointType& localPoint,
                     Dune::DenseMatrix< LocalMatrixImp >& ret) const
  {
    // some stuff
    const DomainType globalPoint = intersection.geometry().global(localPoint);
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    const RangeType zeroRange(0);
    const JacobianRangeType zeroJacobainRange(0);

//    // info
//    Stuff::Common::print(globalPoint, "globalPoint");
//    Stuff::Common::print(localPointEntity, "localPointEntity");
//    Stuff::Common::print(unitOuterNormal, "unitOuterNormal");
//    Stuff::Common::print(ret, "ret");

    // evaluate ansatz basefunctionset
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector< RangeType > localAnsatzBaseFunctionSetEvaluations(rows, zeroRange);
    std::vector< JacobianRangeType > localAnsatzBaseFunctionSetGradients(rows, zeroJacobainRange);
    localAnsatzBaseFunctionSet.evaluate(localPointEntity, localAnsatzBaseFunctionSetEvaluations);
//    Stuff::Common::print(localAnsatzBaseFunctionSetEvaluations, "localAnsatzBaseFunctionSetEvaluations");
    localAnsatzBaseFunctionSet.jacobian(localPointEntity, localAnsatzBaseFunctionSetGradients);
//    Stuff::Common::print(localAnsatzBaseFunctionSetGradients, "localAnsatzBaseFunctionSetGradients");

    // evaluate test basefunctionset
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector< RangeType > localTestBaseFunctionSetEvaluations(cols, zeroRange);
    std::vector< JacobianRangeType > localTestBaseFunctionSetGradients(cols, zeroJacobainRange);
    localTestBaseFunctionSet.evaluate(localPointEntity, localTestBaseFunctionSetEvaluations);
//    Stuff::Common::print(localTestBaseFunctionSetEvaluations, "localTestBaseFunctionSetEvaluations");
    localTestBaseFunctionSet.jacobian(localPointEntity, localTestBaseFunctionSetGradients);
//    Stuff::Common::print(localTestBaseFunctionSetGradients, "localTestBaseFunctionSetGradients");

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_->evaluate(globalPoint, functionValue);
//    Stuff::Common::print(functionValue, "functionValue");

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = penaltyFactor_ / std::pow(intersection.geometry().volume(), 1.0);
//    Stuff::Common::print(penaltyParameter, "penaltyParameter");

    // do loop over all ansatz and test basefunctions
    assert(ret.rows() == rows);
    assert(ret.cols() == cols);
    // loop over all ansatz functions
    for (unsigned int i = 0; i < rows; ++i) {
      // get row of ret matrix
      typename Dune::DenseMatrix< LocalMatrixImp >::row_reference retRow = ret[i];
      // loop over all test function
      for (unsigned int j = 0; j < cols; ++j) { {
//          std::cout << "i = " << i << ", j = " << j << std::endl;
//          Stuff::Common::print(retRow[j], "retRow[j]", std::cout, "  ");
          const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetGradients[j][0] * unitOuterNormal;
//          Stuff::Common::print(gradientTimesNormal, "gradientTimesNormal", std::cout, "  ");
          const RangeType gradientTimesNormalTimesEvaluation = localAnsatzBaseFunctionSetEvaluations[i] * gradientTimesNormal;
//          Stuff::Common::print(gradientTimesNormalTimesEvaluation, "gradientTimesNormalTimesEvaluation", std::cout, "  ");
//          Stuff::Common::print(-1.0 * functionValue * gradientTimesNormalTimesEvaluation, "-1.0 * functionValue * gradientTimesNormalTimesEvaluation", std::cout, "  ");
          retRow[j] += -1.0 * functionValue * gradientTimesNormalTimesEvaluation;
//          Stuff::Common::print(retRow[j], "retRow[j]", std::cout, "  ");
        } {
          const RangeFieldType normalTimesGradient = unitOuterNormal * localAnsatzBaseFunctionSetGradients[i][0];
//          Stuff::Common::print(normalTimesGradient, "normalTimesGradient", std::cout, "  ");
          const RangeType evalautionTimesNormalTimesGradient = localTestBaseFunctionSetEvaluations[j] * normalTimesGradient;
//          Stuff::Common::print(evalautionTimesNormalTimesGradient, "evalautionTimesNormalTimesGradient", std::cout, "  ");
//          Stuff::Common::print(-1.0 * functionValue * evalautionTimesNormalTimesGradient, "-1.0 * functionValue * evalautionTimesNormalTimesGradient", std::cout, "  ");
          retRow[j] += -1.0 * functionValue * evalautionTimesNormalTimesGradient;
//          Stuff::Common::print(retRow[j], "retRow[j]", std::cout, "  ");
        } {
          const RangeFieldType evalautionTimesEvaluation = localAnsatzBaseFunctionSetEvaluations[i] * localTestBaseFunctionSetEvaluations[j];
//          Stuff::Common::print(evalautionTimesEvaluation, "evalautionTimesEvaluation", std::cout, "  ");
//          Stuff::Common::print(penaltyParameter * evalautionTimesEvaluation, "penaltyParameter * evalautionTimesEvaluation", std::cout, "  ");
          retRow[j] += penaltyParameter * evalautionTimesEvaluation;
//          Stuff::Common::print(retRow[j], "retRow[j]", std::cout, "  ");
      } } // loop over all test function
    } // loop over all ansatz functions
//    Stuff::Common::print(ret, "ret");
  } // void evaluateLocal(...) const

private:
  Dirichlet(const ThisType&);
  ThisType& operator=(const ThisType&);

  const Dune::shared_ptr< const InducingFunctionType > inducingFunction_;
  const unsigned int order_;
  const RangeFieldType penaltyFactor_;
}; // end class Dirichlet

} // namespace IPDGfluxes

} // namespace Binary

} // namespace Local

} // namespace Evaluation

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_IPDGFLUXES_HH
