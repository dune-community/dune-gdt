#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_IPDGFLUXES_HH

// dune-helper-tools includes
#include <dune/helper-tools/function/runtime.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace Evaluation {

namespace Local {

namespace Unary {

/**
  \todo   Should be parameterized with the inducing and the dirichlet boundary function.
  \todo   Add neumann.
  \todo   Add penalty parameter.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$ and \f$v\f$ live in.
  **/
template <class FunctionSpaceImp>
class IPDGFluxes
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef IPDGFluxes<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> DirichletFunctionType;

  //! constructor, takes the inducing and the dirichlet functions expressions as runtime parameters
  IPDGFluxes(const std::string inducingExpression = "[1.0;1.0;1.0]",
             const std::string dirichletExpression = "[0.0;0.0;0.0]", const int order = 1)
    : inducingFunction_(inducingExpression)
    , dirichletFunction_(dirichletExpression)
    , order_(std::max(0, order))
  {
  }

  //! copy constructor
  IPDGFluxes(const ThisType& other)
    : inducingFunction_(other.inducingFunction())
    , dirichletFunction_(other.dirichletFunction())
    , order_(other.order())
  {
  }

  //! returns the inducing function
  InducingFunctionType inducingFunction() const
  {
    return inducingFunction_;
  }

  DirichletFunctionType dirichletFunction() const
  {
    return dirichletFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    \attention  Assumes ret to be cleared, since we do multiple +=
    **/
  template <class LocalTestBaseFunctionSetType, class IntersectionType, class LocalVectorType>
  void evaluateLocal(const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const IntersectionType& intersection,
                     const DomainType& localPoint, LocalVectorType& ret) const
  {
    // some stuff
    const DomainType globalPoint     = intersection.geometry().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal();

    // evaluate inducing and dirichlet function
    RangeType inducingFunctionEvaluation(0.0);
    inducingFunction_.evaluate(globalPoint, inducingFunctionEvaluation);
    RangeType dirichletFunctionEvaluation(0.0);
    inducingFunction_.evaluate(globalPoint, dirichletFunctionEvaluation);

    // evaluate test basefunctionset
    const unsigned int size = localTestBaseFunctionSet.size();
    std::vector<RangeType> localTestBaseFunctionSetEvaluations(size, RangeType(0.0));
    localTestBaseFunctionSet.evaluate(localPoint, localTestBaseFunctionSetEvaluations);
    std::vector<JacobianRangeType> localTestBaseFunctionSetGradients(size, JacobianRangeType(0.0));
    localTestBaseFunctionSet.jacobian(localPoint, localTestBaseFunctionSetGradients);

    // evaluate penalty parameter
    const RangeFieldType penaltyParameter = 1.0 / std::pow(intersection.geometry().volume(), 1.0);

    // do loop over all basis functions
    assert(ret.size() >= size);
    for (unsigned int i = 0; i < size; ++i) {
      {
        const RangeType evaluationTimesDirichlet = localTestBaseFunctionSetEvaluations[i] * dirichletFunctionEvaluation;
        ret[i]                                   = penaltyParameter * evaluationTimesDirichlet;
      }
      {
        const RangeFieldType gradientTimesNormal = localTestBaseFunctionSetGradients[i][0] * unitOuterNormal;
        ret[i] += -1.0 * inducingFunctionEvaluation * gradientTimesNormal * dirichletFunctionEvaluation;
      }
    }
  } // end method evaluateLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const InducingFunctionType inducingFunction_;
  const DirichletFunctionType dirichletFunction_;
  const unsigned int order_;
}; // end class IPDGFluxes

} // end namespace Unary

} // end namespace Local

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_IPDGFLUXES_HH
