#ifndef DUNE_FUNCTIONALS_LOCALEVALUATION_WRAPPER_HH
#define DUNE_FUNCTIONALS_LOCALEVALUATION_WRAPPER_HH

namespace Dune {

namespace Functionals {

namespace LocalEvaluation {

template <class BaseLocalEvaluationImp, class DiscreteFunctionImp>
class Wrapper
{
public:
  typedef BaseLocalEvaluationImp BaseLocalEvaluationType;

  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename BaseLocalEvaluationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  Wrapper(const BaseLocalEvaluationType& baseLocalEvaluation, const DiscreteFunctionType& discreteFunction)
    : baseLocalEvaluation_(baseLocalEvaluation)
    , discreteFunction_(discreteFunction)
  {
  }

  template <class LocalTestFunctionType>
  const RangeFieldType evaluate(const LocalTestFunctionType& localTestFunction, const DomainType& localPoint) const
  {
    // get local inducing function
    typedef typename DiscreteFunctionType::LocalFunctionType InducingLocalFunctionType;

    this does not work :

        const InducingLocalFunctionType inducingLocalFunction =
            discreteFunction_.localFunction(localTestFunction.entity());

    // evaluate base local evalaution
    return baseLocalEvaluation_.evaluate(inducingLocalFunction, localTestFunction, localPoint);
  }

private:
  const BaseLocalEvaluationType& baseLocalEvaluation_;
  const DiscreteFunctionType& discreteFunction_;

}; // end class Wrapper

} // end namespace LocalEvaluation

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FUNCTIONALS_LOCALEVALUATION_WRAPPER_HH
