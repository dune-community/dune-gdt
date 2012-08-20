#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_SCALE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_SCALE_HH

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Evaluation {

namespace Local {

namespace Unary {

/**
  \brief  This represents the operation \f$fv\f$.

          \f$f\f$ is a given right hand side (in this case 1) and \f$v\f$ may be a local function, i.e. a
          testfunction.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$ and \f$v\f$ live in.
  **/
template <class FunctionSpaceImp,
          class InducingFunctionImp =
              Dune::Stuff::Function::Expression<typename FunctionSpaceImp::DomainFieldType, FunctionSpaceImp::DimDomain,
                                                typename FunctionSpaceImp::RangeFieldType, FunctionSpaceImp::DimRange>>
class Scale
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef Scale<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef InducingFunctionImp InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Scale(const Dune::shared_ptr<const InducingFunctionType> inducingFunction, const unsigned int order)
    : inducingFunction_(inducingFunction)
    , order_(order)
  {
  }

  //  Scale(const Dune::ParameterTree& paramTree)
  //    : inducingFunction_(new InducingFunctionType(paramTree)),
  //      order_(paramTree.get("order", 0))
  //  {
  //  }

  //! returns the inducing function
  Dune::shared_ptr<const InducingFunctionType> inducingFunction() const
  {
    return inducingFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    \brief      Evaluates \f$f(x)v(x)\f$ for a given local point \f$x\f$.

    \tparam     LocalTestFunctionType
                Type of the local function \f$v\f$, i.e. Dune::LocalFunction.
    \tparam     LocalPointType
                Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    \param[in]  localTestFunction
                The local function \f$v\f$.
    \param[in]  localPoint
                The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
                element.
    \return     \f$f(x)v(x)\f$
    **/
  template <class LocalTestBaseFunctionSetType, class LocalVectorType>
  void evaluateLocal(const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const DomainType& localPoint,
                     LocalVectorType& ret) const
  {
    // get global point
    const DomainType globalPoint = localTestBaseFunctionSet.entity().geometry().global(localPoint);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_->evaluate(globalPoint, functionValue);

    // evaluate set of local basis functions
    const unsigned int size = localTestBaseFunctionSet.size();
    std::vector<RangeType> valuesLocalBaseFunctionSet(size, RangeType(0.0));
    localTestBaseFunctionSet.evaluate(localPoint, valuesLocalBaseFunctionSet);

    // do loop over all basis functions
    assert(ret.size() == size);
    for (unsigned int i = 0; i < size; ++i) {
      ret[i] = functionValue * valuesLocalBaseFunctionSet[i];
    }
  } // end method evalauteLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  //! copy constructor
  Scale(const ThisType&);

  const Dune::shared_ptr<const InducingFunctionType> inducingFunction_;
  const unsigned int order_;
}; // end class Product

} // end namespace Unary

} // end namespace Local

} // end namespace Evaluation

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_UNARY_SCALE_HH
