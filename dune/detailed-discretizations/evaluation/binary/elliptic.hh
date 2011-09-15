#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_ELLIPTIC_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_ELLIPTIC_HH

// dune-helper-tools includes
#include <dune/helper-tools/function/runtime.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace Evaluation {

namespace Binary {

/**
  \brief  This represents the operation \f$a\nabla u \nabla v\f$.

          \f$a\f$ is a given scalar function (in this case 1) and \f$u\f$ and \f$u\f$ may be local functions, i.e.
          ansatz- and testfunctions.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$, \f$u\f$ and \f$v\f$ live in.
  **/
template <class FunctionSpaceImp>
class Elliptic
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef Elliptic<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Elliptic(const std::string expression = "[1.0;1.0;1.0]", const int order = 1)
    : inducingFunction_(expression)
    , order_(std::max(0, order))
  {
  }

  //! copy constructor
  Elliptic(const Elliptic& other)
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
    * \brief      Evaluates \f$a(x)\nabla u(x) \nabla v(x)\f$ for a given local point \f$x\f$.
    *
    * \tparam     LocalAnsatzFunctionType
    *             Type of the local ansatz function \f$u\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalTestFunctionType
    *             Type of the local test function \f$v\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalPointType
    *             Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    * \param[in]  localAnsatzFunction
    *             The local function \f$u\f$.
    * \param[in]  localTestFunction
    *             The local function \f$v\f$.
    * \param[in]  localPoint
    *             The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
    *             element.
    * \return     \f$a(x)\nabla u(x) \nabla v(x)\f$
    **/
  template <class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType, class LocalMatrixType>
  void evaluate(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const DomainType& localPoint,
                LocalMatrixType& ret) const
  {
    // get global point
    const DomainType globalPoint = localAnsatzBaseFunctionSet.entity().geometry().global(localPoint);

    // evaluate first gradient
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalAnsatzBaseFunctionSet(rows, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSet.jacobian(localPoint, gradientLocalAnsatzBaseFunctionSet);

    // evaluate second gradient
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalTestBaseFunctionSet(cols, JacobianRangeType(0.0));
    localTestBaseFunctionSet.jacobian(localPoint, gradientLocalTestBaseFunctionSet);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_.evaluate(globalPoint, functionValue);

    // do loop over all ansatz and test basefunctions
    assert(ret.rows() == rows);
    assert(ret.cols() == cols);
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
        const RangeFieldType gradientProduct =
            gradientLocalAnsatzBaseFunctionSet[i][0] * gradientLocalTestBaseFunctionSet[j][0];
        ret[i][j] = functionValue * gradientProduct;
      }
    }
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const InducingFunctionType inducingFunction_;
  unsigned int order_;
}; // end class Elliptic

} // end namespace Binary

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_ELLIPTIC_HH
