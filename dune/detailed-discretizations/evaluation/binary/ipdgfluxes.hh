#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_IPDGFLUXES_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_IPDGFLUXES_HH

namespace Dune {

namespace DetailedDiscretizations {

namespace Evaluation {

namespace Binary {

namespace IPDGFluxes {

/**
  \brief  This represents the operation \f$a\nabla u \nabla v\f$.

          \f$a\f$ is a given scalar function (in this case 1) and \f$u\f$ and \f$u\f$ may be local functions, i.e.
          ansatz- and testfunctions.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$, \f$u\f$ and \f$v\f$ live in.
  **/
template <class FunctionSpaceImp>
class JumpMeanPenalty
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef JumpMeanPenalty<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  JumpMeanPenalty(const std::string expression = "[1.0;1.0;1.0]", const int order = 1)
    : inducingFunction_(expression)
    , order_(std::max(0, order))
  {
  }

  //! copy constructor
  JumpMeanPenalty(const ThisType& other)
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

  template <class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType, class LocalMatrixType>
  void evaluate(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const DomainType& localPointEntity,
                const DomainType& localPointNeighbour, const DomainType& unitOuterNormal, LocalMatrixType& ret) const
  {
    // get global point
    const DomainType globalPoint = localAnsatzBaseFunctionSet.entity().geometry().global(localPointEntity);

    // evaluate first gradient
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalAnsatzBaseFunctionSet(rows, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSet.jacobian(localPointEntity, gradientLocalAnsatzBaseFunctionSet);

    // evaluate second gradient
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalTestBaseFunctionSet(cols, JacobianRangeType(0.0));
    localTestBaseFunctionSet.jacobian(localPointNeighbour, gradientLocalTestBaseFunctionSet);

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
}; // end class JumpMeanPenalty

} // end namespace IPDGFluxes

} // end namespace Binary

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_BINARY_IPDGFLUXES_HH
