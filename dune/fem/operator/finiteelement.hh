#ifndef DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH

// dune fem includes
#include <dune/fem/function/common/function.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/common/localmatrix.hh>
#include <dune/fem/common/localvector.hh>
#include <dune/fem/common/localbasefunction.hh>
#include <dune/fem/functional/finiteelement.hh>

namespace Dune {

//! dune-fem-functionals
namespace Functionals {

//! Contains several operators.
namespace Operator {

template <class OperatorLocalOperationImp, class InducingFunctionImp>
class FunctionalLocalOperation
{
public:
  typedef OperatorLocalOperationImp OperatorLocalOperationType;

  typedef InducingFunctionImp InducingFunctionType;

  FunctionalLocalOperation(const OperatorLocalOperationType& operatorLocalOperation,
                           const InducingFunctionType& inducingFunction)
    : operatorLocalOperation_(operatorLocalOperation)
    , inducingFunction_(inducingFunction)
  {
  }

  template <class LocalFunctionType, class LocalPointType>
  double operate(const LocalFunctionType& localFunction, const LocalPointType& localPoint) const
  {

    // some types we will need
    typedef typename LocalFunctionType::EntityType EntityType;

    typedef typename InducingFunctionType::LocalFunctionType InducingLocalFunctionType;

    // entity
    const EntityType& entity = localFunction.entity();

    // local function of the inducing function
    const InducingLocalFunctionType inducingLocalFunction = inducingFunction_.localFunction(entity);

    // evaluate the original operation
    const double ret = operatorLocalOperation_.operate(inducingLocalFunction, localFunction, localPoint);

    // return
    return ret;
  }

private:
  const OperatorLocalOperationType operatorLocalOperation_;
  const InducingFunctionType& inducingFunction_;


}; // end of class FunctionalLocalOperation


template <class DiscreteAnsatzFunctionSpaceImp, class DiscreteTestFunctionSpaceImp, class LocalOperationImp>
class FiniteElement
{
public:
  typedef DiscreteAnsatzFunctionSpaceImp DiscreteAnsatzFunctionSpaceType;

  typedef DiscreteTestFunctionSpaceImp DiscreteTestFunctionSpaceType;

  typedef LocalOperationImp LocalOperationType;

  typedef typename DiscreteAnsatzFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::Functionals::Common::LocalMatrix<RangeFieldType> LocalMatrixType;

  typedef Dune::Functionals::Common::LocalVector<RangeFieldType> LocalVectorType;

  typedef typename DiscreteAnsatzFunctionSpaceType::BaseFunctionSetType AnsatzBaseFunctionSetType;

  typedef typename DiscreteTestFunctionSpaceType::BaseFunctionSetType TestBaseFunctionSetType;

  typedef typename DiscreteAnsatzFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteAnsatzFunctionSpaceType::IteratorType EntityIteratorType;

  typedef typename EntityIteratorType::Entity EntityType;

  typedef typename EntityType::Geometry EntityGeometryType;

  typedef CachingQuadrature<GridPartType, 0> EntityQuadratureType;

  typedef Dune::Functionals::Common::LocalBaseFunctionProvider<DiscreteAnsatzFunctionSpaceType>
      LocalAnsatzBaseFunctionProviderType;

  typedef Dune::Functionals::Common::LocalBaseFunctionProvider<DiscreteTestFunctionSpaceType>
      LocalTestBaseFunctionProviderType;

  typedef typename LocalAnsatzBaseFunctionProviderType::LocalBaseFunctionType LocalAnsatzBaseFunctionType;

  typedef typename LocalTestBaseFunctionProviderType::LocalBaseFunctionType LocalTestBaseFunctionType;

  FiniteElement(const DiscreteAnsatzFunctionSpaceType& ansatzSpace, const DiscreteTestFunctionSpaceType& testSpace,
                const LocalOperationType& localOperation)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
    , localOperation_(localOperation)
    , localAnsatzBaseFunctionProvider_(ansatzSpace)
    , localTestBaseFunctionProvider_(testSpace)
  {
  }

  ~FiniteElement()
  {
  }

  const DiscreteAnsatzFunctionSpaceType& ansatzSpace() const
  {
    return ansatzSpace_;
  }

  const DiscreteTestFunctionSpaceType& testSpace() const
  {
    return testSpace_;
  }

  LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  LocalMatrixType applyLocal(const EntityType& entity) const
  {
    const unsigned numberOfLocalAnsatzDoFs = ansatzSpace_.baseFunctionSet(entity).numBaseFunctions();
    const unsigned numberOfLocalTestDoFs   = testSpace_.baseFunctionSet(entity).numBaseFunctions();

    // init return matrix
    LocalMatrixType ret(numberOfLocalAnsatzDoFs, numberOfLocalTestDoFs);

    // do loop over all local ansatz DoFs
    for (unsigned int i = 0; i < numberOfLocalAnsatzDoFs; ++i) {
      // do loop over all local test DoFs
      for (unsigned int j = 0; j < numberOfLocalTestDoFs; ++j) {
        const LocalAnsatzBaseFunctionType localAnsatzBaseFunction_i =
            localAnsatzBaseFunctionProvider_.provide(entity, i);
        const LocalTestBaseFunctionType localTestBaseFunction_j = localTestBaseFunctionProvider_.provide(entity, j);

        const RangeFieldType operator_i_j = localOperation_.operate(localAnsatzBaseFunction_i, localTestBaseFunction_j);

        // set local matrix
        ret[i][j] = operator_i_j;

      } // done loop over all local test DoFs

    } // done loop over all local ansatz DoFs

    return ret;
  }

private:
  const DiscreteAnsatzFunctionSpaceType& ansatzSpace_;
  const DiscreteTestFunctionSpaceType& testSpace_;
  const LocalOperationType localOperation_;
  const LocalAnsatzBaseFunctionProviderType localAnsatzBaseFunctionProvider_;
  const LocalTestBaseFunctionProviderType localTestBaseFunctionProvider_;

}; // end class FiniteElement

} // end of namespace Operator

} // end of namespace Functionals

} // end of namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
