#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH

// std includes
#include <vector>

// local includes
#include "vector.hh"

namespace Dune {

namespace DetailedDiscretizations {

namespace Assembler {

namespace Local {

namespace Codim0 {

template <class LocalOperatorImp>
class Matrix
{
public:
  typedef LocalOperatorImp LocalOperatorType;

  typedef Matrix<LocalOperatorType> ThisType;

  typedef typename LocalOperatorType::RangeFieldType RangeFieldType;

  template <class InducingDiscreteFunctionType>
  class LocalVectorAssembler
  {
  private:
    typedef
        typename LocalOperatorType::template LocalFunctional<InducingDiscreteFunctionType>::Type InducingFunctionalType;

  public:
    typedef Dune::DetailedDiscretizations::Assembler::Local::Codim0::Vector<InducingFunctionalType> Type;
  };

  //! constructor
  Matrix(const LocalOperatorType localOperator)
    : localOperator_(localOperator)
  {
  }

private:
  //! copy constructor
  Matrix(const ThisType& other)
    : localOperator_(other.localOperator())
  {
  }

public:
  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

  template <class InducingDiscreteFunctionType>
  typename LocalVectorAssembler<InducingDiscreteFunctionType>::Type
  localVectorAssembler(const InducingDiscreteFunctionType& inducingDiscreteFunction) const
  {
    typedef typename LocalVectorAssembler<InducingDiscreteFunctionType>::Type LocalVectorAssemblerType;

    return LocalVectorAssemblerType(localOperator_.localFunctional(inducingDiscreteFunction));
  }

  std::vector<unsigned int> numTmpObjectsRequired() const
  {
    std::vector<unsigned int> ret(2, 0);
    // we require 1 tmp matrix in this local assembler
    ret[0] = 1;
    // the operator itself requires that much local matrices
    ret[1] = localOperator_.numTmpObjectsRequired();
    return ret;
  }

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class SystemMatrixType, class LocalMatrixType>
  void assembleLocal(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                     SystemMatrixType& systemMatrix,
                     std::vector<std::vector<LocalMatrixType>>& tmpLocalMatricesContainer) const
  {
    // get the local basefunctionsets
    typedef typename AnsatzSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalAnsatzBaseFunctionSetType;

    const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSet = ansatzSpace.baseFunctionSet().local(entity);

    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTestBaseFunctionSetType;

    const LocalTestBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local(entity);

    // check tmp local matrices
    assert(tmpLocalMatricesContainer.size() > 1);
    std::vector<LocalMatrixType>& tmpLocalMatrices = tmpLocalMatricesContainer[0];
    if (tmpLocalMatrices.size() < 1) {
      tmpLocalMatrices.resize(
          1, LocalMatrixType(ansatzSpace.map().maxLocalSize(), testSpace.map().maxLocalSize(), RangeFieldType(0.0)));
    }

    // write local operator application to tmpLocalMatrix
    localOperator_.applyLocal(
        localAnsatzBaseFunctionSet, localTestBaseFunctionSet, tmpLocalMatrices[0], tmpLocalMatricesContainer[1]);

    // write local matrix to global
    addToMatrix(ansatzSpace, testSpace, entity, tmpLocalMatrices[0], systemMatrix);
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class LocalMatrixType, class SystemMatrixType>
  void addToMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                   const LocalMatrixType& localMatrix, SystemMatrixType& systemMatrix) const
  {
    for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
      for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
        const unsigned int globalI = ansatzSpace.map().toGlobal(entity, i);
        const unsigned int globalJ = testSpace.map().toGlobal(entity, j);

        systemMatrix[globalI][globalJ] += localMatrix[i][j];
      }
    }
  } // end method addToMatrix

  const LocalOperatorType localOperator_;

}; // end class Matrix

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH
