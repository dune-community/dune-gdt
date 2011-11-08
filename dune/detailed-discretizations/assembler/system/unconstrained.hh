#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH

// std includes
#include <vector>

// dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace Assembler {

namespace System {

template <class AnsatzFunctionSpaceImp, class TestFunctionSpaceImp = AnsatzFunctionSpaceImp>
class Unconstrained
{
public:
  typedef AnsatzFunctionSpaceImp AnsatzFunctionSpaceType;

  typedef TestFunctionSpaceImp TestFunctionSpaceType;

  typedef Unconstrained<AnsatzFunctionSpaceImp, TestFunctionSpaceImp> ThisType;

  //! constructor
  Unconstrained(const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
  {
  }

  //! constructor
  Unconstrained(const AnsatzFunctionSpaceType& ansatzSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(ansatzSpace)
  {
  }

private:
  //! copy constructor
  Unconstrained(const ThisType& other)
    : ansatzSpace_(other.ansatzSpace())
    , testSpace_(other.testSpace())
  {
  }

public:
  const AnsatzFunctionSpaceType& ansatzSpace()
  {
    return ansatzSpace_;
  }

  const TestFunctionSpaceType& testSpace()
  {
    return testSpace_;
  }

  template <class LocalMatrixAssemblerType, class MatrixType, class LocalVectorAssemblerType, class VectorType>
  void assembleSystem(const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& systemMatrix,
                      const LocalVectorAssemblerType& localVectorAssembler, VectorType& systemVector) const
  {
    // some types
    typedef typename AnsatzFunctionSpaceType::GridElementIteratorType GridElementIteratorType;

    typedef typename AnsatzFunctionSpaceType::GridElementType GridElementType;

    typedef typename AnsatzFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef Dune::DynamicMatrix<RangeFieldType> LocalMatrixType;

    typedef Dune::DynamicVector<RangeFieldType> LocalVectorType;

    // common tmp storage for all entities
    std::vector<unsigned int> numberTmpMatrices = localMatrixAssembler.numTmpObjectsRequired();
    std::vector<LocalMatrixType> tmpLocalAssemblerMatrices(
        numberTmpMatrices[0],
        LocalMatrixType(ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<LocalMatrixType> tmpLocalOperatorMatrices(
        numberTmpMatrices[1],
        LocalMatrixType(ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<std::vector<LocalMatrixType>> tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);

    std::vector<unsigned int> numberTmpVectors = localVectorAssembler.numTmpObjectsRequired();
    std::vector<LocalVectorType> tmpLocalAssemblerVectors(
        numberTmpVectors[0], LocalVectorType(testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<LocalVectorType> tmpLocalFunctionalVectors(
        numberTmpVectors[1], LocalVectorType(testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<std::vector<LocalVectorType>> tmpLocalVectorsContainer;
    tmpLocalVectorsContainer.push_back(tmpLocalAssemblerVectors);
    tmpLocalVectorsContainer.push_back(tmpLocalFunctionalVectors);

    // do gridwalk to assemble
    const GridElementIteratorType lastElement = ansatzSpace_.gridElementEnd();
    for (GridElementIteratorType elementIterator = ansatzSpace_.gridElementBegin(); elementIterator != lastElement;
         ++elementIterator) {
      const GridElementType& element = *elementIterator;

      localMatrixAssembler.assembleLocal(ansatzSpace_, testSpace_, element, systemMatrix, tmpLocalMatricesContainer);
      localVectorAssembler.assembleLocal(testSpace_, element, systemVector, tmpLocalVectorsContainer);

    } // done gridwalk to assemble

  } // end method assembleSystem

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;

}; // end class Unconstrained

} // end namespace System

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH
