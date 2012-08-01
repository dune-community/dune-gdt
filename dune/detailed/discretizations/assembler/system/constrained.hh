#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_AFFINE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_AFFINE_HH

// std includes
#include <vector>

// dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace System {

template <class AnsatzFunctionSpaceImp, class TestFunctionSpaceImp = AnsatzFunctionSpaceImp>
class Constrained
{
public:
  typedef AnsatzFunctionSpaceImp AnsatzFunctionSpaceType;

  typedef TestFunctionSpaceImp TestFunctionSpaceType;

  typedef Constrained<AnsatzFunctionSpaceImp, TestFunctionSpaceImp> ThisType;

  //! constructor
  Constrained(const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
  {
  }

  //! constructor
  Constrained(const AnsatzFunctionSpaceType& ansatzSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(ansatzSpace)
  {
  }

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
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;

    typedef typename AnsatzFunctionSpaceType::IteratorType EntityIteratorType;

    typedef typename AnsatzFunctionSpaceType::EntityType EntityType;

    typedef typename AnsatzFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef Dune::DynamicMatrix<RangeFieldType> LocalMatrixType;

    typedef Dune::DynamicVector<RangeFieldType> LocalVectorType;

    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;

    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;

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

    // do first gridwalk to assemble
    const EntityIteratorType lastEntity = ansatzSpace_.end();
    for (EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != lastEntity; ++entityIterator) {
      const EntityType& entity = *entityIterator;

      localMatrixAssembler.assembleLocal(ansatzSpace_, testSpace_, entity, systemMatrix, tmpLocalMatricesContainer);
      localVectorAssembler.assembleLocal(testSpace_, entity, systemVector, tmpLocalVectorsContainer);

    } // done first gridwalk to assemble

    // do second gridwalk, to apply constraints
    const ConstraintsType& constraints = ansatzSpace_.constraints();

    for (EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != lastEntity; ++entityIterator) {
      const EntityType& entity = *entityIterator;

      const LocalConstraintsType& localConstraints = constraints.local(entity);

      applyLocalMatrixConstraints(localConstraints, systemMatrix);
      applyLocalVectorConstraints(localConstraints, systemVector);

    } // done second gridwalk, to apply constraints

  } // end method assembleSystem

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  //! copy constructor
  Constrained(const ThisType&);

  template <class LocalConstraintsType, class MatrixType>
  void applyLocalMatrixConstraints(const LocalConstraintsType& localConstraints, MatrixType& matrix) const
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j) {
        matrix.set(localConstraints.rowDofs(i), localConstraints.columnDofs(j), localConstraints.localMatrix(i, j));
      }
    }
  } // end applyLocalMatrixConstraints

  template <class LocalConstraintsType, class VectorType>
  void applyLocalVectorConstraints(const LocalConstraintsType& localConstraints, VectorType& vector) const
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      vector.set(localConstraints.rowDofs(i), 0.0);
    }
  } // end applyLocalVectorConstraints

  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;
}; // end class Constrained

} // end namespace System

} // end namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_AFFINE_HH
