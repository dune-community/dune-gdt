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
  void assemble(const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& systemMatrix,
                const LocalVectorAssemblerType& localVectorAssembler, VectorType& systemVector) const
  {
    assembleSystem(localMatrixAssembler, systemMatrix, localVectorAssembler, systemVector);
    applyConstraints(systemMatrix, systemVector);
  } // void assemble()

  template <class LocalMatrixAssemblerType, class MatrixType, class LocalVectorAssemblerType, class VectorType>
  void assembleSystem(const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& systemMatrix,
                      const LocalVectorAssemblerType& localVectorAssembler, VectorType& systemVector) const
  {
    // some types
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
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

    // walk the grid to assemble
    for (EntityIteratorType entityIterator = ansatzSpace_.gridPart().template begin<0>();
         entityIterator != ansatzSpace_.gridPart().template end<0>();
         ++entityIterator) {
      const EntityType& entity = *entityIterator;
      localMatrixAssembler.assembleLocal(ansatzSpace_, testSpace_, entity, systemMatrix, tmpLocalMatricesContainer);
      localVectorAssembler.assembleLocal(testSpace_, entity, systemVector, tmpLocalVectorsContainer);
    } // walk the grid to assemble
  } // void assembleSystem()

  template <class MatrixType, class VectorType>
  void applyConstraints(MatrixType& matrix, VectorType& vector) const
  {
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    // walk the grid to apply constraints
    const ConstraintsType& constraints = ansatzSpace_.constraints();
    for (EntityIteratorType entityIterator = ansatzSpace_.gridPart().template begin<0>();
         entityIterator != ansatzSpace_.gridPart().template end<0>();
         ++entityIterator) {
      const EntityType& entity                     = *entityIterator;
      const LocalConstraintsType& localConstraints = constraints.local(entity);
      applyLocalMatrixConstraints(localConstraints, matrix);
      applyLocalVectorConstraints(localConstraints, vector);

    } // walk the grid to apply constraints

  } // void applyConstraints()


private:
  ThisType& operator=(const ThisType&);
  Constrained(const ThisType&);

  template <class LocalConstraintsType, class MatrixType>
  void applyLocalMatrixConstraints(const LocalConstraintsType& localConstraints, MatrixType& matrix) const
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      const unsigned int rowDof = localConstraints.rowDofs(i);
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j) {
        matrix.set(rowDof, localConstraints.columnDofs(j), localConstraints.localMatrix(i, j));
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
