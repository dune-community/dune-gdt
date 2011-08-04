#ifndef DUNE_FUNCTIONALS_ASSEMBLER_SYSTEM_AFFINE_HH
#define DUNE_FUNCTIONALS_ASSEMBLER_SYSTEM_AFFINE_HH

// dune-functionals includes
#include <dune/functionals/common/localmatrix.hh>
#include <dune/functionals/common/localvector.hh>

namespace Dune {

namespace Functionals {

namespace Assembler {

namespace System {

template <class AnsatzFunctionSpaceType, class TestFunctionSpaceType = AnsatzFunctionSpaceType>
class Affine
{
public:
  //! constructor
  Affine(const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
  {
  }

  //! constructor
  Affine(const AnsatzFunctionSpaceType& ansatzSpace)
    : ansatzSpace_(ansatzSpace)
    , testSpace_(ansatzSpace)
  {
  }

  template <class LocalMatrixAssemblerType, class MatrixType, class LocalVectorAssemblerType, class VectorType,
            class AffineShiftVectorType>
  void assembleSystem(const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& matrix,
                      const LocalVectorAssemblerType& localVectorAssembler, VectorType& vector,
                      AffineShiftVectorType& affineShiftVector)
  {
    // some types
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;

    typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;

    typedef typename EntityIteratorType::Entity EntityType;

    typedef typename AnsatzFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef Dune::Functionals::Common::LocalMatrix<RangeFieldType> LocalMatrixType;

    typedef Dune::Functionals::Common::LocalVector<RangeFieldType> LocalVectorType;

    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;

    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;

    typedef typename AnsatzFunctionSpaceType::AffineShiftType AffineShiftType;

    //    // vector assembler for the affine shift
    //    typedef typename LocalMatrixAssemblerType::template LocalVectorAssembler< AffineShiftType >::Type
    //      LocalAffineShiftVectorAssemblerType;

    //    const LocalAffineShiftVectorAssemblerType localAffineShiftVectorAssembler =
    //    localMatrixAssembler.localVectorAssembler( ansatzSpace_.affineShift() );

    // common storage for all entities
    LocalMatrixType tmpLocalMatrix(ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize());
    LocalVectorType tmpLocalVector(testSpace_.map().maxLocalSize());

    // do first gridwalk to assemble
    const EntityIteratorType lastEntity = ansatzSpace_.gridPart().template end<0>();
    for (EntityIteratorType entityIterator = ansatzSpace_.gridPart().template begin<0>(); entityIterator != lastEntity;
         ++entityIterator) {
      const EntityType& entity = *entityIterator;

      localMatrixAssembler.assembleLocal(ansatzSpace_, testSpace_, entity, matrix, tmpLocalMatrix);
      localVectorAssembler.assembleLocal(testSpace_, entity, vector, tmpLocalVector);
      //      localAffineShiftVectorAssembler.assembleLocal( testSpace_, entity, affineShiftVector, localVector );

    } // done first gridwalk to assemble

    const ConstraintsType constraints = ansatzSpace_.constraints();

    // do second gridwalk, to apply constraints
    for (EntityIteratorType entityIterator = ansatzSpace_.gridPart().template begin<0>(); entityIterator != lastEntity;
         ++entityIterator) {
      const EntityType& entity = *entityIterator;

      const LocalConstraintsType& localConstraints = constraints.local(entity);

      applyLocalMatrixConstraints(localConstraints, matrix);
      applyLocalVectorConstraints(localConstraints, vector);

    } // done second gridwalk, to apply constraints

  } // end method assembleSystem

private:
  template <class LocalConstraintsType, class MatrixType>
  void applyLocalMatrixConstraints(const LocalConstraintsType& localConstraints, MatrixType& matrix)
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j) {
        matrix[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)] = localConstraints.localMatrix(i, j);
      }
    }
  } // end applyLocalMatrixConstraints

  template <class LocalConstraintsType, class VectorType>
  void applyLocalVectorConstraints(const LocalConstraintsType& localConstraints, VectorType& vector)
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      vector[localConstraints.rowDofs(i)] = 0.0;
    }
  } // end applyLocalVectorConstraints

  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;

}; // end class Affine

} // end namespace System

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMBLER_SYSTEM_AFFINE_HH
