#ifndef DUNE_FUNCTIONALS_ASSEMBLER_GENERIC_HH
#define DUNE_FUNCTIONALS_ASSEMBLER_GENERIC_HH

// dune-functionals includes
#include <dune/functionals/common/localmatrix.hh>

namespace Dune
{

namespace Functionals
{

namespace Assembler
{

template< class AnsatzFunctionSpaceType, class TestFunctionSpaceType = AnsatzFunctionSpaceType >
class System
{
public:

  //! constructor
  System( const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( testSpace )
  {
  }

  //! constructor
  System( const AnsatzFunctionSpaceType& ansatzSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( ansatzSpace )
  {
  }

  template< class LocalMatrixAssemblerType, class MatrixPtrType >
  void assemble( const LocalMatrixAssemblerType localMatrixAssembler, MatrixPtrType matrixPtr )
  {
    // some types
    typedef typename AnsatzFunctionSpaceType::IteratorType
      EntityIteratorType;

    typedef typename AnsatzFunctionSpaceType::EntityType
      EntityType;

    typedef typename AnsatzFunctionSpaceType::RangeFieldType
      RangeFieldType;

    typedef Dune::Functionals::Common::LocalMatrix< RangeFieldType >
      LocalMatrixType;

    typedef typename AnsatzFunctionSpaceType::ConstraintsType
      ConstraintsType;

    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;

    // common storage for all entities
    LocalMatrixType localMatrix( ansatzSpace_.numMaxDoFs(), testSpace_.numMaxDoFs() );

    // do first gridwalk to assemble
    const EntityIteratorType behindLastEntity = ansatzSpace_.end();
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != behindLastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;

      localMatrixAssembler.assembleLocal( entity, localMatrix );

      addToMatrix( entity, localMatrix, matrixPtr );

    } // done first gridwalk to assemble

    const ConstraintsType constraints = ansatzSpace_.constraints();

    // do second gridwalk, to apply constraints
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != behindLastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;

      const LocalConstraintsType& localConstraints = constraints.local( entity );

      applyLocalMatrixConstraints( localConstraints, matrixPtr );

    } // done second gridwalk, to apply constraints

    // apply constraints

  } // end method assemble

private:

  template< class EntityType, class LocalMatrixType, class MatrixPtrType >
  void addToMatrix( const EntityType& entity, const LocalMatrixType& localMatrix, MatrixPtrType& matrixPtr )
  {
    for( int i = 0; i < ansatzSpace_.baseFunctionSet( entity ).numBaseFunctions(); ++i )
    {
      for( int j = 0; j < ansatzSpace_.baseFunctionSet( entity ).numBaseFunctions(); ++j )
      {
        const int globalI = ansatzSpace_.mapToGlobal( entity, i );
        const int globalJ = testSpace_.mapToGlobal( entity, j );

        matrixPtr->operator[](globalI)[globalJ] += localMatrix[i][j];
      }
    }
  } // end method addToMatrix

  template< class LocalConstraintsType, class MatrixPtrType >
  void applyLocalMatrixConstraints( const LocalConstraintsType& localConstraints, MatrixPtrType& matrixPtr )
  {
    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i )
    {
      for( unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j )
      {
        matrixPtr->operator[](localConstraints.rowDofs(i))[localConstraints.columnDofs(j)]
          = localConstraints.localMatrix(i,j);
      }
    }
  } // end applyLocalMatrixConstraints

  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;

}; // end class Generic

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMBLER_GENERIC_HH
