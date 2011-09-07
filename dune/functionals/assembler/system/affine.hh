#ifndef DUNE_FUNCTIONALS_ASSEMBLER_SYSTEM_AFFINE_HH
#define DUNE_FUNCTIONALS_ASSEMBLER_SYSTEM_AFFINE_HH

// dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune
{

namespace Functionals
{

namespace Assembler
{

namespace System
{

template< class AnsatzFunctionSpaceImp, class TestFunctionSpaceImp = AnsatzFunctionSpaceImp >
class Affine
{
public:

  typedef AnsatzFunctionSpaceImp
    AnsatzFunctionSpaceType;

  typedef TestFunctionSpaceImp
    TestFunctionSpaceType;

  typedef Affine< AnsatzFunctionSpaceImp, TestFunctionSpaceImp >
    ThisType;

  //! constructor
  Affine( const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( testSpace )
  {
  }

  //! constructor
  Affine( const AnsatzFunctionSpaceType& ansatzSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( ansatzSpace )
  {
  }

private:
  //! copy constructor
  Affine( const ThisType& other )
    : ansatzSpace_( other.ansatzSpace() ),
      testSpace_( other.testSpace() )
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

  template< class LocalMatrixAssemblerType, class MatrixType,
            class LocalVectorAssemblerType, class VectorType,
            class AffineShiftVectorType >
  void assembleSystem( const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& systemMatrix,
                       const LocalVectorAssemblerType& localVectorAssembler, VectorType& systemVector,
                       AffineShiftVectorType& affineShiftVector ) const
  {
    // some types
    typedef typename AnsatzFunctionSpaceType::GridPartType
      GridPartType;

    typedef typename AnsatzFunctionSpaceType::IteratorType
      EntityIteratorType;

    typedef typename AnsatzFunctionSpaceType::EntityType
      EntityType;

    typedef typename AnsatzFunctionSpaceType::RangeFieldType
      RangeFieldType;

    typedef Dune::DynamicMatrix< RangeFieldType >
      LocalMatrixType;

    typedef Dune::DynamicVector< RangeFieldType >
      LocalVectorType;

    typedef typename AnsatzFunctionSpaceType::ConstraintsType
      ConstraintsType;

    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;

    typedef typename AnsatzFunctionSpaceType::AffineShiftType
      AffineShiftType;

    // vector assembler for the affine shift
    typedef typename LocalMatrixAssemblerType::template LocalVectorAssembler< AffineShiftType >::Type
      LocalAffineShiftVectorAssemblerType;

    const LocalAffineShiftVectorAssemblerType localAffineShiftVectorAssembler = localMatrixAssembler.localVectorAssembler( ansatzSpace_.affineShift() );

    // common storage for all entities
    LocalMatrixType tmpLocalMatrix( ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType( 0.0 ) );
    LocalVectorType tmpLocalVector( testSpace_.map().maxLocalSize(), RangeFieldType( 0.0 ) );

    // do first gridwalk to assemble
    const EntityIteratorType lastEntity = ansatzSpace_.end();
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != lastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;

      localMatrixAssembler.assembleLocal( ansatzSpace_, testSpace_, entity, systemMatrix, tmpLocalMatrix );
      localVectorAssembler.assembleLocal( testSpace_, entity, systemVector, tmpLocalVector );
      localAffineShiftVectorAssembler.assembleLocal( testSpace_, entity, affineShiftVector, tmpLocalVector );

    } // done first gridwalk to assemble

    const ConstraintsType& constraints = ansatzSpace_.constraints();

    // do second gridwalk, to apply constraints
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != lastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;

      const LocalConstraintsType& localConstraints = constraints.local( entity );

      applyLocalMatrixConstraints( localConstraints, systemMatrix );
      applyLocalVectorConstraints( localConstraints, systemVector );

    } // done second gridwalk, to apply constraints

  } // end method assembleSystem

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  template< class LocalConstraintsType, class MatrixType >
  void applyLocalMatrixConstraints( const LocalConstraintsType& localConstraints, MatrixType& matrix ) const
  {
    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i )
    {
      for( unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j )
      {
        matrix[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)]
          = localConstraints.localMatrix(i,j);
      }
    }
  } // end applyLocalMatrixConstraints

  template< class LocalConstraintsType, class VectorType >
  void applyLocalVectorConstraints( const LocalConstraintsType& localConstraints, VectorType& vector ) const
  {
    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i )
    {
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
