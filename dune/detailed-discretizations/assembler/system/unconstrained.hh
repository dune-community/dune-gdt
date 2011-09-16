#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH

// std includes
#include <vector>

// dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune
{

namespace DetailedDiscretizations
{

namespace Assembler
{

namespace System
{

template< class AnsatzFunctionSpaceImp, class TestFunctionSpaceImp = AnsatzFunctionSpaceImp >
class Unconstrained
{
public:

  typedef AnsatzFunctionSpaceImp
    AnsatzFunctionSpaceType;

  typedef TestFunctionSpaceImp
    TestFunctionSpaceType;

  typedef Unconstrained< AnsatzFunctionSpaceImp, TestFunctionSpaceImp >
    ThisType;

  //! constructor
  Unconstrained( const AnsatzFunctionSpaceType& ansatzSpace, const TestFunctionSpaceType& testSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( testSpace )
  {
  }

  //! constructor
  Unconstrained( const AnsatzFunctionSpaceType& ansatzSpace )
    : ansatzSpace_( ansatzSpace ),
      testSpace_( ansatzSpace )
  {
  }

private:
  //! copy constructor
  Unconstrained( const ThisType& other )
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
            class LocalVectorAssemblerType, class VectorType >
  void assembleSystem( const LocalMatrixAssemblerType& localMatrixAssembler, MatrixType& systemMatrix,
                       const LocalVectorAssemblerType& localVectorAssembler, VectorType& systemVector ) const
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

    // common storage for all entities
    std::vector< LocalMatrixType > tmpLocalMatrices( 4, LocalMatrixType( ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType( 0.0 ) ) );
    LocalVectorType tmpLocalVector( testSpace_.map().maxLocalSize(), RangeFieldType( 0.0 ) );

    // do gridwalk to assemble
    const EntityIteratorType lastEntity = ansatzSpace_.end();
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != lastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;

      localMatrixAssembler.assembleLocal( ansatzSpace_, testSpace_, entity, systemMatrix, tmpLocalMatrix );
      localVectorAssembler.assembleLocal( testSpace_, entity, systemVector, tmpLocalVector );

    } // done gridwalk to assemble

  } // end method assembleSystem

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;

}; // end class Unconstrained

} // end namespace System

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_UNCONSTRAINED_HH
