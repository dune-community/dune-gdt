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
  void assemble( const LocalMatrixAssemblerType localMatrixAssembler, MatrixPtrType matrix )
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

    // storage
    LocalMatrixType localMatrix( ansatzSpace_.numMaxDoFs(), testSpace_.numMaxDoFs() );

    // do gridwalk
    const EntityIteratorType behindLastEntity = ansatzSpace_.end();
    for( EntityIteratorType entityIterator = ansatzSpace_.begin(); entityIterator != behindLastEntity; ++entityIterator )
    {
      // stuff
      const EntityType& entity = *entityIterator;

      // assemble local
      localMatrixAssembler.assembleLocal( entity, localMatrix );

      // add to matrix

    } // done gridwalk


  } // end method assemble

private:
  const AnsatzFunctionSpaceType& ansatzSpace_;
  const TestFunctionSpaceType& testSpace_;

}; // end class Generic

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMBLER_GENERIC_HH
