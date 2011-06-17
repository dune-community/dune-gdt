#ifndef DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH
#define DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH

// dune-functionals includes
#include <dune/functionals/common/localmatrix.hh>

namespace Dune
{

namespace Functionals
{

namespace Assembler
{

namespace Local
{

namespace Codim0
{

template< class LocalOperatorImp >
class Matrix
{
public:

  typedef LocalOperatorImp
    LocalOperatorType;

  typedef typename LocalOperatorType::RangeFieldType
    RangeFieldType;

  //! constructor
  Matrix( const LocalOperatorType& localOperator )
    : localOperator_( localOperator )
  {
  }

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

  template< class AnsatzSpaceType,
            class TestSpaceType,
            class EntityType,
            class MatrixType,
            class LocalMatrixType = Dune::Functionals::Common::LocalMatrix< RangeFieldType > >
  void assembleLocal( const AnsatzSpaceType& ansatzSpace,
                      const TestSpaceType& testSpace,
                      const EntityType& entity,
                      MatrixType& matrix,
                      LocalMatrixType localMatrix ) const
  {
    // write local operator application to tmpLocalMatrix
    localOperator_.applyLocal( ansatzSpace.localBaseFunctionSet( entity ), testSpace.localBaseFunctionSet( entity ), localMatrix );

    // write local matrix to global
    addToMatrix( ansatzSpace, testSpace, entity, localMatrix, matrix );
  }

private:

  template< class AnsatzSpaceType,
            class TestSpaceType,
            class EntityType,
            class LocalMatrixType,
            class MatrixType >
  void addToMatrix( const AnsatzSpaceType& ansatzSpace,
                    const TestSpaceType& testSpace,
                    const EntityType& entity,
                    const LocalMatrixType& localMatrix,
                    MatrixType& matrix ) const
  {
    for( int i = 0; i < ansatzSpace.baseFunctionSet( entity ).numBaseFunctions(); ++i )
    {
      for( int j = 0; j < testSpace.baseFunctionSet( entity ).numBaseFunctions(); ++j )
      {
        const int globalI = ansatzSpace.mapToGlobal( entity, i );
        const int globalJ = testSpace.mapToGlobal( entity, j );

        matrix[globalI][globalJ] += localMatrix[i][j];
      }
    }
  } // end method addToMatrix

  const LocalOperatorType& localOperator_;

}; // end class Matrix

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH
