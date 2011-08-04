#ifndef DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH
#define DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_MATRIX_HH

// local includes
#include "vector.hh"

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

  typedef Matrix< LocalOperatorType >
    ThisType;

  typedef typename LocalOperatorType::RangeFieldType
    RangeFieldType;

  template< class InducingDiscreteFunctionType >
  class LocalVectorAssembler
  {
  private:
    typedef typename LocalOperatorType::template LocalFunctional< InducingDiscreteFunctionType >::Type
      InducingFunctionalType;

  public:
    typedef Dune::Functionals::Assembler::Local::Codim0::Vector< InducingFunctionalType >
      Type;
  };

  //! constructor
  Matrix( const LocalOperatorType& localOperator )
    : localOperator_( localOperator )
  {
  }

private:
  //! copy constructor
  Matrix( const ThisType& other )
    : localOperator_( other.localOperator() )
  {
  }

public:
  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

  template< class InducingDiscreteFunctionType >
  typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type localVectorAssembler( const InducingDiscreteFunctionType& inducingDiscreteFunction ) const
  {
    typedef typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type
      LocalVectorAssemblerType;

    return LocalVectorAssemblerType( localOperator_.localFunctional( inducingDiscreteFunction ) );
  }

  template< class AnsatzSpaceType,
            class TestSpaceType,
            class EntityType,
            class MatrixType,
            class LocalMatrixType >
  void assembleLocal( const AnsatzSpaceType& ansatzSpace,
                      const TestSpaceType& testSpace,
                      const EntityType& entity,
                      MatrixType& matrix,
                      LocalMatrixType& tmpLocalMatrix ) const
  {
    // write local operator application to tmpLocalMatrix
    localOperator_.applyLocal( ansatzSpace.baseFunctionSet().local( entity ), testSpace.baseFunctionSet().local( entity ), tmpLocalMatrix );

    // write local matrix to global
    addToMatrix( ansatzSpace, testSpace, entity, tmpLocalMatrix, matrix );
  }

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

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
    for( unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local( entity ).size(); ++i )
    {
      for( unsigned int j = 0; j < testSpace.baseFunctionSet().local( entity ).size(); ++j )
      {
        const unsigned int globalI = ansatzSpace.map().toGlobal( entity, i );
        const unsigned int globalJ = testSpace.map().toGlobal( entity, j );

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
