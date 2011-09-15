#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH

// local includes
//#include "vector.hh"

namespace Dune {

namespace DetailedDiscretizations {

namespace Assembler {

namespace Local {

namespace Codim1 {

template <class InnerLocalOperatorImp, class DirichletLocalOperatorImp, class NeumannLocalOperatorImp>
class Matrix
{
public:
  typedef InnerLocalOperatorImp InnerLocalOperatorType;

  typedef DirichletLocalOperatorImp DirichletLocalOperatorType;

  typedef NeumannLocalOperatorImp NeumannLocalOperatorType;

  typedef Matrix<InnerLocalOperatorType, DirichletLocalOperatorType, NeumannLocalOperatorType> ThisType;

  typedef typename InnerLocalOperatorType::RangeFieldType RangeFieldType;

  //  template< class InducingDiscreteFunctionType >
  //  class LocalVectorAssembler
  //  {
  //  private:
  //    typedef typename LocalOperatorType::template LocalFunctional< InducingDiscreteFunctionType >::Type
  //      InducingFunctionalType;

  //  public:
  //    typedef Dune::Functionals::Assembler::Local::Codim0::Vector< InducingFunctionalType >
  //      Type;
  //  };

  //! constructor
  Matrix(const InnerLocalOperatorType innerLocalOperator, const DirichletLocalOperatorType dirichletLocalOperator,
         const NeumannLocalOperatorType neumannLocalOperator)
    : innerLocalOperator_(innerLocalOperator)
    , dirichletLocalOperator_(dirichletLocalOperator)
    , neumannLocalOperator_(neumannLocalOperator)
  {
  }

private:
  //! copy constructor
  Matrix(const ThisType& other)
    : localOperator_(other.localOperator())
  {
  }

public:
  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

  //  template< class InducingDiscreteFunctionType >
  //  typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type
  //    localVectorAssembler( const InducingDiscreteFunctionType& inducingDiscreteFunction ) const
  //  {
  //    typedef typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type
  //      LocalVectorAssemblerType;

  //    return LocalVectorAssemblerType( localOperator_.localFunctional( inducingDiscreteFunction ) );
  //  }

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class MatrixType, class LocalMatrixType>
  void assembleLocal(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                     MatrixType& matrix, LocalMatrixType& tmpLocalMatrix) const
  {
    typedef typename AnsatzSpaceType::GridPartType GridPartType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    typedef typename IntersectionType::EntityPointer EntityPointerType;

    const GridPartType& gridPart = ansatzSpace.gridPart();

    const IntersectionIteratorType lastIntersection = gridPart.iend(entity);

    for (IntersectionIteratorType intIt = entity.ibegin(); intIt != lastIntersection; ++intIt) {
      const IntersectionType& intersection = *intIt;

      // if inner intersection
      if (intersection.neighbor() && !intersection.boundary()) {
        const EntityPointerType neighbourPtr = intersection.outside();
        const EntityType& neighbour          = *neighbourPtr;

        innerLocalOperator_.applyLocal(intersection,
                                       ansatzSpace.baseFunctionSet().local(entity),
                                       testSpace.baseFunctionSet().local(neighbour),
                                       tmpLocalMatrix);

        // write local matrix to global
        addToMatrix(ansatzSpace, testSpace, entity, neighbour, tmpLocalMatrix, matrix);


      } else if (!intersection.neighbor() && intersection.boundary()) // else boundary
      {
        const unsigned int boundaryId = intersection.boundaryId();

        // if dirichlet
        if (boundaryId == 2) {

        } else if (boundaryId == 3) // else neumann
        {
        }

      } // end if boundary intersection
    }
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class LocalMatrixType, class MatrixType>
  void addToMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                   const EntityType& neighbour, const LocalMatrixType& localMatrix, MatrixType& matrix) const
  {
    for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
      for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
        const unsigned int globalI = ansatzSpace.map().toGlobal(entity, i);
        const unsigned int globalJ = testSpace.map().toGlobal(neighbour, j);

        matrix[globalI][globalJ] += localMatrix[i][j];
      }
    }
  } // end method addToMatrix

  const InnerLocalOperatorType innerLocalOperator_;
  const DirichletLocalOperatorType dirichletLocalOperator_;
  const NeumannLocalOperatorType neumannLocalOperator_;

}; // end class Matrix

} // end namespace Codim1

} // end namespace Local

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
