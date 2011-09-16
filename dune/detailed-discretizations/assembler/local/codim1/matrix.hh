#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH

// std includes
#include <vector>

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
                     MatrixType& matrix, std::vector<LocalMatrixType>& tmpLocalMatrices) const
  {
    // get the local basefunction sets
    typedef typename AnsatzSpaceType::BaseFunctionSet::Local LocalAnsatzBaseFunctionSetType;

    const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSetEntity = ansatzSpace.baseFunctionSet().local(entity);
    const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSetNeighbour =
        ansatzSpace.baseFunctionSet().local(neighbour);

    typedef typename TestSpaceType::BaseFunctionSet::Local LocalTesBaseFunctionSetType;

    const LocalTesBaseFunctionSetType localTesBaseFunctionSetEntity    = testSpace.baseFunctionSet().local(entity);
    const LocalTesBaseFunctionSetType localTesBaseFunctionSetNeighbour = testSpace.baseFunctionSet().local(neighbour);

    // check, if we have enough tmp matrices
    if (tmpLocalMatrices.size() < 4) {
      tmpLocalMatrices.resize(
          4, LocalMatrixType(ansatzSpace.map().maxLocalSize(), testSpace.map().maxLocalSize(), RangeFieldType(0.0)));
    }

    // some types
    typedef typename AnsatzSpaceType::GridPartType GridPartType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    typedef typename IntersectionType::EntityPointer EntityPointerType;

    const GridPartType& gridPart = ansatzSpace.gridPart();

    const IntersectionIteratorType lastIntersection = gridPart.iend(entity);

    // do loop over all intersections
    for (IntersectionIteratorType intIt = entity.ibegin(); intIt != lastIntersection; ++intIt) {
      const IntersectionType& intersection = *intIt;

      // if inner intersection
      if (intersection.neighbor() && !intersection.boundary()) {
        const EntityPointerType neighbourPtr = intersection.outside();
        const EntityType& neighbour          = *neighbourPtr;

        innerLocalOperator_.applyLocal(localAnsatzBaseFunctionSetEntity,
                                       localAnsatzBaseFunctionSetNeighbour,
                                       localTesBaseFunctionSetEntity,
                                       localTesBaseFunctionSetNeighbour,
                                       intersection,
                                       tmpLocalMatrices[0],
                                       tmpLocalMatrices[1],
                                       tmpLocalMatrices[2],
                                       tmpLocalMatrices[3]);

        // write local matrix to global
        addToMatrix(ansatzSpace, testSpace, entity, neighbour, tmpLocalMatrix, matrix);

      } // end if inner intersection
      else if (!intersection.neighbor() && intersection.boundary()) // if boundary intersection
      {
        const unsigned int boundaryId = intersection.boundaryId();

        // if dirichlet boundary intersection
        if (boundaryId == 2) {

        } // end if dirichlet boundary intersection
        else if (boundaryId == 3) // if neumann boundary intersection
        {

        } // end if neumann boundary intersection
      } // end if boundary intersection
    } // done loop over all intersections
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class LocalMatrixType, class MatrixType>
  void addToMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                   const EntityType& neighbour, const LocalMatrixType& localMatrix, MatrixType& matrix) const
  {
    unsigned int rows = ansatzSpace.baseFunctionSet().local(entity).size();
    unsigned int cols = testSpace.baseFunctionSet().local(entity).size();
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
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
