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

/**
  \todo Add neumann boundary treatment
  \todo When adding neumann: think of numTmpObjectsRequired()!
  \todo Add vector assembler
  \todo Add penalty parameter
  **/
template <class LocalInnerOperatorImp, class LocalDirichletOperatorImp /*, class LocalNeumannOperatorImp*/>
class Matrix
{
public:
  typedef LocalInnerOperatorImp LocalInnerOperatorType;

  typedef LocalDirichletOperatorImp LocalDirichletOperatorType;

  //  typedef LocalNeumannOperatorImp
  //    LocalNeumannOperatorType;

  typedef Matrix<LocalInnerOperatorType, LocalDirichletOperatorType /*, LocalNeumannOperatorType*/> ThisType;

  typedef typename LocalInnerOperatorType::RangeFieldType RangeFieldType;

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
  Matrix(const LocalInnerOperatorType localInnerOperator, const LocalDirichletOperatorType localDirichletOperator /*,
          const LocalNeumannOperatorType localNeumannOperator*/)
    : localInnerOperator_(localInnerOperator)
    , localDirichletOperator_(localDirichletOperator) /*,
     localNeumannOperator_( localNeumannOperator )*/
  {
  }

private:
  //! copy constructor
  Matrix(const ThisType& other)
    : localInnerOperator_(other.localInnerOperator())
    , localDirichletOperator_(other.localDirichletOperator()) /*,
   localNeumannOperator_( other.localNeumannOperator() )*/
  {
  }

public:
  const LocalInnerOperatorType& localInnerOperator() const
  {
    return localInnerOperator_;
  }

  const LocalDirichletOperatorType& localDirichletOperator() const
  {
    return localDirichletOperator_;
  }

  //  const LocalNeumannOperatorType& localNeumannOperator() const
  //  {
  //    return localNeumannOperator_;
  //  }

  //  template< class InducingDiscreteFunctionType >
  //  typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type
  //    localVectorAssembler( const InducingDiscreteFunctionType& inducingDiscreteFunction ) const
  //  {
  //    typedef typename LocalVectorAssembler< InducingDiscreteFunctionType >::Type
  //      LocalVectorAssemblerType;

  //    return LocalVectorAssemblerType( localOperator_.localFunctional( inducingDiscreteFunction ) );
  //  }

  /**
    \todo Add neumann treatment here!
    **/
  std::vector<unsigned int> numTmpObjectsRequired() const
  {
    std::vector<unsigned int> ret(2, 0);
    // we require 4 tmp matrix in this local assembler
    ret[0] = 4;
    // the operator itself requires that much local matrices
    ret[1] = std::max(localInnerOperator_.numTmpObjectsRequired(), localDirichletOperator_.numTmpObjectsRequired());
    return ret;
  }

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class SystemMatrixType, class LocalMatrixType>
  void assembleLocal(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                     SystemMatrixType& systemMatrix,
                     std::vector<std::vector<LocalMatrixType>>& tmpLocalMatricesContainer) const
  {
    // get the local basefunction sets
    typedef typename AnsatzSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalAnsatzBaseFunctionSetType;

    const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSetEn = ansatzSpace.baseFunctionSet().local(entity);

    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTesBaseFunctionSetType;

    const LocalTesBaseFunctionSetType localTesBaseFunctionSetEn = testSpace.baseFunctionSet().local(entity);

    // check tmp local matrices
    assert(tmpLocalMatricesContainer.size() > 1);
    std::vector<LocalMatrixType>& tmpLocalMatrices = tmpLocalMatricesContainer[0];
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
        // get neighbouring entity
        const EntityPointerType neighbourPtr = intersection.outside();
        const EntityType& neighbour          = *neighbourPtr;

        // do visit only once
        if (gridPart.indexSet().index(entity) < gridPart.indexSet().index(neighbour)) {
          // get neighbouring local basefunction sets
          const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSetNe =
              ansatzSpace.baseFunctionSet().local(neighbour);
          const LocalTesBaseFunctionSetType localTesBaseFunctionSetNe = testSpace.baseFunctionSet().local(neighbour);

          localInnerOperator_.applyLocal(localAnsatzBaseFunctionSetEn,
                                         localAnsatzBaseFunctionSetNe,
                                         localTesBaseFunctionSetEn,
                                         localTesBaseFunctionSetNe,
                                         intersection,
                                         tmpLocalMatrices[0],
                                         tmpLocalMatrices[1],
                                         tmpLocalMatrices[2],
                                         tmpLocalMatrices[3],
                                         tmpLocalMatricesContainer[1]);

          // write local matrix to global (see below)
          addToMatrix(ansatzSpace, testSpace, entity, entity, tmpLocalMatrices[0], systemMatrix);
          addToMatrix(ansatzSpace, testSpace, entity, neighbour, tmpLocalMatrices[1], systemMatrix);
          addToMatrix(ansatzSpace, testSpace, neighbour, entity, tmpLocalMatrices[2], systemMatrix);
          addToMatrix(ansatzSpace, testSpace, neighbour, neighbour, tmpLocalMatrices[3], systemMatrix);
        } // done visit only once

      } // end if inner intersection
      else if (!intersection.neighbor() && intersection.boundary()) // if boundary intersection
      {
        const unsigned int boundaryId = intersection.boundaryId();

        //        // if dirichlet boundary intersection
        //        if( boundaryId == 2 )
        //        {

        localDirichletOperator_.applyLocal(localAnsatzBaseFunctionSetEn,
                                           localTesBaseFunctionSetEn,
                                           intersection,
                                           tmpLocalMatrices[0],
                                           tmpLocalMatricesContainer[1]);

        addToMatrix(ansatzSpace, testSpace, entity, entity, tmpLocalMatrices[0], systemMatrix);

        //        } // end if dirichlet boundary intersection
        //        else if( boundaryId == 3 ) // if neumann boundary intersection
        //        {

        //        } // end if neumann boundary intersection

      } // end if boundary intersection
    } // done loop over all intersections
  } // end method assembleLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class LocalMatrixType, class SystemMatrixType>
  void addToMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& ansatzEntity,
                   const EntityType& testEntity, const LocalMatrixType& localMatrix,
                   SystemMatrixType& systemMatrix) const
  {
    unsigned int rows = ansatzSpace.baseFunctionSet().local(ansatzEntity).size();
    unsigned int cols = testSpace.baseFunctionSet().local(testEntity).size();
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
        const unsigned int globalI = ansatzSpace.map().toGlobal(ansatzEntity, i);
        const unsigned int globalJ = testSpace.map().toGlobal(testEntity, j);

        systemMatrix[globalI][globalJ] += localMatrix[i][j];
      }
    }
  } // end method addToMatrix

  const LocalInnerOperatorType localInnerOperator_;
  const LocalDirichletOperatorType localDirichletOperator_;
  //  const LocalNeumannOperatorType localNeumannOperator;

}; // end class Matrix

} // end namespace Codim1

} // end namespace Local

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
