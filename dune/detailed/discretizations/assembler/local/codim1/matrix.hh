#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH

// std includes
#include <vector>

// local includes
//#include "vector.hh"

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Local {

namespace Codim1 {

template< class LocalOperatorImp >
class Inner
{
public:
  typedef LocalOperatorImp LocalOperatorType;

  typedef Inner< LocalOperatorType > ThisType;

  typedef typename LocalOperatorType::RangeFieldType RangeFieldType;

  Inner(const LocalOperatorType& localOperator)
    : localOperator_( localOperator )
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const unsigned int numTmpObjectsRequired_ = 4;

public:
  std::vector< unsigned int > numTmpObjectsRequired() const
  {
    std::vector< unsigned int > ret(2, 0);
    ret[0] = numTmpObjectsRequired_;
    ret[1] = localOperator_.numTmpObjectsRequired();
    return ret;
  } // std::vector< unsigned int > numTmpObjectsRequired() const

  template< class IntersectionType,
            class InnerAnsatzSpaceType,
            class InnerTestSpaceType,
            class OuterAnsatzSpaceType,
            class OuterTestSpaceType,
            class MatrixBackendType,
            class LocalMatrixType >
  void assembleLocal(const IntersectionType& intersection,
                     const InnerAnsatzSpaceType& innerAnsatzSpace,
                     const InnerTestSpaceType& innerTestSpace,
                     const OuterAnsatzSpaceType& outerAnsatzSpace,
                     const OuterTestSpaceType& outerTestSpace,
                     MatrixBackendType& innerInnerMatrix,
                     MatrixBackendType& outerOuterMatrix,
                     MatrixBackendType& innerOuterMatrix,
                     MatrixBackendType& outerInnerMatrix,
                     std::vector< std::vector< LocalMatrixType > >& tmpLocalMatricesContainer) const
  {
    // preparations
    assert(intersection.neighbor() && !intersection.boundary());
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    typedef typename IntersectionType::Entity EntityType;
    typedef typename InnerAnsatzSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType InnerAnsatzBaseFunctionSetType;
    typedef typename InnerTestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType InnerTestBaseFunctionSetType;
    typedef typename OuterAnsatzSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType OuterAnsatzBaseFunctionSetType;
    typedef typename OuterTestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType OuterTestBaseFunctionSetType;
    // get inside entity and basefunctionsets
    const EntityPointerType insideEntityPtr = intersection.inside();
    const EntityType& insideEntity = *insideEntityPtr;
    const InnerAnsatzBaseFunctionSetType innerAnsatzBaseFunctionSet = innerAnsatzSpace.baseFunctionSet().local(insideEntity);
    const InnerTestBaseFunctionSetType innerTestBaseFunctionSet = innerTestSpace.baseFunctionSet().local(insideEntity);
    // get outside neighbor and basefunctionsets
    const EntityPointerType outsideNeighborPtr = intersection.outside();
    const EntityType& outsideNeighbor = *outsideNeighborPtr;
    const OuterAnsatzBaseFunctionSetType outerAnsatzBaseFunctionSet = outerAnsatzSpace.baseFunctionSet().local(outsideNeighbor);
    const OuterTestBaseFunctionSetType outerTestBaseFunctionSet = outerTestSpace.baseFunctionSet().local(outsideNeighbor);
    // ensure enough tmp local matrices
    assert(tmpLocalMatricesContainer.size() > 1);
    std::vector< LocalMatrixType >& tmpLocalMatrices = tmpLocalMatricesContainer[0];
    if( tmpLocalMatrices.size() < numTmpObjectsRequired_ ) {
      tmpLocalMatrices.resize(numTmpObjectsRequired_, LocalMatrixType(std::max(innerAnsatzSpace.map().maxLocalSize(), outerAnsatzSpace.map().maxLocalSize()),
                                                                      std::max(innerTestSpace.map().maxLocalSize(), outerTestSpace.map().maxLocalSize()),
                                                                      RangeFieldType(0.0)));
    } // ensure enough tmp local matrices
    // apply local operator
    localOperator_.applyLocal(innerAnsatzBaseFunctionSet,
                              innerTestBaseFunctionSet,
                              outerAnsatzBaseFunctionSet,
                              outerTestBaseFunctionSet,
                              intersection,
                              tmpLocalMatrices[0], // inside/inside
                              tmpLocalMatrices[1], // outside/outside
                              tmpLocalMatrices[2], // inside/outside
                              tmpLocalMatrices[3], // outside/inside
                              tmpLocalMatricesContainer[1]);
    // write local matrices to global (see below)
    addToMatrix(innerAnsatzSpace, innerTestSpace, insideEntity,    insideEntity,    tmpLocalMatrices[0], innerInnerMatrix);
    addToMatrix(outerAnsatzSpace, outerTestSpace, outsideNeighbor, outsideNeighbor, tmpLocalMatrices[1], outerOuterMatrix);
    addToMatrix(innerAnsatzSpace, outerTestSpace, insideEntity,    outsideNeighbor, tmpLocalMatrices[2], innerOuterMatrix);
    addToMatrix(outerAnsatzSpace, innerTestSpace, outsideNeighbor, insideEntity,    tmpLocalMatrices[3], outerInnerMatrix);
  } // void assembleLocal() const

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  template< class AnsatzSpaceType,
            class TestSpaceType,
            class EntityType,
            class LocalMatrixType,
            class SystemMatrixType >
  void addToMatrix( const AnsatzSpaceType& ansatzSpace,
                    const TestSpaceType& testSpace,
                    const EntityType& ansatzEntity,
                    const EntityType& testEntity,
                    const LocalMatrixType& localMatrix,
                    SystemMatrixType& systemMatrix ) const
  {
    unsigned int rows = ansatzSpace.baseFunctionSet().local( ansatzEntity ).size();
    unsigned int cols = testSpace.baseFunctionSet().local( testEntity ).size();
    for( unsigned int i = 0; i < rows; ++i )
    {
      for( unsigned int j = 0; j < cols; ++j )
      {
        const unsigned int globalI = ansatzSpace.map().toGlobal( ansatzEntity, i );
        const unsigned int globalJ = testSpace.map().toGlobal( testEntity, j );

        systemMatrix.add(globalI, globalJ, localMatrix[i][j]);
      }
    }
  } // end method addToMatrix

  const LocalOperatorType& localOperator_;
}; // end class Inner

template< class LocalOperatorImp >
class Boundary
{
public:
  typedef LocalOperatorImp LocalOperatorType;

  typedef Boundary< LocalOperatorType > ThisType;

  typedef typename LocalOperatorType::RangeFieldType RangeFieldType;

  Boundary(const LocalOperatorType& localOperator)
    : localOperator_( localOperator )
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const unsigned int numTmpObjectsRequired_ = 4;

public:
  std::vector< unsigned int > numTmpObjectsRequired() const
  {
    std::vector< unsigned int > ret(2, 0);
    ret[0] = numTmpObjectsRequired_;
    ret[1] = localOperator_.numTmpObjectsRequired();
    return ret;
  } // std::vector< unsigned int > numTmpObjectsRequired() const

  template< class IntersectionType,
            class AnsatzSpaceType,
            class TestSpaceType,
            class MatrixBackendType,
            class LocalMatrixType >
  void assembleLocal(const IntersectionType& intersection,
                     const AnsatzSpaceType& ansatzSpace,
                     const TestSpaceType& testSpace,
                     MatrixBackendType& matrix,
                     std::vector< std::vector< LocalMatrixType > >& tmpLocalMatricesContainer) const
  {
    // preparations
    assert(intersection.boundary());
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    typedef typename IntersectionType::Entity EntityType;
    typedef typename AnsatzSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalAnsatzBaseFunctionSetType;
    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTestBaseFunctionSetType;
    // get inside entity and basefunctionsets
    const EntityPointerType entityPtr = intersection.inside();
    const EntityType& entity = *entityPtr;
    const LocalAnsatzBaseFunctionSetType localAnsatzBaseFunctionSet = ansatzSpace.baseFunctionSet().local(entity);
    const LocalTestBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local(entity);
    // ensure enough tmp local matrices
    assert(tmpLocalMatricesContainer.size() > 1);
    std::vector< LocalMatrixType >& tmpLocalMatrices = tmpLocalMatricesContainer[0];
    if (tmpLocalMatrices.size() < numTmpObjectsRequired_) {
      tmpLocalMatrices.resize(numTmpObjectsRequired_,
                              LocalMatrixType(ansatzSpace.map().maxLocalSize(),
                                              testSpace.map().maxLocalSize(),
                                              RangeFieldType(0.0)));
    } // ensure enough tmp local matrices
    // apply local operator
    localOperator_.applyLocal(localAnsatzBaseFunctionSet,
                              localTestBaseFunctionSet,
                              intersection,
                              tmpLocalMatrices[0],
                              tmpLocalMatricesContainer[1]);
    // write local matrices to global (see below)
    addToMatrix(ansatzSpace, testSpace, entity, entity, tmpLocalMatrices[0], matrix);
  } // void assembleLocal() const

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  template< class AnsatzSpaceType,
            class TestSpaceType,
            class EntityType,
            class LocalMatrixType,
            class SystemMatrixType >
  void addToMatrix( const AnsatzSpaceType& ansatzSpace,
                    const TestSpaceType& testSpace,
                    const EntityType& ansatzEntity,
                    const EntityType& testEntity,
                    const LocalMatrixType& localMatrix,
                    SystemMatrixType& systemMatrix ) const
  {
    unsigned int rows = ansatzSpace.baseFunctionSet().local( ansatzEntity ).size();
    unsigned int cols = testSpace.baseFunctionSet().local( testEntity ).size();
    for( unsigned int i = 0; i < rows; ++i )
    {
      for( unsigned int j = 0; j < cols; ++j )
      {
        const unsigned int globalI = ansatzSpace.map().toGlobal( ansatzEntity, i );
        const unsigned int globalJ = testSpace.map().toGlobal( testEntity, j );

        systemMatrix.add(globalI, globalJ, localMatrix[i][j]);
      }
    }
  } // end method addToMatrix

  const LocalOperatorType& localOperator_;
}; // end class Boundary

} // end namespace Codim1

} // end namespace Local

} // end namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_MATRIX_HH
