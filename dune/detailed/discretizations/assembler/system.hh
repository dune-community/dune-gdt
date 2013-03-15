#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH

#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Assembler {

template< class TestFunctionSpaceImp, class AnsatzFunctionSpaceImp = TestFunctionSpaceImp >
class System
{
public:
  typedef TestFunctionSpaceImp TestFunctionSpaceType;

  typedef AnsatzFunctionSpaceImp AnsatzFunctionSpaceType;

  typedef System< TestFunctionSpaceImp, AnsatzFunctionSpaceImp > ThisType;

private:
  typedef typename TestFunctionSpaceType::GridViewType GridViewType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;

  typedef typename TestFunctionSpaceType::FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;

  typedef Dune::DynamicVector< RangeFieldType > LocalVectorType;

  typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;

  typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;

  class LocalMatrixAssemblerApplication
  {
  public:
    virtual ~LocalMatrixAssemblerApplication(){}

    virtual void apply(const TestFunctionSpaceType& /*_testSpace*/,
                       const AnsatzFunctionSpaceType& /*_ansatzSpace*/,
                       const EntityType& /*_entity*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/) const = 0;

    virtual std::vector< unsigned int > numTmpObjectsRequired() const = 0;
  }; // class LocalOperatorApplication

  template< class LocalMatrixAssemblerType, class MatrixType >
  class LocalMatrixAssemblerApplicationWrapper
    : public LocalMatrixAssemblerApplication
  {
  public:
    LocalMatrixAssemblerApplicationWrapper(const std::shared_ptr< const LocalMatrixAssemblerType > _localMatrixAssembler,
                                           std::shared_ptr< MatrixType > _matrix)
      : localMatrixAssembler_(_localMatrixAssembler)
      , matrix_(_matrix)
    {}

    virtual void apply(const TestFunctionSpaceType& _testSpace,
                       const AnsatzFunctionSpaceType& _ansatzSpace,
                       const EntityType& _entity,
                       LocalMatricesContainerType& _localMatricesContainer) const
    {
      localMatrixAssembler_->assembleLocal(_testSpace, _ansatzSpace, _entity, *matrix_, _localMatricesContainer);
    } // virtual void applyLocal(...) const

    virtual std::vector< unsigned int > numTmpObjectsRequired() const
    {
      return localMatrixAssembler_->numTmpObjectsRequired();
    } // virtual std::vector< unsigned int > numTmpObjectsRequired() const

  private:
    const std::shared_ptr< const LocalMatrixAssemblerType > localMatrixAssembler_;
    std::shared_ptr< MatrixType > matrix_;
  }; // class LocalMatrixAssemblerApplicationWrapper

  class LocalVectorAssemblerApplication
  {
  public:
    virtual ~LocalVectorAssemblerApplication(){}

    virtual void apply(const TestFunctionSpaceType& /*_testSpace*/,
                       const EntityType& /*_entity*/,
                       LocalVectorsContainerType& /*_localMatricesContainer*/) const = 0;

    virtual std::vector< unsigned int > numTmpObjectsRequired() const = 0;
  }; // class LocalOperatorApplication

  template< class LocalVectorAssemblerType, class VectorType >
  class LocalVectorAssemblerApplicationWrapper
    : public LocalVectorAssemblerApplication
  {
  public:
    LocalVectorAssemblerApplicationWrapper(const std::shared_ptr< const LocalVectorAssemblerType > _localVectorAssembler,
                                           std::shared_ptr< VectorType > _vector)
      : localVectorAssembler_(_localVectorAssembler)
      , vector_(_vector)
    {}

    virtual void apply(const TestFunctionSpaceType& _testSpace,
                       const EntityType& _entity,
                       LocalVectorsContainerType& _localVectorsContainer) const
    {
      localVectorAssembler_->assembleLocal(_testSpace, _entity, *vector_, _localVectorsContainer);
    } // virtual void applyLocal(...) const

    virtual std::vector< unsigned int > numTmpObjectsRequired() const
    {
      return localVectorAssembler_->numTmpObjectsRequired();
    } // virtual std::vector< unsigned int > numTmpObjectsRequired() const

  private:
    const std::shared_ptr< const LocalVectorAssemblerType > localVectorAssembler_;
    std::shared_ptr< VectorType > vector_;
  }; // class LocalMatrixAssemblerApplicationWrapper

public:
  System(const TestFunctionSpaceType& _testSpace, const AnsatzFunctionSpaceType& _ansatzSpace)
    :  testSpace_(_testSpace)
    ,  ansatzSpace_(_ansatzSpace)
  {}

  System(const TestFunctionSpaceType& _testSpace)
    :  testSpace_(_testSpace)
    ,  ansatzSpace_(_testSpace)
  {}

  ~System()
  {
    for (auto& localMatrixAssembler: localMatrixAssemblers_)
      delete localMatrixAssembler;
    for (auto& localVectorAssembler: localVectorAssemblers_)
      delete localVectorAssembler;
  }

  const TestFunctionSpaceType& testSpace()
  {
    return testSpace_;
  }

  const AnsatzFunctionSpaceType& ansatzSpace()
  {
    return ansatzSpace_;
  }

  template< class LocalMatrixAssemblerType, class MatrixType >
  void addLocalMatrixAssembler(const std::shared_ptr< const LocalMatrixAssemblerType > _localMatrixAssembler,
                               std::shared_ptr< MatrixType > _matrix)
  {
    typedef LocalMatrixAssemblerApplicationWrapper< LocalMatrixAssemblerType, MatrixType > WrapperType;
    WrapperType* wrapper = new WrapperType(_localMatrixAssembler, _matrix);
    localMatrixAssemblers_.push_back(wrapper);
  }

  template< class LocalVectorAssemblerType, class VectorType >
  void addLocalVectorAssembler(const std::shared_ptr< const LocalVectorAssemblerType > _localVectorAssembler,
                               std::shared_ptr< VectorType > _vector)
  {
    typedef LocalVectorAssemblerApplicationWrapper< LocalVectorAssemblerType, VectorType > WrapperType;
    WrapperType* wrapper = new WrapperType(_localVectorAssembler, _vector);
    localVectorAssemblers_.push_back(wrapper);
  }

  void assemble() const
  {
    // common tmp storage for all entities
    // * for the matrix assemblers
    std::vector< unsigned int > numberOfTmpMatricesNeeded(2, 0);
    for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii) {
      const std::vector< unsigned int > tmp = localMatrixAssemblers_[ii]->numTmpObjectsRequired();
      numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
      numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
    }
    std::vector< LocalMatrixType > tmpLocalAssemblerMatrices( numberOfTmpMatricesNeeded[0],
                                                              LocalMatrixType(testSpace_.map().maxLocalSize(),
                                                                              ansatzSpace_.map().maxLocalSize(),
                                                                              RangeFieldType(0)));
    std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                            LocalMatrixType(testSpace_.map().maxLocalSize(),
                                                                            ansatzSpace_.map().maxLocalSize(),
                                                                            RangeFieldType(0)));
    std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // * for the vector assemblers
    std::vector< unsigned int > numberOfTmpVectorsNeeded(2, 0);
    for (unsigned int ii = 0; ii < localVectorAssemblers_.size(); ++ii) {
      const std::vector< unsigned int > tmp = localVectorAssemblers_[ii]->numTmpObjectsRequired();
      numberOfTmpVectorsNeeded[0] = std::max(numberOfTmpVectorsNeeded[0], tmp[0]);
      numberOfTmpVectorsNeeded[1] = std::max(numberOfTmpVectorsNeeded[1], tmp[1]);
    }
    std::vector< LocalVectorType > tmpLocalAssemblerVectors(numberOfTmpVectorsNeeded[0],
                                                            LocalVectorType(testSpace_.map().maxLocalSize(),
                                                                            RangeFieldType(0)));
    std::vector< LocalVectorType > tmpLocalFunctionalVectors( numberOfTmpVectorsNeeded[1],
                                                              LocalVectorType(testSpace_.map().maxLocalSize(),
                                                                              RangeFieldType(0)));
    std::vector< std::vector< LocalVectorType > > tmpLocalVectorsContainer;
    tmpLocalVectorsContainer.push_back(tmpLocalAssemblerVectors);
    tmpLocalVectorsContainer.push_back(tmpLocalFunctionalVectors);

    // walk the grid
    typedef typename GridViewType::template Codim< 0 >::Iterator EntityIteratorType;
    for(EntityIteratorType entityIt = ansatzSpace_.gridView().template begin< 0 >();
        entityIt != ansatzSpace_.gridView().template end< 0 >();
        ++entityIt ) {
      const EntityType& entity = *entityIt;
      // assemble local matrices
      for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii)
        localMatrixAssemblers_[ii]->apply(testSpace_, ansatzSpace_, entity, tmpLocalMatricesContainer);
      // assemble local vectors
      for (unsigned int ii = 0; ii < localVectorAssemblers_.size(); ++ii)
        localVectorAssemblers_[ii]->apply(testSpace_, entity, tmpLocalVectorsContainer);
    } // walk the grid
  } // void assemble() const

  template< class MatrixType, class VectorType >
  void applyConstraints(MatrixType& matrix, VectorType& vector) const
  {
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType EntityIteratorType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    // walk the grid to apply constraints
    const ConstraintsType& constraints = testSpace_.constraints();
    for(EntityIteratorType entityIterator = testSpace_.gridPart().template begin< 0 >();
        entityIterator != testSpace_.gridPart().template end< 0 >();
        ++entityIterator ) {
      const EntityType& entity = *entityIterator;
      const LocalConstraintsType& localConstraints = constraints.local(entity);
      applyLocalMatrixConstraints(localConstraints, matrix);
      applyLocalVectorConstraints(localConstraints, vector);
    } // walk the grid to apply constraints
  } // void applyConstraints(MatrixType& matrix, VectorType& vector) const

  template< class MatrixType >
  void applyMatrixConstraints(MatrixType& matrix) const
  {
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType EntityIteratorType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    // walk the grid to apply constraints
    const ConstraintsType& constraints = testSpace_.constraints();
    for(EntityIteratorType entityIterator = testSpace_.gridPart().template begin< 0 >();
        entityIterator != testSpace_.gridPart().template end< 0 >();
        ++entityIterator ) {
      const EntityType& entity = *entityIterator;
      const LocalConstraintsType& localConstraints = constraints.local(entity);
      applyLocalMatrixConstraints(localConstraints, matrix);
    } // walk the grid to apply constraints
  } // void applyMatrixConstraints(MatrixType& matrix) const

  template< class VectorType >
  void applyVectorConstraints(VectorType& vector) const
  {
    typedef typename AnsatzFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType EntityIteratorType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename AnsatzFunctionSpaceType::ConstraintsType ConstraintsType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    // walk the grid to apply constraints
    const ConstraintsType& constraints = testSpace_.constraints();
    for(EntityIteratorType entityIterator = testSpace_.gridPart().template begin< 0 >();
        entityIterator != testSpace_.gridPart().template end< 0 >();
        ++entityIterator ) {
      const EntityType& entity = *entityIterator;
      const LocalConstraintsType& localConstraints = constraints.local(entity);
      applyLocalVectorConstraints(localConstraints, vector);
    } // walk the grid to apply constraints
  } // void applyVectorConstraints(VectorType& vector) const

private:
  System(const ThisType&);
  ThisType& operator=( const ThisType& );

  template< class LocalConstraintsType, class MatrixType >
  void applyLocalMatrixConstraints(const LocalConstraintsType& localConstraints, MatrixType& matrix) const
  {
    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
      const unsigned int rowDof = localConstraints.rowDofs(i);
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); ++j) {
        matrix.set(rowDof, localConstraints.columnDofs(j), localConstraints.localMatrix(i,j));
      }
    }
  } // void applyLocalMatrixConstraints(...)

  template< class LocalConstraintsType, class VectorType >
  void applyLocalVectorConstraints( const LocalConstraintsType& localConstraints, VectorType& vector ) const
  {
    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i )
    {
        vector.set(localConstraints.rowDofs(i), 0.0);
    }
  } // void applyLocalVectorConstraints(...)

  const TestFunctionSpaceType& testSpace_;
  const AnsatzFunctionSpaceType& ansatzSpace_;
  std::vector< LocalMatrixAssemblerApplication* > localMatrixAssemblers_;
  std::vector< LocalVectorAssemblerApplication* > localVectorAssemblers_;
}; // class System

} // namespace Assembler
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH
