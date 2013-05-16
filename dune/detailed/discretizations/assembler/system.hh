#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH

#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/la/container/interface.hh>
#include <dune/stuff/la/container/pattern.hh>

#include <dune/detailed/discretizations/space/interface.hh>
#include <dune/detailed/discretizations/space/constraints.hh>

#include "local/codim0.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {

template< class TestSpaceImp, class AnsatzSpaceImp >
class SystemAssembler
{
public:
  typedef SpaceInterface< typename TestSpaceImp::Traits >   TestSpaceType;
  typedef SpaceInterface< typename AnsatzSpaceImp::Traits > AnsatzSpaceType;

private:
  typedef typename TestSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
  typedef Dune::DynamicVector< RangeFieldType > LocalVectorType;
  typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;
  typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;
  typedef std::vector< Dune::DynamicVector< size_t > > IndicesContainer;

  class LocalCodim0MatrixAssemblerApplication
  {
  public:
    virtual ~LocalCodim0MatrixAssemblerApplication(){}

    virtual void apply(const TestSpaceType& /*_testSpace*/,
                       const AnsatzSpaceType& /*_ansatzSpace*/,
                       const EntityType& /*_entity*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/,
                       IndicesContainer& /*indicesContainer*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class L, class M >
  class LocalCodim0MatrixAssemblerWrapper
    : public LocalCodim0MatrixAssemblerApplication
  {
  public:
    LocalCodim0MatrixAssemblerWrapper(const LocalAssemblerCodim0Matrix< L >& localAssembler,
                                      Dune::Stuff::LA::Container::MatrixInterface< M >& matrix)
      : localMatrixAssembler_(localAssembler)
      , matrix_(matrix)
    {}

    virtual void apply(const TestSpaceType& testSpace,
                       const AnsatzSpaceType& ansatzSpace,
                       const EntityType& entity,
                       LocalMatricesContainerType& localMatricesContainer,
                       IndicesContainer& indicesContainer) const
    {
      localMatrixAssembler_.assembleLocal(testSpace, ansatzSpace, entity, matrix_, localMatricesContainer, indicesContainer);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const
    {
      return localMatrixAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssemblerCodim0Matrix< L >& localMatrixAssembler_;
    Dune::Stuff::LA::Container::MatrixInterface< M >& matrix_;
  }; // class LocalCodim0MatrixAssemblerWrapper

  class LocalCodim0VectorAssemblerApplication
  {
  public:
    virtual ~LocalCodim0VectorAssemblerApplication(){}

    virtual void apply(const TestSpaceType& /*_testSpace*/,
                       const EntityType& /*_entity*/,
                       LocalVectorsContainerType& /*_localVectorsContainer*/,
                       Dune::DynamicVector< size_t >& /*indices*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class L, class V >
  class LocalCodim0VectorAssemblerWrapper
    : public LocalCodim0VectorAssemblerApplication
  {
  public:
    LocalCodim0VectorAssemblerWrapper(const LocalAssemblerCodim0Vector< L >& localAssembler,
                                      Dune::Stuff::LA::Container::VectorInterface< V >& vector)
      : localVectorAssembler_(localAssembler)
      , vector_(vector)
    {}

    virtual void apply(const TestSpaceType& testSpace,
                       const EntityType& entity,
                       LocalVectorsContainerType& localVectorsContainer,
                       Dune::DynamicVector< size_t >& indices) const
    {
      localVectorAssembler_.assembleLocal(testSpace, entity, vector_, localVectorsContainer, indices);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const
    {
      return localVectorAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssemblerCodim0Vector< L >& localVectorAssembler_;
    Dune::Stuff::LA::Container::VectorInterface< V >& vector_;
  }; // class LocalCodim0VectorAssemblerWrapper

public:
  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz)
    :  testSpace_(test)
    ,  ansatzSpace_(ansatz)
  {}

  SystemAssembler(const TestSpaceType& test)
    :  testSpace_(test)
    ,  ansatzSpace_(test)
  {}

  ~SystemAssembler()
  {
    for (auto& localCodim0MatrixAssembler: localCodim0MatrixAssemblers_)
      delete localCodim0MatrixAssembler;
    for (auto& localCodim0VectorAssembler: localCodim0VectorAssemblers_)
      delete localCodim0VectorAssembler;
  }

  const TestSpaceType& testSpace()
  {
    return testSpace_;
  }

  const AnsatzSpaceType& ansatzSpace()
  {
    return ansatzSpace_;
  }

  template< class L, class M >
  void addLocalMatrixAssembler(const LocalAssemblerCodim0Matrix< L >& localAssembler,
                               Dune::Stuff::LA::Container::MatrixInterface< M >& matrix)
  {
    assert(matrix.rows() == testSpace_.mapper().size());
    assert(matrix.cols() == ansatzSpace_.mapper().size());
    localCodim0MatrixAssemblers_.push_back(new LocalCodim0MatrixAssemblerWrapper< L, M >(localAssembler, matrix));
  }

  template< class L, class V >
  void addLocalVectorAssembler(const LocalAssemblerCodim0Vector< L >& localAssembler,
                               Dune::Stuff::LA::Container::VectorInterface< V >& vector)
  {
    assert(vector.size() == int(testSpace_.mapper().size()));
    localCodim0VectorAssemblers_.push_back(new LocalCodim0VectorAssemblerWrapper< L, V >(localAssembler, vector));
  }

  void assemble() const
  {
    // common tmp storage for all entities
    // * for the matrix assemblers
    std::vector< size_t > numberOfTmpMatricesNeeded(2, 0);
    for (auto& localCodim0MatrixAssembler : localCodim0MatrixAssemblers_) {
      const auto tmp = localCodim0MatrixAssembler->numTmpObjectsRequired();
      assert(tmp.size() == 2);
      numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
      numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
    }
    const size_t maxLocalSize = std::max(testSpace_.mapper().maxNumDofs(), ansatzSpace_.mapper().maxNumDofs());
    std::vector< LocalMatrixType > tmpLocalAssemblerMatrices( numberOfTmpMatricesNeeded[0],
                                                              LocalMatrixType(maxLocalSize,
                                                                              maxLocalSize,
                                                                              RangeFieldType(0)));
    std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                            LocalMatrixType(maxLocalSize,
                                                                            maxLocalSize,
                                                                            RangeFieldType(0)));
    std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // * for the vector assemblers
    std::vector< size_t > numberOfTmpVectorsNeeded(2, 0);
    for (auto& localCodim0VectorAssembler : localCodim0VectorAssemblers_) {
      const auto tmp = localCodim0VectorAssembler->numTmpObjectsRequired();
      assert(tmp.size() == 2);
      numberOfTmpVectorsNeeded[0] = std::max(numberOfTmpVectorsNeeded[0], tmp[0]);
      numberOfTmpVectorsNeeded[1] = std::max(numberOfTmpVectorsNeeded[1], tmp[1]);
    }
    std::vector< LocalVectorType > tmpLocalAssemblerVectors(numberOfTmpVectorsNeeded[0],
                                                            LocalVectorType(maxLocalSize,
                                                                            RangeFieldType(0)));
    std::vector< LocalVectorType > tmpLocalFunctionalVectors( numberOfTmpVectorsNeeded[1],
                                                              LocalVectorType(maxLocalSize,
                                                                              RangeFieldType(0)));
    std::vector< std::vector< LocalVectorType > > tmpLocalVectorsContainer;
    tmpLocalVectorsContainer.push_back(tmpLocalAssemblerVectors);
    tmpLocalVectorsContainer.push_back(tmpLocalFunctionalVectors);
    // * for the global indices
    std::vector< Dune::DynamicVector< size_t > > tmpIndices = {
        Dune::DynamicVector< size_t >(maxLocalSize)
      , Dune::DynamicVector< size_t >(maxLocalSize)
    };

    // walk the grid
    const auto entityEndIt = testSpace_.gridPart().template end< 0 >();
    for(auto entityIt = testSpace_.gridPart().template begin< 0 >(); entityIt != entityEndIt; ++entityIt ) {
      const EntityType& entity = *entityIt;
      // assemble local matrices
      for (auto& localCodim0MatrixAssembler : localCodim0MatrixAssemblers_) {
        localCodim0MatrixAssembler->apply(testSpace_, ansatzSpace_, entity, tmpLocalMatricesContainer, tmpIndices);
      }
      // assemble local vectors
      for (auto& localCodim0VectorAssembler : localCodim0VectorAssemblers_) {
        localCodim0VectorAssembler->apply(testSpace_, entity, tmpLocalVectorsContainer, tmpIndices[0]);
      }
    } // walk the grid
  } // void assemble() const

  template< class ConstraintsType, class M, class V >
  void applyConstraints(ConstraintsType& constraints,
                        Dune::Stuff::LA::Container::MatrixInterface< M >& matrix,
                        Dune::Stuff::LA::Container::VectorInterface< V >& vector) const
  {
    // walk the grid
    const auto entityEndIt = testSpace_.gridPart().template end< 0 >();
    for (auto entityIt = testSpace_.gridPart().template begin< 0 >(); entityIt != entityEndIt; ++entityIt) {
      const EntityType& entity = *entityIt;
      testSpace_.localConstraints(entity, constraints);
      applyLocalMatrixConstraints(constraints, matrix);
      applyLocalVectorConstraints(constraints, vector);
    } // walk the grid to apply constraints
  } // void applyConstraints(...) const

private:
  template< class M >
  void applyLocalMatrixConstraints(const Constraints::LocalDefault< RangeFieldType >& localConstraints,
                                   Dune::Stuff::LA::Container::MatrixInterface< M >& matrix) const
  {
    for (size_t ii = 0; ii < localConstraints.rows(); ++ii) {
      const size_t row = localConstraints.globalRow(ii);
      for (size_t jj = 0; jj < localConstraints.cols(); ++jj) {
        matrix.set(row, localConstraints.globalCol(jj), localConstraints.value(ii, jj));
      }
    }
  } // void applyLocalMatrixConstraints(...)

  template< class V >
  void applyLocalVectorConstraints(const Constraints::LocalDefault< RangeFieldType >& localConstraints,
                                   Dune::Stuff::LA::Container::VectorInterface< V >& vector ) const
  {
    for (size_t ii = 0; ii < localConstraints.rows(); ++ii) {
      vector.set(localConstraints.globalRow(ii), RangeFieldType(0));
    }
  } // void applyLocalVectorConstraints(...)

  const TestSpaceType& testSpace_;
  const AnsatzSpaceType& ansatzSpace_;
  std::vector< LocalCodim0MatrixAssemblerApplication* > localCodim0MatrixAssemblers_;
  std::vector< LocalCodim0VectorAssemblerApplication* > localCodim0VectorAssemblers_;
}; // class SystemAssembler

} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_SYSTEM_HH
