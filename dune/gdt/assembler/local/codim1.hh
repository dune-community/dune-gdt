// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH

#include <vector>

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/walker/functors.hh>
#include <dune/stuff/common/tmp-storage.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/localfunctional/interface.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalAssembler {


template< class LocalOperatorImp >
class Codim1CouplingMatrix
{
  static_assert(std::is_base_of< LocalOperator::Codim1CouplingInterface< typename LocalOperatorImp::Traits >,
                                 LocalOperatorImp >::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1CouplingInterface!");
public:
  typedef LocalOperatorImp LocalOperatorType;

  explicit Codim1CouplingMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 4;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class TE, size_t TEd, size_t TEr, size_t TErC,
            class AE, size_t AEd, size_t AEr, size_t AErC,
            class TN, size_t TNd, size_t TNr, size_t TNrC,
            class AN, size_t ANd, size_t ANr, size_t ANrC,
            class IntersectionType, class MEE, class MNN, class MEN, class MNE, class R >
  void assembleLocal(const SpaceInterface< TE, TEd, TEr, TErC >& testSpaceEntity,
                     const SpaceInterface< AE, AEd, AEr, AErC >& ansatzSpaceEntity,
                     const SpaceInterface< TN, TNd, TNr, TNrC >& testSpaceNeighbor,
                     const SpaceInterface< AN, ANd, ANr, ANrC >& ansatzSpaceNeighbor,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< MEE, R >& entityEntityMatrix,
                     Dune::Stuff::LA::MatrixInterface< MNN, R >& neighborNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface< MEN, R >& entityNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface< MNE, R >& neighborEntityMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 4);
    // get and clear matrix
    auto& localEntityEntityMatrix = tmpLocalMatricesContainer[0][0];
    auto& localNeighborNeighborMatrix = tmpLocalMatricesContainer[0][1];
    auto& localEntityNeighborMatrix = tmpLocalMatricesContainer[0][2];
    auto& localNeighborEntityMatrix = tmpLocalMatricesContainer[0][3];
    localEntityEntityMatrix *= 0.0;
    localNeighborNeighborMatrix *= 0.0;
    localEntityNeighborMatrix *= 0.0;
    localNeighborEntityMatrix *= 0.0;
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // get entities
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    const auto neighborPtr = intersection.outside();
    const auto& neighbor = *neighborPtr;
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpaceEntity.base_function_set(entity), ansatzSpaceEntity.base_function_set(entity),
                         testSpaceNeighbor.base_function_set(neighbor), ansatzSpaceNeighbor.base_function_set(neighbor),
                         intersection,
                         localEntityEntityMatrix,
                         localNeighborNeighborMatrix,
                         localEntityNeighborMatrix,
                         localNeighborEntityMatrix,
                         tmpOperatorMatrices);
    // write local matrices to global
    const size_t rowsEn = testSpaceEntity.mapper().numDofs(entity);
    const size_t colsEn = ansatzSpaceEntity.mapper().numDofs(entity);
    const size_t rowsNe = testSpaceNeighbor.mapper().numDofs(neighbor);
    const size_t colsNe = ansatzSpaceNeighbor.mapper().numDofs(neighbor);
    auto& globalRowsEn = tmpIndicesContainer[0];
    auto& globalColsEn = tmpIndicesContainer[1];
    auto& globalRowsNe = tmpIndicesContainer[2];
    auto& globalColsNe = tmpIndicesContainer[3];
    assert(globalRowsEn.size() >= rowsEn);
    assert(globalColsEn.size() >= colsEn);
    assert(globalRowsNe.size() >= rowsNe);
    assert(globalColsNe.size() >= colsNe);
    testSpaceEntity.mapper().globalIndices(entity, globalRowsEn);
    ansatzSpaceEntity.mapper().globalIndices(entity, globalColsEn);
    testSpaceNeighbor.mapper().globalIndices(neighbor, globalRowsNe);
    ansatzSpaceNeighbor.mapper().globalIndices(neighbor, globalColsNe);
    assert(localEntityEntityMatrix.rows() >= rowsEn);
    assert(localEntityEntityMatrix.cols() >= colsEn);
    assert(localNeighborNeighborMatrix.rows() >= rowsNe);
    assert(localNeighborNeighborMatrix.cols() >= colsNe);
    assert(localEntityNeighborMatrix.rows() >= rowsEn);
    assert(localEntityNeighborMatrix.cols() >= colsNe);
    assert(localNeighborEntityMatrix.rows() >= rowsNe);
    assert(localNeighborEntityMatrix.cols() >= colsEn);
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      const auto& localEntityEntityMatrixRow = localEntityEntityMatrix[ii];
      const auto& localEntityNeighborMatrixRow = localEntityNeighborMatrix[ii];
      const size_t globalII = globalRowsEn[ii];
      for (size_t jj = 0; jj < colsEn; ++jj) {
        const size_t globalJJ = globalColsEn[jj];
        entityEntityMatrix.add_to_entry(globalII, globalJJ, localEntityEntityMatrixRow[jj]);
      }
      for (size_t jj = 0; jj < colsNe; ++jj) {
        const size_t globalJJ = globalColsNe[jj];
        entityNeighborMatrix.add_to_entry(globalII, globalJJ, localEntityNeighborMatrixRow[jj]);
      }
    }
    for (size_t ii = 0; ii < rowsNe; ++ii) {
      const auto& localNeighborEntityMatrixRow = localNeighborEntityMatrix[ii];
      const auto& localNeighborNeighborMatrixRow = localNeighborNeighborMatrix[ii];
      const size_t globalII = globalRowsNe[ii];
      for (size_t jj = 0; jj < colsEn; ++jj) {
        const size_t globalJJ = globalColsEn[jj];
        neighborEntityMatrix.add_to_entry(globalII, globalJJ, localNeighborEntityMatrixRow[jj]);
      }
      for (size_t jj = 0; jj < colsNe; ++jj) {
        const size_t globalJJ = globalColsNe[jj];
        neighborNeighborMatrix.add_to_entry(globalII, globalJJ, localNeighborNeighborMatrixRow[jj]);
      }
    }
  } // void assembleLocal(...) const

  template< class T, size_t Td, size_t Tr, size_t TrC,
            class A, size_t Ad, size_t Ar, size_t ArC,
            class IntersectionType, class M, class R >
  void assembleLocal(const SpaceInterface< T, Td, Tr, TrC >& testSpace,
                     const SpaceInterface< A, Ad, Ar, ArC >& ansatzSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< M, R >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    assembleLocal(testSpace, ansatzSpace, testSpace, ansatzSpace,
                  intersection,
                  systemMatrix, systemMatrix, systemMatrix, systemMatrix,
                  tmpLocalMatricesContainer,
                  tmpIndicesContainer);
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1CouplingMatrix


template< class LocalOperatorImp >
class Codim1BoundaryMatrix
{
  static_assert(std::is_base_of< LocalOperator::Codim1BoundaryInterface< typename LocalOperatorImp::Traits >,
                                 LocalOperatorImp >::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1BoundaryInterface!");
public:
  typedef LocalOperatorImp LocalOperatorType;

  explicit Codim1BoundaryMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class T, size_t Td, size_t Tr, size_t TrC,
            class A, size_t Ad, size_t Ar, size_t ArC,
            class IntersectionType, class M, class R >
  void assembleLocal(const SpaceInterface< T, Td, Tr, TrC >& testSpace,
                     const SpaceInterface< A, Ad, Ar, ArC >& ansatzSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< M, R >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix
    auto& localMatrix = tmpLocalMatricesContainer[0][0];
    localMatrix *= 0.0;
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpace.base_function_set(entity), ansatzSpace.base_function_set(entity),
                         intersection,
                         localMatrix, tmpOperatorMatrices);
    // write local matrices to global
    const size_t rows = testSpace.mapper().numDofs(entity);
    const size_t cols = ansatzSpace.mapper().numDofs(entity);
    auto& globalRows = tmpIndicesContainer[0];
    auto& globalCols = tmpIndicesContainer[1];
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    assert(localMatrix.size() >= rows);
    assert(localMatrix.size() >= cols);
    testSpace.mapper().globalIndices(entity, globalRows);
    ansatzSpace.mapper().globalIndices(entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localMatrixRow = localMatrix[ii];
      const size_t globalII = globalRows[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t globalJJ = globalCols[jj];
        systemMatrix.add_to_entry(globalII, globalJJ, localMatrixRow[jj]);
      }
    }
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1BoundaryMatrix


template< class LocalFunctionalImp >
class Codim1Vector
{
  static_assert(std::is_base_of< LocalFunctional::Codim1Interface< typename LocalFunctionalImp::Traits >,
                                 LocalFunctionalImp >::value,
                "LocalFunctionalImp has to be derived from LocalFunctional::Codim1Interface!");
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  explicit Codim1Vector(const LocalFunctionalType& fu)
    : localFunctional_(fu)
  {}

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired()};
  }

  template< class T, size_t d, size_t r, size_t rC, class IntersectionType, class V, class R >
  void assembleLocal(const SpaceInterface< T, d, r, rC >& testSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::VectorInterface< V, R >& systemVector,
                     std::vector< std::vector< Dune::DynamicVector< R > > >& tmpLocalVectorsContainer,
                     Dune::DynamicVector< size_t >& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalVectorsContainer.size() >= 2);
    assert(tmpLocalVectorsContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorsContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    auto& localVector = tmpLocalVectorsContainer[0][0];
    localVector *= 0.0;
    auto& tmpFunctionalVectors = tmpLocalVectorsContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    // apply local functional (results are in localVector)
    localFunctional_.apply(testSpace.base_function_set(entity), intersection, localVector, tmpFunctionalVectors);
    // write local vectors to global
    const size_t size = testSpace.mapper().numDofs(entity);
    assert(tmpIndicesContainer.size() >= size);
    assert(localVector.size() >= size);
    testSpace.mapper().globalIndices(entity, tmpIndicesContainer);
    for (size_t ii = 0; ii < size; ++ii) {
      const size_t globalII = tmpIndicesContainer[ii];
      systemVector.add_to_entry(globalII, localVector[ii]);
    }
  } // void assembleLocal(...) const

private:
  const LocalFunctionalType& localFunctional_;
}; // class Codim1Vector


template< class GridViewImp, class LocalOperatorType, class TestFunctionType, class AnsatzFunctionType, class FieldType >
class Codim1CouplingOperatorAccumulateFunctor
  : public Stuff::Grid::Functor::Codim1< GridViewImp >
{
  static_assert(std::is_base_of< LocalOperator::Codim1CouplingInterface< typename LocalOperatorType::Traits >,
                                 LocalOperatorType >::value,
                "LocalOperatorType has to be derived from LocalOperator::Codim1CouplingInterface!");
  static_assert(Stuff::is_localizable_function< TestFunctionType >::value,
                "TestFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< AnsatzFunctionType >::value,
                "AnsatzFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");

  typedef Stuff::Grid::Functor::Codim1< GridViewImp > BaseType;
  typedef DSC::TmpMatricesStorage< FieldType > TmpMatricesProviderType;
public:
  typedef typename BaseType::GridViewType     GridViewType;
  typedef typename BaseType::EntityType       EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  Codim1CouplingOperatorAccumulateFunctor(const GridViewType& grd_vw,
                                          const LocalOperatorType& local_op,
                                          const TestFunctionType& test_function,
                                          const AnsatzFunctionType& ansatz_function)
    : grid_view_(grd_vw)
    , local_operator_(local_op)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , result_(0)
    , finalized_(false)
  {
    // can not use make_unique here, at least clang does not get it
    tmp_storage_ = std::unique_ptr< TmpMatricesProviderType >(new TmpMatricesProviderType(
        {4, local_operator_.numTmpObjectsRequired()}, 1, 1));
  }

  virtual ~Codim1CouplingOperatorAccumulateFunctor() = default;

  FieldType compute_locally(const IntersectionType& intersection,
                            const EntityType& inside_entity,
                            const EntityType& outside_entity)
  {
#ifndef NDEBUG
    auto logger = DSC::TimedLogger().get("gdt.assembler.local.codim1couplingoperatoraccumulatefunctor");
#endif
    auto& tmp_storage = tmp_storage_->matrices();
    assert(tmp_storage.size() >= 2);
    assert(tmp_storage[0].size() >= 4);
    auto& local_operator_result_en_en = tmp_storage[0][0];
    auto& local_operator_result_ne_ne = tmp_storage[0][1];
    auto& local_operator_result_en_ne = tmp_storage[0][2];
    auto& local_operator_result_ne_en = tmp_storage[0][3];
    auto& tmp_matrices                = tmp_storage[1];
    // get the local functions
    const auto local_test_function_en    = test_function_.local_function(inside_entity);
    const auto local_test_function_ne    = test_function_.local_function(outside_entity);
    const auto local_ansatz_function_en = ansatz_function_.local_function(inside_entity);
    const auto local_ansatz_function_ne = ansatz_function_.local_function(outside_entity);
    // apply the local operator
    this->local_operator_.apply(*local_test_function_en,
                                *local_ansatz_function_en,
                                *local_ansatz_function_ne,
                                *local_test_function_ne,
                                intersection,
                                local_operator_result_en_en,
                                local_operator_result_ne_ne,
                                local_operator_result_en_ne,
                                local_operator_result_ne_en,
                                tmp_matrices);
    assert(local_operator_result_en_en.rows() >= 1);
    assert(local_operator_result_en_en.cols() >= 1);
    assert(local_operator_result_ne_ne.rows() >= 1);
    assert(local_operator_result_ne_ne.cols() >= 1);
    assert(local_operator_result_en_ne.rows() >= 1);
    assert(local_operator_result_en_ne.cols() >= 1);
    assert(local_operator_result_ne_en.rows() >= 1);
    assert(local_operator_result_ne_en.cols() >= 1);
#ifndef NDEBUG
    logger.debug() << intersection.geometry().center() << ": "
                   << local_operator_result_en_en[0][0]
                    + local_operator_result_ne_ne[0][0]
                    + local_operator_result_en_ne[0][0]
                    + local_operator_result_ne_en[0][0] << std::endl;
#endif
    return local_operator_result_en_en[0][0]
         + local_operator_result_ne_ne[0][0]
         + local_operator_result_en_ne[0][0]
         + local_operator_result_ne_en[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override
  {
    *result_ += compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  FieldType result() const
  {
    if (!finalized_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong, "Call finalize() first!");
    return finalized_result_;
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  DS::PerThreadValue< FieldType > result_;
  std::unique_ptr< TmpMatricesProviderType > tmp_storage_;
  bool finalized_;
  FieldType finalized_result_;
}; // class Codim1CouplingOperatorAccumulateFunctor


template< class GridViewImp, class LocalOperatorType, class TestFunctionType, class AnsatzFunctionType, class FieldType >
class Codim1BoundaryOperatorAccumulateFunctor
  : public Stuff::Grid::Functor::Codim1< GridViewImp >
{
  static_assert(std::is_base_of< LocalOperator::Codim1BoundaryInterface< typename LocalOperatorType::Traits >,
                                 LocalOperatorType >::value,
                "LocalOperatorType has to be derived from LocalOperator::Codim1BoundaryInterface!");
  static_assert(Stuff::is_localizable_function< TestFunctionType >::value,
                "TestFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< AnsatzFunctionType >::value,
                "AnsatzFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");

  typedef Stuff::Grid::Functor::Codim1< GridViewImp > BaseType;
  typedef DSC::TmpMatricesStorage< FieldType > TmpMatricesProviderType;
public:
  typedef typename BaseType::GridViewType     GridViewType;
  typedef typename BaseType::EntityType       EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  Codim1BoundaryOperatorAccumulateFunctor(const GridViewType& grd_vw,
                                          const LocalOperatorType& local_op,
                                          const TestFunctionType& test_function,
                                          const AnsatzFunctionType& ansatz_function)
    : grid_view_(grd_vw)
    , local_operator_(local_op)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , result_(0)
    , finalized_(false)
  {
    // can not use make_unique here, at least clang does not get it
    tmp_storage_ = std::unique_ptr< TmpMatricesProviderType >(new TmpMatricesProviderType(
        {1, local_operator_.numTmpObjectsRequired()}, 1, 1));
  }

  virtual ~Codim1BoundaryOperatorAccumulateFunctor() = default;

  FieldType compute_locally(const IntersectionType& intersection,
                            const EntityType& inside_entity,
                            const EntityType& /*outside_entity*/)
  {
#ifndef NDEBUG
    auto logger = DSC::TimedLogger().get("gdt.assembler.local.codim1boundaryoperatoraccumulatefunctor");
#endif
    assert(tmp_storage_->matrices().size() >= 2);
    assert(tmp_storage_->matrices()[0].size() >= 1);
    auto& local_operator_result = tmp_storage_->matrices()[0][0];
    auto& tmp_matrices          = tmp_storage_->matrices()[1];
    // get the local functions
    const auto local_test_function    = test_function_.local_function(inside_entity);
    const auto local_ansatz_function = ansatz_function_.local_function(inside_entity);
    // apply the local operator
    this->local_operator_.apply(*local_test_function,
                                *local_ansatz_function,
                                intersection,
                                local_operator_result,
                                tmp_matrices);
    assert(local_operator_result.rows() >= 1);
    assert(local_operator_result.cols() >= 1);
#ifndef NDEBUG
    logger.debug() << intersection.geometry().center() << ": " << local_operator_result[0][0] << std::endl;
#endif
    return local_operator_result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override
  {
    *result_ += compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  FieldType result() const
  {
    if (!finalized_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong, "Call finalize() first!");
    return finalized_result_;
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  DS::PerThreadValue< FieldType > result_;
  std::unique_ptr< TmpMatricesProviderType > tmp_storage_;
  bool finalized_;
  FieldType finalized_result_;
}; // class Codim1BoundaryOperatorAccumulateFunctor


template< class LocalOperatorImp >
class Codim1CouplingFV
{
  static_assert(std::is_base_of< LocalOperator::Codim1CouplingInterface< typename LocalOperatorImp::Traits >,
                                 LocalOperatorImp >::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1CouplingInterface!");
public:
  typedef LocalOperatorImp LocalOperatorType;

  explicit Codim1CouplingFV(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class SourceSpaceType, class RangeSpaceType, class VectorType, class IntersectionType, class RangeFieldType >
  void assembleLocal(const Dune::GDT::DiscreteFunction< SourceSpaceType, VectorType >& discreteFunction,
                     Dune::GDT::DiscreteFunction< RangeSpaceType, VectorType >& discreteFunctionUpdate,
                     const IntersectionType& intersection,
                     Dune::DynamicMatrix< RangeFieldType >& updateMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< RangeFieldType > > >& tmpLocalMatrices) const
  {
    // check
    const size_t dimRange = discreteFunction.space().dimRange;
    assert(intersection.neighbor());
    assert(updateMatrix.cols() >= 1);
    assert(updateMatrix.rows() >= dimRange);
    //assert(discreteFunction.vector().size() == discreteFunctionUpdate.vector().size());
    //clear matrix
    updateMatrix *= 0.0;
    //get entity and neighbor and local discrete functions
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    const auto neighborPtr = intersection.outside();
    const auto& neighbor = *neighborPtr;
    const auto entityAverage = discreteFunction.local_discrete_function(entity);
    const auto neighborAverage = discreteFunction.local_discrete_function(neighbor);
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(*entityAverage, *entityAverage,
                         *neighborAverage, *neighborAverage,
                         intersection,
                         updateMatrix,
                         updateMatrix,
                         updateMatrix,
                         updateMatrix,
                         tmpLocalMatrices[0]);
    // write value from updateMatrix to discreteFunctionUpdate
    for (size_t kk = 0; kk < dimRange; ++kk)
      discreteFunctionUpdate.local_discrete_function(entity)->vector().add(kk, updateMatrix[kk][0]);
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1CouplingFV

template< class LocalOperatorImp >
class Codim1BoundaryFV
{
  static_assert(std::is_base_of< LocalOperator::Codim1BoundaryInterface< typename LocalOperatorImp::Traits >,
                                 LocalOperatorImp >::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1CouplingInterface!");
public:
  typedef LocalOperatorImp LocalOperatorType;

  explicit Codim1BoundaryFV(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class SourceSpaceType, class RangeSpaceType, class VectorType, class IntersectionType, class RangeFieldType >
  void assembleLocal(const Dune::GDT::DiscreteFunction< SourceSpaceType, VectorType >& discreteFunction,
                     Dune::GDT::DiscreteFunction< RangeSpaceType, VectorType >& discreteFunctionUpdate,
                     const IntersectionType& intersection,
                     Dune::DynamicMatrix< RangeFieldType >& updateMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< RangeFieldType > > >& tmpLocalMatrices) const
  {
    // check
    const size_t dimRange = discreteFunction.space().dimRange;
    assert(updateMatrix.cols() >= 1);
    assert(updateMatrix.rows() >= dimRange);
//    assert(discreteFunction.vector().size() == discreteFunctionUpdate.vector().size());
    //clear matrix
    updateMatrix *= 0.0;
    //get entity and neighbor and local discrete functions
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    const auto entityAverage = discreteFunction.local_discrete_function(entity);
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(*entityAverage, *entityAverage, intersection, updateMatrix, tmpLocalMatrices[0]);
    // write value from updateMatrix to discreteFunctionUpdate
    for (size_t kk = 0; kk < dimRange; ++kk)
      discreteFunctionUpdate.local_discrete_function(entity)->vector().add(kk, updateMatrix[kk][0]);
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1BoundaryFV

} // namespace LocalAssembler
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
