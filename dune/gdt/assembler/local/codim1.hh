// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH

#include <vector>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/walker/functors.hh>

#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/localfunctional/interface.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalAssembler {


template <class LocalOperatorImp>
class Codim1CouplingMatrix
{
  static_assert(std::is_base_of<LocalOperator::Codim1CouplingInterface<typename LocalOperatorImp::Traits>,
                                LocalOperatorImp>::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1CouplingInterface!");

public:
  typedef LocalOperatorImp LocalOperatorType;

  Codim1CouplingMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {
  }

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 4;

public:
  std::vector<size_t> numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template <class TE, class AE, class TN, class AN, class IntersectionType, class MEE, class MNN, class MEN, class MNE,
            class R>
  void assembleLocal(const SpaceInterface<TE>& testSpaceEntity, const SpaceInterface<AE>& ansatzSpaceEntity,
                     const SpaceInterface<TN>& testSpaceNeighbor, const SpaceInterface<AN>& ansatzSpaceNeighbor,
                     const IntersectionType& intersection, Dune::Stuff::LA::MatrixInterface<MEE>& entityEntityMatrix,
                     Dune::Stuff::LA::MatrixInterface<MNN>& neighborNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface<MEN>& entityNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface<MNE>& neighborEntityMatrix,
                     std::vector<std::vector<Dune::DynamicMatrix<R>>>& tmpLocalMatricesContainer,
                     std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 4);
    // get and clear matrix
    Dune::DynamicMatrix<R>& localEntityEntityMatrix     = tmpLocalMatricesContainer[0][0];
    Dune::DynamicMatrix<R>& localNeighborNeighborMatrix = tmpLocalMatricesContainer[0][1];
    Dune::DynamicMatrix<R>& localEntityNeighborMatrix   = tmpLocalMatricesContainer[0][2];
    Dune::DynamicMatrix<R>& localNeighborEntityMatrix   = tmpLocalMatricesContainer[0][3];
    localEntityEntityMatrix *= 0.0;
    localNeighborNeighborMatrix *= 0.0;
    localEntityNeighborMatrix *= 0.0;
    localNeighborEntityMatrix *= 0.0;
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // get entities
    const auto entityPtr   = intersection.inside();
    const auto& entity     = *entityPtr;
    const auto neighborPtr = intersection.outside();
    const auto& neighbor   = *neighborPtr;
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpaceEntity.base_function_set(entity),
                         ansatzSpaceEntity.base_function_set(entity),
                         testSpaceNeighbor.base_function_set(neighbor),
                         ansatzSpaceNeighbor.base_function_set(neighbor),
                         intersection,
                         localEntityEntityMatrix,
                         localNeighborNeighborMatrix,
                         localEntityNeighborMatrix,
                         localNeighborEntityMatrix,
                         tmpOperatorMatrices);
    // write local matrices to global
    const size_t rowsEn                       = testSpaceEntity.mapper().numDofs(entity);
    const size_t colsEn                       = ansatzSpaceEntity.mapper().numDofs(entity);
    const size_t rowsNe                       = testSpaceNeighbor.mapper().numDofs(neighbor);
    const size_t colsNe                       = ansatzSpaceNeighbor.mapper().numDofs(neighbor);
    Dune::DynamicVector<size_t>& globalRowsEn = tmpIndicesContainer[0];
    Dune::DynamicVector<size_t>& globalColsEn = tmpIndicesContainer[1];
    Dune::DynamicVector<size_t>& globalRowsNe = tmpIndicesContainer[2];
    Dune::DynamicVector<size_t>& globalColsNe = tmpIndicesContainer[3];
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
      const auto& localEntityEntityMatrixRow   = localEntityEntityMatrix[ii];
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
      const auto& localNeighborEntityMatrixRow   = localNeighborEntityMatrix[ii];
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

  template <class T, class A, class IntersectionType, class M, class R>
  void assembleLocal(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace,
                     const IntersectionType& intersection, Dune::Stuff::LA::MatrixInterface<M>& systemMatrix,
                     std::vector<std::vector<Dune::DynamicMatrix<R>>>& tmpLocalMatricesContainer,
                     std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const
  {
    assembleLocal(testSpace,
                  ansatzSpace,
                  testSpace,
                  ansatzSpace,
                  intersection,
                  systemMatrix,
                  systemMatrix,
                  systemMatrix,
                  systemMatrix,
                  tmpLocalMatricesContainer,
                  tmpIndicesContainer);
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1CouplingMatrix


template <class LocalOperatorImp>
class Codim1BoundaryMatrix
{
  static_assert(std::is_base_of<LocalOperator::Codim1BoundaryInterface<typename LocalOperatorImp::Traits>,
                                LocalOperatorImp>::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim1BoundaryInterface!");

public:
  typedef LocalOperatorImp LocalOperatorType;

  Codim1BoundaryMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {
  }

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector<size_t> numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template <class T, class A, class IntersectionType, class M, class R>
  void assembleLocal(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace,
                     const IntersectionType& intersection, Dune::Stuff::LA::MatrixInterface<M>& systemMatrix,
                     std::vector<std::vector<Dune::DynamicMatrix<R>>>& tmpLocalMatricesContainer,
                     std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix
    Dune::DynamicMatrix<R>& localMatrix = tmpLocalMatricesContainer[0][0];
    localMatrix *= 0.0;
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity   = *entityPtr;
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpace.base_function_set(entity),
                         ansatzSpace.base_function_set(entity),
                         intersection,
                         localMatrix,
                         tmpOperatorMatrices);
    // write local matrices to global
    const size_t rows                       = testSpace.mapper().numDofs(entity);
    const size_t cols                       = ansatzSpace.mapper().numDofs(entity);
    Dune::DynamicVector<size_t>& globalRows = tmpIndicesContainer[0];
    Dune::DynamicVector<size_t>& globalCols = tmpIndicesContainer[1];
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


template <class LocalFunctionalImp>
class Codim1Vector
{
  static_assert(
      std::is_base_of<LocalFunctional::Codim1Interface<typename LocalFunctionalImp::Traits>, LocalFunctionalImp>::value,
      "LocalFunctionalImp has to be derived from LocalFunctional::Codim1Interface!");

public:
  typedef LocalFunctionalImp LocalFunctionalType;

  Codim1Vector(const LocalFunctionalType& fu)
    : localFunctional_(fu)
  {
  }

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector<size_t> numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired()};
  }

  template <class T, class IntersectionType, class V, class R>
  void assembleLocal(const SpaceInterface<T>& testSpace, const IntersectionType& intersection,
                     Dune::Stuff::LA::VectorInterface<V>& systemVector,
                     std::vector<std::vector<Dune::DynamicVector<R>>>& tmpLocalVectorsContainer,
                     Dune::DynamicVector<size_t>& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalVectorsContainer.size() >= 2);
    assert(tmpLocalVectorsContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorsContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    Dune::DynamicVector<R>& localVector = tmpLocalVectorsContainer[0][0];
    localVector *= 0.0;
    auto& tmpFunctionalVectors = tmpLocalVectorsContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity   = *entityPtr;
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


template <class GridViewImp, class LocalOperatorType, class TestFunctionType, class AnsatzFunctionType, class FieldType>
class Codim1BoundaryOperatorAccumulateFunctor : public Stuff::Grid::Functor::Codim1<GridViewImp>
{
  static_assert(std::is_base_of<LocalOperator::Codim1BoundaryInterface<typename LocalOperatorType::Traits>,
                                LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperator::Codim1BoundaryInterface!");
  static_assert(Stuff::is_localizable_function<TestFunctionType>::value,
                "TestFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<AnsatzFunctionType>::value,
                "AnsatzFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");

  typedef Stuff::Grid::Functor::Codim1<GridViewImp> BaseType;
  typedef DSC::TmpMatricesStorage<FieldType> TmpMatricesProviderType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  Codim1BoundaryOperatorAccumulateFunctor(const GridViewType& grd_vw, const LocalOperatorType& local_op,
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
    tmp_storage_ = std::unique_ptr<TmpMatricesProviderType>(
        new TmpMatricesProviderType({1, local_operator_.numTmpObjectsRequired()}, 1, 1));
  }

  virtual ~Codim1BoundaryOperatorAccumulateFunctor() = default;

  FieldType compute_locally(const IntersectionType& intersection, const EntityType& inside_entity,
                            const EntityType& /*outside_entity*/)
  {
    assert(tmp_storage_->matrices().size() >= 2);
    assert(tmp_storage_->matrices()[0].size() >= 1);
    auto& local_operator_result = tmp_storage_->matrices()[0][0];
    auto& tmp_matrices          = tmp_storage_->matrices()[1];
    // get the local functions
    const auto local_test_function   = test_function_.local_function(inside_entity);
    const auto local_ansatz_function = ansatz_function_.local_function(inside_entity);
    // apply the local operator
    this->local_operator_.apply(
        *local_test_function, *local_ansatz_function, intersection, local_operator_result, tmp_matrices);
    assert(local_operator_result.rows() >= 1);
    assert(local_operator_result.cols() >= 1);
    return local_operator_result[0][0];
  } // ... compute_locally(...)


  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity) DS_OVERRIDE
  {
    *result_ += compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() DS_OVERRIDE
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_        = true;
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
  DS::PerThreadValue<FieldType> result_;
  std::unique_ptr<TmpMatricesProviderType> tmp_storage_;
  bool finalized_;
  FieldType finalized_result_;
}; // class Codim1BoundaryOperatorAccumulateFunctor


} // namespace LocalAssembler
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
