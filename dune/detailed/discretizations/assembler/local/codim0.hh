#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/la/container/interface.hh>

#include <dune/detailed/discretizations/localoperator/codim0.hh>
#include <dune/detailed/discretizations/space/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template <class LocalOperatorImp>
class LocalAssemblerCodim0Matrix;


template <class LocalOperatorImp>
class LocalAssemblerCodim0MatrixTraits
{
public:
  typedef LocalAssemblerCodim0Matrix<LocalOperatorImp> derived_type;
  typedef LocalOperatorCodim0Interface<typename LocalOperatorImp::Traits> LocalOperatorType;
}; // class LocalAssemblerCodim0MatrixTraits


template <class LocalOperatorImp>
class LocalAssemblerCodim0Matrix
{
public:
  typedef LocalAssemblerCodim0MatrixTraits<LocalOperatorImp> Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  LocalAssemblerCodim0Matrix(const LocalOperatorType& op)
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

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of testSpace
   *  \tparam A
   *  \tparam EntityType
   *  \tparam M
   *  \tparam L
   */
  template <class T, class A, class EntityType, class M, class R>
  void assembleLocal(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace, const EntityType& entity,
                     Dune::Stuff::LA::Container::MatrixInterface<M>& systemMatrix,
                     std::vector<std::vector<Dune::DynamicMatrix<R>>>& tmpLocalMatricesContainer,
                     std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 1);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix
    auto& localMatrix = tmpLocalMatricesContainer[0][0];
    Dune::Stuff::Common::clear(localMatrix);
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // apply local operator (result is in localMatrix)
    localOperator_.apply(
        testSpace.baseFunctionSet(entity), ansatzSpace.baseFunctionSet(entity), localMatrix, tmpOperatorMatrices);
    // write local matrix to global
    auto& globalRows  = tmpIndicesContainer[0];
    auto& globalCols  = tmpIndicesContainer[1];
    const size_t rows = testSpace.mapper().numDofs(entity);
    const size_t cols = ansatzSpace.mapper().numDofs(entity);
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    testSpace.mapper().globalIndices(entity, globalRows);
    ansatzSpace.mapper().globalIndices(entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localRow = localMatrix[ii];
      for (size_t jj = 0; jj < cols; ++jj)
        systemMatrix.set(globalRows[ii], globalCols[jj], localRow[jj]);
    } // write local matrix to global
  } // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
}; // class LocalAssemblerCodim0Matrix

} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH
