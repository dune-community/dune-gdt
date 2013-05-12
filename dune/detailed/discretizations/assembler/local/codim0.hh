#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/la/container/interface.hh>

#include <dune/detailed/discretizations/localoperator/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim0.hh>
#include <dune/detailed/discretizations/space/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits
template< class LocalOperatorImp >
class LocalAssemblerCodim0Matrix;


template< class LocalOperatorImp >
class LocalAssemblerCodim0MatrixTraits
{
public:
  typedef LocalAssemblerCodim0Matrix< LocalOperatorImp > derived_type;
  typedef LocalOperatorCodim0Interface< typename LocalOperatorImp::Traits > LocalOperatorType;
}; // class LocalAssemblerCodim0MatrixTraits


template< class LocalOperatorImp >
class LocalAssemblerCodim0Matrix
{
public:
  typedef LocalAssemblerCodim0MatrixTraits< LocalOperatorImp > Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  LocalAssemblerCodim0Matrix(const LocalOperatorType& op)
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
    return { numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired() };
  }

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of testSpace
   *  \tparam A           Traits of the SpaceInterface implementation, representing the type of ansatzSpace
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam M           Traits of the Dune::Stuff::LA::Container::MatrixInterface implementation, representing the type of systemMatrix
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class T, class A, class EntityType, class M, class R >
  void assembleLocal(const SpaceInterface< T >& testSpace,
                     const SpaceInterface< A >& ansatzSpace,
                     const EntityType& entity,
                     Dune::Stuff::LA::Container::MatrixInterface< M >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
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
    localOperator_.apply(testSpace.baseFunctionSet(entity),
                         ansatzSpace.baseFunctionSet(entity),
                         localMatrix,
                         tmpOperatorMatrices);
    // write local matrix to global
    auto& globalRows = tmpIndicesContainer[0];
    auto& globalCols = tmpIndicesContainer[1];
    const size_t rows = testSpace.mapper().numDofs(entity);
    const size_t cols = ansatzSpace.mapper().numDofs(entity);
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    testSpace.mapper().globalIndices(entity, globalRows);
    ansatzSpace.mapper().globalIndices(entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localRow = localMatrix[ii];
      for (size_t jj = 0; jj < cols; ++jj)
        systemMatrix.add(globalRows[ii], globalCols[jj], localRow[jj]);
    } // write local matrix to global
  } // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
}; // class LocalAssemblerCodim0Matrix


// forward, to be used in the traits
template< class LocalOperatorImp >
class LocalAssemblerCodim0Vector;


template< class LocalFunctionalImp >
class LocalAssemblerCodim0VectorTraits
{
public:
  typedef LocalAssemblerCodim0Vector< LocalFunctionalImp > derived_type;
  typedef LocalFunctionalCodim0Interface< typename LocalFunctionalImp::Traits > LocalFunctionalType;
}; // class LocalAssemblerCodim0MatrixTraits


template< class LocalFunctionalImp >
class LocalAssemblerCodim0Vector
{
public:
  typedef LocalAssemblerCodim0VectorTraits< LocalFunctionalImp > Traits;
  typedef typename Traits::LocalFunctionalType LocalFunctionalType;

  LocalAssemblerCodim0Vector(const LocalFunctionalType& func)
    : localFunctional_(func)
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
    return { numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired() };
  }

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of testSpace
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam V           Traits of the Dune::Stuff::LA::Container::VectorInterface implementation, representing the type of systemVector
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class T, class EntityType, class V, class R >
  void assembleLocal(const SpaceInterface< T >& testSpace,
                     const EntityType& entity,
                     Dune::Stuff::LA::Container::VectorInterface< V >& systemVector,
                     std::vector< std::vector< Dune::DynamicVector< R > > >& tmpLocalVectorContainer,
                     Dune::DynamicVector< size_t >& tmpIndices) const
  {
    // check
    assert(tmpLocalVectorContainer.size() >= 1);
    assert(tmpLocalVectorContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    auto& localVector = tmpLocalVectorContainer[0][0];
    Dune::Stuff::Common::clear(localVector);
    auto& tmpFunctionalVectors = tmpLocalVectorContainer[1];
    // apply local functional (result is in localVector)
    localFunctional_.apply(testSpace.baseFunctionSet(entity),
                           localVector,
                           tmpFunctionalVectors);
    // write local vector to global
    const size_t size = testSpace.mapper().numDofs(entity);
    assert(tmpIndices.size() >= size);
    testSpace.mapper().globalIndices(entity, tmpIndices);
    for (size_t ii = 0; ii < size; ++ii) {
      systemVector.add(tmpIndices[ii], localVector[ii]);
    } // write local matrix to global
  } // ... assembleLocal(...)

private:
  const LocalFunctionalType& localFunctional_;
}; // class LocalAssemblerCodim0Vector


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_HH
