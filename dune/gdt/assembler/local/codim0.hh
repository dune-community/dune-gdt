// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_CODIM0_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_CODIM0_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/localfunctional/interface.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalAssembler {


// forward, to be used in the traits
template< class LocalOperatorImp >
class Codim0Matrix;


template< class LocalOperatorImp >
class Codim0MatrixTraits
{
public:
  typedef Codim0Matrix< LocalOperatorImp > derived_type;
  typedef LocalOperator::Codim0Interface< typename LocalOperatorImp::Traits > LocalOperatorType;
}; // class LocalAssemblerCodim0MatrixTraits


template< class LocalOperatorImp >
class Codim0Matrix
{
public:
  typedef Codim0MatrixTraits< LocalOperatorImp > Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  Codim0Matrix(const LocalOperatorType& op)
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
                     Dune::Stuff::LA::MatrixInterface< M >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    // check
    assert(tmpLocalMatricesContainer.size() >= 1);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix
    Dune::DynamicMatrix< R >& localMatrix = tmpLocalMatricesContainer[0][0];
    Dune::Stuff::Common::clear(localMatrix);
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // apply local operator (result is in localMatrix)
    localOperator_.apply(testSpace.base_function_set(entity),
                         ansatzSpace.base_function_set(entity),
                         localMatrix,
                         tmpOperatorMatrices);
    // write local matrix to global
    Dune::DynamicVector< size_t >& globalRows = tmpIndicesContainer[0];
    Dune::DynamicVector< size_t >& globalCols = tmpIndicesContainer[1];
    const size_t rows = testSpace.mapper().numDofs(entity);
    const size_t cols = ansatzSpace.mapper().numDofs(entity);
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    testSpace.mapper().globalIndices(entity, globalRows);
    ansatzSpace.mapper().globalIndices(entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localRow = localMatrix[ii];
      const size_t globalII = globalRows[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t globalJJ = globalCols[jj];
        systemMatrix.add_to_entry(globalII, globalJJ, localRow[jj]);
      }
    } // write local matrix to global
  } // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
}; // class LocalAssemblerCodim0Matrix


// forward, to be used in the traits
template< class LocalOperatorImp >
class Codim0Vector;


template< class LocalFunctionalImp >
class Codim0VectorTraits
{
public:
  typedef Codim0Vector< LocalFunctionalImp > derived_type;
  typedef LocalFunctional::Codim0Interface< typename LocalFunctionalImp::Traits > LocalFunctionalType;
}; // class LocalAssemblerCodim0MatrixTraits


template< class LocalFunctionalImp >
class Codim0Vector
{
public:
  typedef Codim0VectorTraits< LocalFunctionalImp > Traits;
  typedef typename Traits::LocalFunctionalType LocalFunctionalType;

  Codim0Vector(const LocalFunctionalType& func)
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
                     Dune::Stuff::LA::VectorInterface< V >& systemVector,
                     std::vector< std::vector< Dune::DynamicVector< R > > >& tmpLocalVectorContainer,
                     Dune::DynamicVector< size_t >& tmpIndices) const
  {
    // check
    assert(tmpLocalVectorContainer.size() >= 2);
    assert(tmpLocalVectorContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    auto& localVector = tmpLocalVectorContainer[0][0];
    Dune::Stuff::Common::clear(localVector);
    auto& tmpFunctionalVectors = tmpLocalVectorContainer[1];
    // apply local functional (result is in localVector)
    localFunctional_.apply(testSpace.base_function_set(entity),
                           localVector,
                           tmpFunctionalVectors);
    // write local vector to global
    const size_t size = testSpace.mapper().numDofs(entity);
    assert(tmpIndices.size() >= size);
    testSpace.mapper().globalIndices(entity, tmpIndices);
    for (size_t ii = 0; ii < size; ++ii) {
      systemVector.add_to_entry(tmpIndices[ii], localVector[ii]);
    } // write local matrix to global
  } // ... assembleLocal(...)

private:
  const LocalFunctionalType& localFunctional_;
}; // class Codim0Vector


} // namespace LocalAssembler
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_CODIM0_HH
