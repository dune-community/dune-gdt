// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_VIRTUAL_CODIM0_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_VIRTUAL_CODIM0_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>
#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/walker/functors.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/localfunctional/interface.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalAssembler {


template< class LocalOperatorImp >
class VirtualCodim0Matrix
{
  static_assert(std::is_base_of< LocalOperator::Codim0Interface< typename LocalOperatorImp::Traits >,
                                 LocalOperatorImp >::value,
                "LocalOperatorImp has to be derived from LocalOperator::Codim0Interface!");
public:
  typedef LocalOperatorImp LocalOperatorType;

  explicit VirtualCodim0Matrix(const LocalOperatorType& op)
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
   *  \tparam *d          dimDomain of testSpace (* == T) or ansatzSpace (* == A)
   *  \tparam *r          dimRange of testSpace (* == T) or ansatzSpace (* == A)
   *  \tparam *rC         dimRangeCols of testSpace (* == T) or ansatzSpace (* == A)
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam M           Traits of the Dune::Stuff::LA::Container::MatrixInterface implementation, representing the type of systemMatrix
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class T, size_t Td, size_t Tr, size_t TrC, class A, size_t Ad, size_t Ar, size_t ArC, class EntityType, class M, class R >
  void assembleLocal(const SpaceInterface< T, Td, Tr, TrC >& testSpace,
                     const SpaceInterface< A, Ad, Ar, ArC >& ansatzSpace,
                     const EntityType& entity,
                     Dune::Stuff::LA::MatrixInterface< M, R >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {

  } // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
}; // class VirtualCodim0Matrix


template< class LocalFunctionalImp >
class VirtualCodim0Vector
{
  static_assert(std::is_base_of< LocalFunctional::Codim0Interface< typename LocalFunctionalImp::Traits >,
                                 LocalFunctionalImp >::value,
                "LocalFunctionalImp has to be derived from LocalFunctional::Codim0Interface!");
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  explicit VirtualCodim0Vector(const LocalFunctionalType& func)
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
   *  \tparam d           dimDomain of testSpace
   *  \tparam r           dimRange of testSpace
   *  \tparam rC          dimRangeCols of testSpace
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam V           Traits of the Dune::Stuff::LA::Container::VectorInterface implementation, representing the type of systemVector
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class T, size_t d, size_t r, size_t rC, class EntityType, class V, class R >
  void assembleLocal(const SpaceInterface< T, d, r, rC >& testSpace,
                     const EntityType& entity,
                     Dune::Stuff::LA::VectorInterface< V, R >& systemVector,
                     std::vector< std::vector< Dune::DynamicVector< R > > >& tmpLocalVectorContainer,
                     Dune::DynamicVector< size_t >& tmpIndices) const
  {
    // check
    assert(tmpLocalVectorContainer.size() >= 2);
    assert(tmpLocalVectorContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    auto& localVector = tmpLocalVectorContainer[0][0];
    localVector *= 0.0;
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
}; // class VirtualCodim0Vector


} // namespace LocalAssembler
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_VIRTUAL_CODIM0_HH
