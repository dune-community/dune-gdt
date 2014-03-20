// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_PDELAB_HH
#define DUNE_GDT_MAPPER_PDELAB_HH

#include <dune/common/dynvector.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {


//// forward, to be used in the traits and to allow for specialization
template< class PdelabSpaceImp >
class PdelabWrapper;


template< class PdelabSpaceImp >
class PdelabWrapperTraits
{
public:
  typedef PdelabWrapper< PdelabSpaceImp > derived_type;
  typedef PdelabSpaceImp BackendType;
};


template< class PdelabSpaceType >
class PdelabWrapper
  : public MapperInterface< PdelabWrapperTraits< PdelabSpaceType > >
{
  typedef MapperInterface< PdelabWrapperTraits< PdelabSpaceType > > InterfaceType;
public:
  typedef PdelabWrapperTraits< PdelabSpaceType > Traits;
  typedef typename Traits::BackendType           BackendType;
private:
  typedef PDELab::LocalFunctionSpace< BackendType, PDELab::TrialSpaceTag > PdeLabLFSType;

public:
  PdelabWrapper(const BackendType& space)
    : backend_(space)
    , lfs_(backend_)
  {}

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  template< class EntityType >
  size_t numDofs(const EntityType& entity) const
  {
    lfs_.bind(entity);
    return lfs_.size();
  }

  size_t maxNumDofs() const
  {
    return backend_.maxLocalSize();
  }

  template< class EntityType >
  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    lfs_.bind(entity);
    // some checks
    const size_t numLocalDofs = numDofs(entity);
    if (ret.size() < numLocalDofs)
      ret.resize(numLocalDofs);
    // compute
    for (size_t ii = 0; ii < numLocalDofs; ++ii)
      ret[ii] = mapToGlobal(entity, ii);
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  template< class EntityType >
  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    lfs_.bind(entity);
    const size_t numLocalDofs = numDofs(entity);
    assert(localIndex < numLocalDofs);
    return lfs_.dofIndex(localIndex).entityIndex()[1];
  } // ... mapToGlobal(...)

private:
  const BackendType& backend_;
  mutable PdeLabLFSType lfs_;
}; // class PdelabWrapper


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_PDELAB_HH
