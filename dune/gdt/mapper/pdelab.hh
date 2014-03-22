// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_PDELAB_HH
#define DUNE_GDT_MAPPER_PDELAB_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_PDELAB
# include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#endif

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class PdelabSpaceImp, int p = 1, int r = 1, int rR = 1 >
class PdelabPkQk
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "Not (yet) implemented for these dimensions!");
};


template< class PdelabSpaceImp, int p = 1, int r = 1, int rR = 1 >
class PdelabPkQkTraits
{
public:
  typedef PdelabPkQk< PdelabSpaceImp > derived_type;
  typedef PdelabSpaceImp BackendType;
};


template< class PdelabSpaceImp >
class PdelabPkQk< PdelabSpaceImp, 1, 1, 1 >
  : public MapperInterface< PdelabPkQkTraits< PdelabSpaceImp, 1, 1, 1 > >
{
  typedef MapperInterface< PdelabPkQkTraits< PdelabSpaceImp, 1, 1, 1 > > InterfaceType;
public:
  typedef PdelabPkQkTraits< PdelabSpaceImp, 1, 1, 1 > Traits;
  typedef typename Traits::BackendType                BackendType;
private:
  typedef PDELab::LocalFunctionSpace< BackendType, PDELab::TrialSpaceTag > PdeLabLFSType;

public:
  PdelabPkQk(const BackendType& pdelab_space)
    : backend_(pdelab_space)
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
  } // ... numDofs(...)

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
}; // class PdelabPkQk


#else // HAVE_DUNE_PDELAB


template< class PdelabSpaceImp, int p = 1, int r = 1, int rR = 1 >
class PdelabPkQk
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_PDELAB_HH
