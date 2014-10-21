// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_PDELAB_HH
#define DUNE_GDT_MAPPER_PDELAB_HH

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_PDELAB
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/common/parallel/threadmanager.hh>
#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits
template< class PdelabSpaceImp >
class SimplePdelabWrapper;


template< class PdelabSpaceImp >
class SimplePdelabWrapperTraits
{
public:
  typedef SimplePdelabWrapper< PdelabSpaceImp > derived_type;
  typedef PdelabSpaceImp BackendType;
  typedef typename BackendType::Element EntityType;
};


template< class PdelabSpaceImp >
class SimplePdelabWrapper
  : public MapperInterface< SimplePdelabWrapperTraits< PdelabSpaceImp > >
{
  typedef MapperInterface< SimplePdelabWrapperTraits< PdelabSpaceImp > > InterfaceType;
public:
  typedef SimplePdelabWrapperTraits< PdelabSpaceImp >  Traits;
  typedef typename Traits::BackendType                BackendType;
private:
  typedef PDELab::LocalFunctionSpace< BackendType, PDELab::TrialSpaceTag > PdeLabLFSType;

public:
  explicit SimplePdelabWrapper(const BackendType& pdelab_space)
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

  void globalIndices(const typename Traits::EntityType& entity, Dune::DynamicVector< size_t >& ret) const
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
    assert(localIndex < lfs_.size());
    return lfs_.dofIndex(localIndex).entityIndex()[1];
  } // ... mapToGlobal(...)

private:
  const BackendType& backend_;
  mutable PdeLabLFSType lfs_;
}; // class SimplePdelabWrapper


#else // HAVE_DUNE_PDELAB


template< class PdelabSpaceImp >
class SimplePdelabWrapper
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_PDELAB_HH
