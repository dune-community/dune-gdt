// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_INTERFACE_HH
#define DUNE_GDT_MAPPER_INTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynvector.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class MapperInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BackendType  BackendType;

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  size_t size() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  size_t maxNumDofs() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().maxNumDofs());
    return asImp().maxNumDofs();
  }

  template< class EntityType >
  size_t numDofs(const EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs(entity));
    return asImp().numDofs(entity);
  }

  template< class EntityType >
  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().globalIndices(entity, ret));
  }

  template< class EntityType >
  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity, localIndex));
    return asImp().mapToGlobal(entity, localIndex);
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class MapperInterface


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_MAPPER_INTERFACE_HH
