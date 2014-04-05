// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_FINITEVOLUME_HH
#define DUNE_GDT_MAPPER_FINITEVOLUME_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include "../../mapper/interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {


// forward, to be used in the traits
template <class GridViewImp>
class FiniteVolume;


template <class GridViewImp>
class FiniteVolumeTraits
{
public:
  typedef GridViewImp GridViewType;
  typedef FiniteVolume<GridViewType> derived_type;
  typedef typename GridViewImp::IndexSet BackendType;
};


template <class GridViewImp>
class FiniteVolume : public MapperInterface<FiniteVolumeTraits<GridViewImp>>
{
  typedef MapperInterface<FiniteVolumeTraits<GridViewImp>> InterfaceType;

public:
  typedef FiniteVolumeTraits<GridViewImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::BackendType BackendType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  FiniteVolume(const GridViewType& grid_view)
    : backend_(grid_view.indexSet())
  {
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size(0);
  }

  size_t numDofs(const EntityType& /*entity*/) const
  {
    return 1;
  }

  size_t maxNumDofs() const
  {
    return 1;
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    if (ret.size() < 1)
      ret.resize(1);
    ret[0] = mapToGlobal(entity, 0);
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex == 0);
    return backend_.index(entity);
  } // ... mapToGlobal(...)

private:
  const BackendType& backend_;
}; // class FiniteVolume


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_FINITEVOLUME_HH
