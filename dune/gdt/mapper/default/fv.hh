// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_DEFAULT_FV_HH
#define DUNE_GDT_MAPPER_DEFAULT_FV_HH

#include <type_traits>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/stuff/common/debug.hh>

#include "../../mapper/interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {


// forward, to be used in the traits
template< class GridViewImp, int rangeDim = 1, int rangeDimCols = 1 >
class FiniteVolume;


template< class GridViewImp, int rangeDim, int rangeDimCols >
class FiniteVolumeTraits
{
  static_assert(rangeDim >= 1, "Really?");
  static_assert(rangeDimCols >= 1, "Really?");
public:
  typedef GridViewImp GridViewType;
  typedef FiniteVolume< GridViewType, rangeDim, rangeDimCols> derived_type;
  typedef typename GridViewImp::IndexSet BackendType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
};


template< class GridViewImp >
class FiniteVolume< GridViewImp, 1, 1 >
  : public MapperInterface< FiniteVolumeTraits< GridViewImp, 1, 1 > >
{
  typedef MapperInterface< FiniteVolumeTraits< GridViewImp, 1, 1 > > InterfaceType;
public:
  typedef FiniteVolumeTraits< GridViewImp, 1, 1 > Traits;
  typedef typename Traits::GridViewType           GridViewType;
  typedef typename Traits::BackendType            BackendType;
  typedef typename Traits::EntityType             EntityType;

  FiniteVolume(const GridViewType& grid_view)
    : backend_(grid_view.indexSet())
  {}

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

  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    if (ret.size() < 1)
      ret.resize(1);
    ret[0] = mapToGlobal(entity, 0);
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& UNUSED_UNLESS_DEBUG(localIndex)) const
  {
    assert(localIndex == 0);
    return backend_.index(entity);
  } // ... mapToGlobal(...)

private:
  const BackendType& backend_;
}; // class FiniteVolume< ..., 1, 1 >


template< class GridViewImp, int rangeDim >
class FiniteVolume< GridViewImp, rangeDim, 1 >
  : public MapperInterface< FiniteVolumeTraits< GridViewImp, rangeDim, 1 > >
{
  typedef MapperInterface< FiniteVolumeTraits< GridViewImp, rangeDim, 1 > > InterfaceType;
  static const unsigned int dimRange = rangeDim;
public:
  typedef FiniteVolumeTraits< GridViewImp, rangeDim, 1 >  Traits;
  typedef typename Traits::GridViewType                   GridViewType;
  typedef typename Traits::BackendType                    BackendType;
  typedef typename Traits::EntityType                     EntityType;

  FiniteVolume(const GridViewType& grid_view)
    : backend_(grid_view.indexSet())
  {}

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return dimRange * backend_.size(0);
  }

  size_t numDofs(const EntityType& /*entity*/) const
  {
    return dimRange;
  }

  size_t maxNumDofs() const
  {
    return dimRange;
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    if (ret.size() < dimRange)
      ret.resize(dimRange);
    const size_t base = dimRange * backend_.index(entity);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = base + ii;
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < dimRange);
    return (dimRange * backend_.index(entity)) + localIndex;
  } // ... mapToGlobal(...)

private:
  const BackendType& backend_;
}; // class FiniteVolume< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_DEFAULT_FV_HH
