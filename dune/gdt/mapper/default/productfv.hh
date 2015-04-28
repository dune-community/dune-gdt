// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_DEFAULT_PRODUCTFV_HH
#define DUNE_GDT_MAPPER_DEFAULT_PRODUCTFV_HH

#include <dune/common/dynvector.hh>

#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/type_utils.hh>

#include "../../mapper/interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {


// forward
template< class GridViewImp, size_t rangeDim = 1, size_t rangeDimCols = 1 >
class ProductFiniteVolume
{
  static_assert(AlwaysFalse< GridViewImp >::value, "Not available for these dimensions!");
};


namespace internal {


template< class GridViewImp, size_t rangeDim, size_t rangeDimCols >
class ProductFiniteVolumeTraits
{
  static_assert(rangeDim >= 1, "Really?");
  static_assert(rangeDimCols >= 1, "Really?");
public:
  typedef GridViewImp GridViewType;
  typedef ProductFiniteVolume< GridViewType, rangeDim, rangeDimCols> derived_type;
  typedef typename GridViewImp::IndexSet BackendType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
};


} // namespace internal

template< class GridViewImp, size_t rangeDim >
class ProductFiniteVolume< GridViewImp, rangeDim, 1 >
  : public MapperInterface< internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > >
{
  typedef MapperInterface< internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > > InterfaceType;
  static const size_t dimRange = rangeDim;
public:
  typedef internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > Traits;
  typedef typename Traits::GridViewType                            GridViewType;
  typedef typename Traits::BackendType                             BackendType;
  typedef typename Traits::EntityType                              EntityType;

  ProductFiniteVolume(const GridViewType& grid_view)
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

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    assert(factor_index < dimRange);
    if (ret.size() < 1)
      ret.resize(1);
    const size_t base = dimRange * backend_.index(entity);
      ret[0] = base + factor_index;
  } // ... globalIndices(...)

  Dune::DynamicVector< size_t > globalIndices(const size_t factor_index, const EntityType& entity) const
  {
    Dune::DynamicVector< size_t > ret(numDofs(entity), 0);
    globalIndices(factor_index, entity, ret);
    return ret;
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < numDofs(entity));
    return (dimRange * backend_.index(entity)) + localIndex;
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& UNUSED_UNLESS_DEBUG(localIndex)) const
  {
    assert(localIndex == 0);
    assert(factor_index < numDofs(entity));
    return (dimRange * backend_.index(entity)) + factor_index;
  }

private:
  const BackendType& backend_;
}; // class ProductFiniteVolume< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_DEFAULT_PRODUCTFV_HH
