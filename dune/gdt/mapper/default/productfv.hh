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
#include "fv.hh"

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
    : public internal::FiniteVolumeTraits< GridViewImp, rangeDim, rangeDimCols >
{
public:
  typedef ProductFiniteVolume< GridViewImp, rangeDim, rangeDimCols> derived_type;
  static const size_t dimRange = rangeDim;
};


} // namespace internal


template< class GridViewImp, size_t rangeDim >
class ProductFiniteVolume< GridViewImp, rangeDim, 1 >
  : public ProductMapperInterface< internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > >
{
  typedef ProductMapperInterface< internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > > BaseType;
  typedef FiniteVolume< GridViewImp, rangeDim, 1 > FiniteVolumeMapperType;
public:
  typedef internal::ProductFiniteVolumeTraits< GridViewImp, rangeDim, 1 > Traits;
  typedef typename Traits::GridViewType GridViewType;
  static const size_t dimRange = Traits::dimRange;
  using typename BaseType::EntityType;
  using typename BaseType::BackendType;

  ProductFiniteVolume(const GridViewType& grid_view)
    : fv_mapper_(grid_view)
  {}

  // These methods are required by the ProductMapperInterface
  Dune::DynamicVector< size_t > globalIndices(const size_t factor_index, const EntityType& entity) const
  {
    Dune::DynamicVector< size_t > ret(numDofs(entity), 0);
    globalIndices(factor_index, entity, ret);
    return ret;
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& UNUSED_UNLESS_DEBUG(localIndex)) const
  {
    assert(localIndex == 0);
    assert(factor_index < numDofs(entity));
    return (dimRange*(backend().index(entity))) + factor_index;
  }

  // The remaining methods are just redirected to the usual finite volume mapper
  const BackendType& backend() const
  {
    return fv_mapper_.backend();
  }

  size_t size() const
  {
    return fv_mapper_.size();
  }

  size_t numDofs(const EntityType& entity) const
  {
    return fv_mapper_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return fv_mapper_.maxNumDofs();
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    return fv_mapper_.globalIndices(entity, ret);
  } // ... globalIndices(...)

  using BaseType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    return fv_mapper_.mapToGlobal(entity, localIndex);
  }

private:
  const FiniteVolumeMapperType fv_mapper_;
}; // class ProductFiniteVolume< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_DEFAULT_PRODUCTFV_HH
