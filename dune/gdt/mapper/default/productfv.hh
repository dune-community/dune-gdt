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
template <class GridViewImp, size_t rangeDim = 1, size_t rangeDimCols = 1>
class ProductFiniteVolume
{
  static_assert(AlwaysFalse<GridViewImp>::value, "Not available for these dimensions!");
};


namespace internal {


template <class GridViewImp, size_t rangeDim, size_t rangeDimCols>
class ProductFiniteVolumeTraits : internal::FiniteVolumeTraits<GridViewImp, rangeDim, rangeDimCols>
{
public:
  typedef ProductFiniteVolume<GridViewType, rangeDim, rangeDimCols> derived_type;
};


} // namespace internal


template <class GridViewImp, size_t rangeDim>
class ProductFiniteVolume<GridViewImp, rangeDim, 1>
    : public ProductMapperInterface<internal::ProductFiniteVolumeTraits<GridViewImp, rangeDim, 1>>,
      public FiniteVolume<GridViewImp, rangeDim, 1>
{
  typedef ProductMapperInterface<internal::ProductFiniteVolumeTraits<GridViewImp, rangeDim, 1>> InterfaceType;
  typedef FiniteVolume<GridViewImp, rangeDim, 1> BaseType;
  static const size_t dimRange = rangeDim;

public:
  typedef internal::ProductFiniteVolumeTraits<GridViewImp, rangeDim, 1> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  ProductFiniteVolume(const GridViewType& grid_view)
    : BaseType(grid_view)
  {
  }

  using BaseType::backend;
  using BaseType::size;
  using BaseType::numDofs;
  using BaseType::maxNumDofs;
  using BaseType::globalIndices;

  Dune::DynamicVector<size_t> globalIndices(const size_t factor_index, const EntityType& entity) const
  {
    Dune::DynamicVector<size_t> ret(numDofs(entity), 0);
    globalIndices(factor_index, entity, ret);
    return ret;
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity,
                     const size_t& UNUSED_UNLESS_DEBUG(localIndex)) const
  {
    assert(localIndex == 0);
    assert(factor_index < numDofs(entity));
    return (dimRange * (backend_.index(entity))) + factor_index;
  }
}; // class ProductFiniteVolume< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_DEFAULT_PRODUCTFV_HH
