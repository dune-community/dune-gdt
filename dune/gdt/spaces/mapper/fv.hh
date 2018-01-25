// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_MAPPER_FV_HH
#define DUNE_GDT_SPACES_MAPPER_FV_HH

#include <type_traits>

#include <dune/common/dynvector.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/xt/common/unused.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, size_t rangeDim = 1, size_t rangeDimCols = 1>
class FvMapper
{
  static_assert(AlwaysFalse<GridLayerImp>::value, "Not available for these dimensions!");
};

template <class GridLayerImp, size_t rangeDim = 1, size_t rangeDimCols = 1>
class FvProductMapper
{
  static_assert(AlwaysFalse<GridLayerImp>::value, "Not available for these dimensions!");
};


namespace internal {


template <class GridLayerImp, size_t rangeDim, size_t rangeDimCols>
class FvMapperTraits
{
  static_assert(rangeDim >= 1, "Really?");
  static_assert(rangeDimCols >= 1, "Really?");

  template <int dim_>
  struct GeometryTypeLayout
  {
    bool contains(const GeometryType& gt) const
    {
      return gt.dim() == dim_;
    }
  };

public:
  typedef GridLayerImp GridLayerType;
  typedef FvMapper<GridLayerType, rangeDim, rangeDimCols> derived_type;
  // just using the index set of the grid layer fails for mixed geometry types
  typedef MultipleCodimMultipleGeomTypeMapper<GridLayerImp, GeometryTypeLayout> BackendType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
};

template <class GridLayerImp, size_t rangeDim, size_t rangeDimCols>
class FvProductMapperTraits : public internal::FvMapperTraits<GridLayerImp, rangeDim, rangeDimCols>
{
public:
  typedef FvProductMapper<GridLayerImp, rangeDim, rangeDimCols> derived_type;
  static const size_t dimRange = rangeDim;
};


} // namespace internal


template <class GridLayerImp, size_t rangeDim>
class FvMapper<GridLayerImp, rangeDim, 1> : public MapperInterface<internal::FvMapperTraits<GridLayerImp, rangeDim, 1>>
{
  typedef MapperInterface<internal::FvMapperTraits<GridLayerImp, rangeDim, 1>> InterfaceType;
  static const size_t dimRange = rangeDim;

public:
  typedef internal::FvMapperTraits<GridLayerImp, rangeDim, 1> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  FvMapper(const GridLayerType& grd_layr)
    : mapper_(new BackendType(grd_layr))
  {
  }

  const BackendType& backend() const
  {
    return *mapper_;
  }

  size_t size() const
  {
    return dimRange * mapper_.size();
  }

  template <int cd, class GridImp, template <int, int, class> class EntityImp>
  size_t numDofs(const Entity<cd, EntityType::dimension, GridImp, EntityImp>& entity) const
  {
    return dimRange;
  }

  size_t maxNumDofs() const
  {
    return dimRange;
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    if (ret.size() < dimRange)
      ret.resize(dimRange);
    const size_t base = dimRange * mapper_->subIndex(entity, 0, 0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = base + ii;
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < dimRange);
    return (dimRange * mapper_->subIndex(entity, 0, 0)) + localIndex;
  }

private:
  const std::shared_ptr<BackendType> mapper_;
}; // class FvMapper< ..., rangeDim, 1 >


template <class GridLayerImp>
class FvMapper<GridLayerImp, 1, 1> : public MapperInterface<internal::FvMapperTraits<GridLayerImp, 1, 1>>
{
  typedef MapperInterface<internal::FvMapperTraits<GridLayerImp, 1, 1>> InterfaceType;

public:
  typedef internal::FvMapperTraits<GridLayerImp, 1, 1> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  FvMapper(const GridLayerType& grd_layr)
    : mapper_(new BackendType(grd_layr))
  {
  }

  const BackendType& backend() const
  {
    return *mapper_;
  }

  size_t size() const
  {
    return mapper_->size();
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

  size_t mapToGlobal(const EntityType& entity, const size_t& DXTC_DEBUG_ONLY(localIndex)) const
  {
    assert(localIndex == 0);
    return mapper_->subIndex(entity, 0, 0);
  }

private:
  const std::shared_ptr<BackendType> mapper_;
}; // class FvMapper< ..., 1, 1 >


template <class GridLayerImp, size_t rangeDim>
class FvProductMapper<GridLayerImp, rangeDim, 1>
    : public ProductMapperInterface<internal::FvProductMapperTraits<GridLayerImp, rangeDim, 1>>
{
  typedef ProductMapperInterface<internal::FvProductMapperTraits<GridLayerImp, rangeDim, 1>> BaseType;
  typedef FvMapper<GridLayerImp, rangeDim, 1> FvMapperMapperType;

public:
  typedef internal::FvProductMapperTraits<GridLayerImp, rangeDim, 1> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  static const size_t dimRange = Traits::dimRange;
  using typename BaseType::EntityType;
  using typename BaseType::BackendType;

  FvProductMapper(const GridLayerType& grd_layr)
    : fv_mapper_(grd_layr)
  {
  }

  // These methods are required by the ProductMapperInterface
  size_t numDofs(const size_t /*factor_index*/, const EntityType& /*entity*/) const
  {
    return 1;
  }

  size_t maxNumDofs(const size_t /*factor_index*/) const
  {
    return 1;
  }

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    if (ret.size() != 1)
      ret.resize(1);
    ret[0] = dimRange * (backend().index(entity)) + factor_index;
  }

  size_t mapToGlobal(const size_t factor_index,
                     const EntityType& entity,
                     const size_t& DXTC_DEBUG_ONLY(local_index_in_factor)) const
  {
    assert(local_index_in_factor == 0);
    assert(factor_index < numDofs(entity));
    return dimRange * (backend().index(entity)) + factor_index;
  }

  size_t mapToLocal(const size_t factor_index,
                    const EntityType& DXTC_DEBUG_ONLY(entity),
                    const size_t& DXTC_DEBUG_ONLY(local_index_in_factor)) const
  {
    assert(local_index_in_factor == 0);
    assert(factor_index < numDofs(entity));
    return factor_index;
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

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    return fv_mapper_.globalIndices(entity, ret);
  } // ... globalIndices(...)

  using BaseType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    return fv_mapper_.mapToGlobal(entity, localIndex);
  }

private:
  const FvMapperMapperType fv_mapper_;
}; // class FvProductMapper< ..., rangeDim, 1 >


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_FV_HH
