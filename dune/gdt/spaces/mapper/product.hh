// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_MAPPER_PRODUCT_HH
#define DUNE_GDT_SPACES_MAPPER_PRODUCT_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/gdt/spaces/basefunctionset/product.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class... MapperTypes>
class DefaultProductMapper;


namespace internal {


template <size_t I = 0>
struct DynamicTupleGetter
{
  template <class EntityType, typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), size_t>::type
  numDofs(const std::tuple<TupleArgs...>& /*tuple*/, const size_t /*factor_index*/, const EntityType& /*entity*/)
  {
    return 0;
  }

  template <class EntityType, typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), size_t>::type
        numDofs(const std::tuple<TupleArgs...>& tuple, const size_t factor_index, const EntityType& entity)
  {
    if (factor_index == 0)
      return std::get<I>(tuple).numDofs(entity);
    return DynamicTupleGetter<I + 1>::numDofs(tuple, factor_index - 1, entity);
  }

  template <typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), size_t>::type
  maxNumDofs(const std::tuple<TupleArgs...>& /*tuple*/, const size_t /*factor_index*/)
  {
    return 0;
  }

  template <typename... TupleArgs>
      static typename std::enable_if < I<sizeof...(TupleArgs), size_t>::type
                                       maxNumDofs(const std::tuple<TupleArgs...>& tuple, const size_t factor_index)
  {
    if (factor_index == 0)
      return std::get<I>(tuple).maxNumDofs();
    return DynamicTupleGetter<I + 1>::maxNumDofs(tuple, factor_index - 1);
  }

  template <typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), size_t>::type
  size(const std::tuple<TupleArgs...>& /*tuple*/, const size_t /*factor_index*/)
  {
    return 0;
  }

  template <typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), size_t>::type size(const std::tuple<TupleArgs...>& tuple, const size_t factor_index)
  {
    if (factor_index == 0)
      return std::get<I>(tuple).size();
    return DynamicTupleGetter<I + 1>::size(tuple, factor_index - 1);
  }

  template <class EntityType, typename... TupleArgs>
  static typename std::enable_if<I == sizeof...(TupleArgs), size_t>::type
  mapToGlobal(const std::tuple<TupleArgs...>& /*tuple*/,
              const size_t /*factor_index*/,
              const EntityType& /*entity*/,
              const size_t& /*local_index_in_factor*/)
  {
    return 0;
  }

  template <class EntityType, typename... TupleArgs>
      static typename std::enable_if
      < I<sizeof...(TupleArgs), size_t>::type mapToGlobal(const std::tuple<TupleArgs...>& tuple,
                                                          const size_t factor_index,
                                                          const EntityType& entity,
                                                          const size_t& local_index_in_factor)
  {
    if (factor_index == 0)
      return std::get<I>(tuple).mapToGlobal(entity, local_index_in_factor);
    return DynamicTupleGetter<I + 1>::mapToGlobal(tuple, factor_index - 1, entity, local_index_in_factor);
  }
};

template <size_t I, class... SpaceTypes>
struct MapperTuplefromSpaceTupleCreator
{
  template <class... MapperTypes>
  static typename std::enable_if<sizeof...(MapperTypes) == sizeof...(SpaceTypes),
                                 typename std::tuple<MapperTypes...>>::type
  create(const std::tuple<SpaceTypes...>& /*spaces*/, const std::tuple<MapperTypes...>& mappers)
  {
    return mappers;
  }

  template <class... MapperTypes>
  static typename std::enable_if<sizeof...(MapperTypes) < sizeof...(SpaceTypes),
                                 typename std::tuple<typename SpaceTypes::MapperType...>>::type
  create(const std::tuple<SpaceTypes...>& spaces, const std::tuple<MapperTypes...>& mappers)
  {
    const auto new_mappers = std::tuple_cat(mappers, std::make_tuple(std::get<I>(spaces).mapper()));
    return MapperTuplefromSpaceTupleCreator<I + 1, SpaceTypes...>::create(spaces, new_mappers);
  }
};


template <class GridLayerImp, class... MapperTypes>
class DefaultProductMapperTraits
{
public:
  typedef GridLayerImp GridLayerType;
  typedef DefaultProductMapper<GridLayerType, MapperTypes...> derived_type;
  static const size_t dimRange = GDT::BaseFunctionSet::internal::SumDimRange<MapperTypes...>::dimRange;
  typedef typename GridLayerType::IndexSet BackendType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
};


} // namespace internal


template <class GridLayerImp, class TupleType>
struct DefaultProductMapperFromTuple;

template <class GridLayerImp, class... MapperTypes>
struct DefaultProductMapperFromTuple<GridLayerImp, std::tuple<MapperTypes...>>
{
  typedef DefaultProductMapper<GridLayerImp, MapperTypes...> type;
};


template <class GridLayerImp, class... MapperTypes>
class DefaultProductMapper
    : public ProductMapperInterface<internal::DefaultProductMapperTraits<GridLayerImp, MapperTypes...>>
{
  typedef ProductMapperInterface<internal::DefaultProductMapperTraits<GridLayerImp, MapperTypes...>> BaseType;

public:
  typedef internal::DefaultProductMapperTraits<GridLayerImp, MapperTypes...> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  static const size_t dimRange = Traits::dimRange;
  using typename BaseType::EntityType;
  using typename BaseType::BackendType;

  DefaultProductMapper(const GridLayerType& grd_layr, const MapperTypes&... mappers)
    : mappers_(std::make_tuple(mappers...))
    , grid_layer_(grd_layr)
  {
  }

  DefaultProductMapper(const GridLayerType& grd_layr, const std::tuple<MapperTypes...>& mappers)
    : mappers_(mappers)
    , grid_layer_(grd_layr)
  {
  }

  template <class... SpaceTypes>
  DefaultProductMapper(const std::tuple<SpaceTypes...>& spaces)
    : DefaultProductMapper(std::get<0>(spaces).grid_layer(),
                           internal::MapperTuplefromSpaceTupleCreator<0, SpaceTypes...>::create(spaces, std::tuple<>()))
  {
  }

  // These methods are required by the ProductMapperInterface
  size_t numDofs(const size_t factor_index, const EntityType& entity) const
  {
    return internal::DynamicTupleGetter<0>::numDofs(mappers_, factor_index, entity);
  }

  size_t maxNumDofs(const size_t factor_index) const
  {
    return internal::DynamicTupleGetter<0>::maxNumDofs(mappers_, factor_index);
  }

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    const size_t num_dofs_entity = numDofs(factor_index, entity);
    if (ret.size() != num_dofs_entity)
      ret.resize(num_dofs_entity);
    for (size_t ii = 0; ii < num_dofs_entity; ++ii)
      ret[ii] = mapToGlobal(factor_index, entity, ii);
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    size_t first_global_index_of_factor = 0;
    for (size_t ii = 0; ii < factor_index; ++ii)
      first_global_index_of_factor += size(ii);
    return first_global_index_of_factor
           + internal::DynamicTupleGetter<0>::mapToGlobal(mappers_, factor_index, entity, local_index_in_factor);
  }

  size_t mapToLocal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    size_t first_local_index_of_factor = 0;
    for (size_t ii = 0; ii < factor_index; ++ii)
      first_local_index_of_factor += numDofs(ii, entity);
    return first_local_index_of_factor + local_index_in_factor;
  }

  const BackendType& backend() const
  {
    return BackendType();
  }

  size_t size(const size_t factor_index) const
  {
    return internal::DynamicTupleGetter<0>::size(mappers_, factor_index);
  }

  size_t size() const
  {
    size_t ret = 0;
    for (size_t ii = 0; ii < sizeof...(MapperTypes); ++ii)
      ret += size(ii);
    return ret;
  }

  size_t numDofs(const EntityType& entity) const
  {
    size_t ret = 0;
    for (size_t ii = 0; ii < sizeof...(MapperTypes); ++ii)
      ret += numDofs(ii, entity);
    return ret;
  }

  size_t maxNumDofs() const
  {
    size_t max_num_dofs = 0;
    const auto it_end = grid_layer_.template end<0>();
    for (auto it = grid_layer_.template begin<0>(); it != it_end; ++it) {
      const auto& entity = *it;
      if (max_num_dofs < numDofs(entity))
        max_num_dofs = numDofs(entity);
    }
    return max_num_dofs;
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    const size_t num_dofs_entity = numDofs(entity);
    if (ret.size() != num_dofs_entity)
      ret.resize(num_dofs_entity);
    for (size_t ii = 0; ii < num_dofs_entity; ++ii)
      ret[ii] = mapToGlobal(entity, ii);
  } // ... globalIndices(...)

  using BaseType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    size_t factor_index = 0;
    size_t local_index_in_factor = localIndex;
    while (local_index_in_factor >= numDofs(factor_index, entity)) {
      local_index_in_factor -= numDofs(factor_index, entity);
      ++factor_index;
    }
    return mapToGlobal(factor_index, entity, local_index_in_factor);
  }

private:
  const std::tuple<MapperTypes...> mappers_;
  const GridLayerType& grid_layer_;
}; // class DefaultProductMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_PRODUCT_HH
