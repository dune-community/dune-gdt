// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BLOCK_HH
#define DUNE_GDT_SPACES_BLOCK_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>

#include "mapper/block.hh"

#include "../../spaces/interface.hh"

namespace Dune {
namespace GDT {


template <class LocalSpaceImp>
class BlockSpace;


namespace internal {


template <class LocalSpaceType>
class BlockSpaceTraits
{
  static_assert(is_space<LocalSpaceType>::value, "LocalSpaceType has to be derived from SpaceInterface!");
  typedef XT::Grid::DD::SubdomainGrid<typename LocalSpaceType::GridLayerType::Grid> DdSubdomainsGridType;

public:
  typedef BlockSpace<LocalSpaceType> derived_type;
  static const int polOrder = LocalSpaceType::polOrder;
  static const bool continuous = false;
  typedef std::vector<LocalSpaceType> BackendType;
  typedef BlockMapper<LocalSpaceType> MapperType;
  typedef typename LocalSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename LocalSpaceType::CommunicatorType CommunicatorType;
  typedef typename DdSubdomainsGridType::GlobalGridPartType GridLayerType;
  typedef typename LocalSpaceType::RangeFieldType RangeFieldType;

  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::part;

  static const bool needs_grid_view = false;
}; // class BlockSpaceTraits


} // namespace internal


template <class LocalSpaceImp>
class BlockSpace : public SpaceInterface<internal::BlockSpaceTraits<LocalSpaceImp>,
                                         LocalSpaceImp::dimDomain,
                                         LocalSpaceImp::dimRange,
                                         LocalSpaceImp::dimRangeCols>
{
  typedef SpaceInterface<internal::BlockSpaceTraits<LocalSpaceImp>,
                         LocalSpaceImp::dimDomain,
                         LocalSpaceImp::dimRange,
                         LocalSpaceImp::dimRangeCols>
      BaseType;
  typedef BlockSpace<LocalSpaceImp> ThisType;

public:
  typedef internal::BlockSpaceTraits<LocalSpaceImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef LocalSpaceImp LocalSpaceType;

  typedef typename BaseType::PatternType PatternType;
  typedef typename BaseType::GridLayerType GridLayerType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::CommunicatorType CommunicatorType;

  typedef XT::Grid::DD::SubdomainGrid<typename XT::Grid::extract_grid<GridLayerType>::type> DdSubdomainsGridType;

  BlockSpace(const DdSubdomainsGridType& dd_grid, const std::shared_ptr<std::vector<LocalSpaceType>> local_spaces)
    : entity_to_subdomain_map_(dd_grid.entityToSubdomainMap())
    , global_grid_part_(new GridLayerType(dd_grid.globalGridPart()))
    , local_spaces_(local_spaces)
    , mapper_(new MapperType(dd_grid, global_grid_part_, local_spaces_))
  {
    if (local_spaces_->size() != dd_grid.size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "You have to provide a local space for each subdomain of the DD subdomains grid!\n"
                     << "  Number of subdomains: "
                     << dd_grid.size()
                     << "\n"
                     << "  Number of local spaces given: "
                     << local_spaces_->size());
  } // BlockSpace(...)

  BlockSpace(const ThisType& other) = default;
  BlockSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return *global_grid_part_;
  }

  GridLayerType& grid_layer()
  {
    return *global_grid_part_;
  }

  const BackendType& backend() const
  {
    return *local_spaces_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    const size_t block = find_block_of(entity);
    return backend()[block].base_function_set(entity);
  }

  template <class ConstraintsType>
  void local_constraints(const EntityType& /*entity*/, ConstraintsType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
  }

  template <class G, class S, size_t d, size_t r, size_t rC>
  PatternType compute_pattern(const GridView<G>& /*local_grid_layer*/,
                              const SpaceInterface<S, d, r, rC>& /*ansatz_space*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
    return PatternType();
  }

  CommunicatorType& communicator() const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this, I probably need my own communicator!");
    return backend()[0].communicator();
  }

private:
  template <class EntityType>
  size_t find_block_of(const EntityType& entity) const
  {
    const auto global_entity_index = global_grid_part_->indexSet().index(entity);
    const auto result = entity_to_subdomain_map_->find(global_entity_index);
#ifndef NDEBUG
    if (result == entity_to_subdomain_map_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "Entity " << global_entity_index << " of the global grid part was not found in the multiscale grid!");
#endif // NDEBUG
    const size_t subdomain = result->second;
#ifndef NDEBUG
    if (subdomain >= local_spaces_->size())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "The DD subdomains grid is corrupted!\nIt reports Entity " << global_entity_index
                                                                            << " to be in subdomain "
                                                                            << subdomain
                                                                            << " while only having "
                                                                            << local_spaces_->size()
                                                                            << " subdomains!");
#endif // NDEBUG
    return subdomain;
  } // ... find_block_of(...)

  const std::shared_ptr<const typename DdSubdomainsGridType::EntityToSubdomainMapType> entity_to_subdomain_map_;
  const std::shared_ptr<GridLayerType> global_grid_part_;
  const std::shared_ptr<std::vector<LocalSpaceType>> local_spaces_;
  const std::shared_ptr<MapperType> mapper_;
}; // class BlockSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BLOCK_HH
