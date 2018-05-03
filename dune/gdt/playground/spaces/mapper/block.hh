// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_MAPPER_BLOCK_HH
#define DUNE_GDT_PLAYGROUND_SPACES_MAPPER_BLOCK_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class LocalSpaceImp>
class BlockMapper;


namespace internal {


template <class LocalSpaceType>
class BlockMapperTraits
{
  static_assert(is_space<LocalSpaceType>::value, "LocalSpaceType has to be derived from SpaceInterface!");

public:
  typedef BlockMapper<LocalSpaceType> derived_type;
  typedef typename LocalSpaceType::EntityType EntityType;
  typedef std::vector<std::shared_ptr<const LocalSpaceType>> BackendType;
}; // class BlockMapperTraits


} // namespace internal


template <class LocalSpaceImp>
class BlockMapper : public MapperInterface<internal::BlockMapperTraits<LocalSpaceImp>>
{
  typedef MapperInterface<internal::BlockMapperTraits<LocalSpaceImp>> BaseType;
  typedef BlockMapper<LocalSpaceImp> ThisType;

public:
  typedef internal::BlockMapperTraits<LocalSpaceImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;
  typedef LocalSpaceImp LocalSpaceType;

  typedef XT::Grid::DD::SubdomainGrid<XT::Grid::extract_grid_t<typename LocalSpaceType::GridLayerType>>
      DdSubdomainsGridType;

private:
  typedef typename DdSubdomainsGridType::GlobalGridViewType GridLayerType;
  typedef typename DdSubdomainsGridType::EntityToSubdomainMapType EntityToSubdomainMapType;

  template <class L, class E>
  class Compute
  {
    static_assert(AlwaysFalse<L>::value, "Not implemented for this kind of entity (only codim 0)!");
  };

  template <class L>
  class Compute<L, typename DdSubdomainsGridType::EntityType>
  {
    typedef typename DdSubdomainsGridType::EntityType Comdim0EntityType;

  public:
    template <int cd, class GridImp, template <int, int, class> class EntityImp>
    static size_t numDofs(const ThisType& self, const Entity<cd, EntityType::dimension, GridImp, EntityImp>& entity)
    {
      const size_t block = find_block_of(self, entity);
      if (self.backend()[block] == nullptr)
        DUNE_THROW(InvalidStateException, "You did not provide a local space for block " << block << "!");
      return self.backend()[block]->mapper().numDofs(entity);
    }

    static void globalIndices(const ThisType& self, const Comdim0EntityType& entity, Dune::DynamicVector<size_t>& ret)
    {
      const size_t block = find_block_of(self, entity);
      if (self.backend()[block] == nullptr)
        DUNE_THROW(InvalidStateException, "You did not provide a local space for block " << block << "!");
      self.backend()[block]->mapper().globalIndices(entity, ret);
      const size_t num_dofs = self.backend()[block]->mapper().numDofs(entity);
      assert(ret.size() >= num_dofs);
      for (size_t ii = 0; ii < num_dofs; ++ii)
        ret[ii] += self.global_start_indices_[block];
    }

    static size_t mapToGlobal(const ThisType& self, const Comdim0EntityType& entity, const size_t& localIndex)
    {
      const size_t block = find_block_of(self, entity);
      if (self.backend()[block] == nullptr)
        DUNE_THROW(InvalidStateException, "You did not provide a local space for block " << block << "!");
      const size_t block_local_index = self.backend()[block]->mapper().mapToGlobal(entity, localIndex);
      return self.global_start_indices_[block] + block_local_index;
    }

  private:
    template <int cd, class GridImp, template <int, int, class> class EntityImp>
    static size_t find_block_of(const ThisType& self,
                                const Entity<cd, EntityType::dimension, GridImp, EntityImp>& entity)
    {
      const auto global_entity_index = self.global_grid_part_->indexSet().index(entity);
      const auto result = self.entity_to_subdomain_map_->find(global_entity_index);
#ifndef NDEBUG
      if (result == self.entity_to_subdomain_map_->end())
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "Entity " << global_entity_index
                             << " of the global grid part was not found in the dd subdomain grid!");
#endif // NDEBUG
      const size_t subdomain = result->second;
#ifndef NDEBUG
      if (subdomain >= self.num_blocks_)
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "The DD subdomains grid is corrupted!\nIt reports Entity " << global_entity_index
                                                                              << " to be in subdomain "
                                                                              << subdomain
                                                                              << " while only having "
                                                                              << self.num_blocks_
                                                                              << " subdomains!");
#endif // NDEBUG
      return subdomain;
    } // ... find_block_of(...)
  }; // class Compute< ..., EntityType >

public:
  BlockMapper(const DdSubdomainsGridType& dd_grid,
              const std::shared_ptr<GridLayerType> grid_layer,
              const std::shared_ptr<std::vector<std::shared_ptr<const LocalSpaceType>>> local_spaces)
    : global_grid_part_(grid_layer)
    , entity_to_subdomain_map_(dd_grid.entityToSubdomainMap())
    , local_spaces_(local_spaces)
    , num_blocks_(local_spaces_->size())
    , size_(0)
    , max_num_dofs_(0)
  {
    if (local_spaces_->size() != dd_grid.size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "You have to provide a local space for each subdomain of the DD subdomains grid!\n"
                     << "  Number of subdomains: "
                     << dd_grid.size()
                     << "\n"
                     << "  Number of local spaces given: "
                     << local_spaces_->size());
    for (size_t bb = 0; bb < num_blocks_; ++bb) {
      if (backend()[bb] != nullptr) {
        max_num_dofs_ = std::max(max_num_dofs_, backend()[bb]->mapper().maxNumDofs());
        global_start_indices_.push_back(size_);
        size_ += backend()[bb]->mapper().size();
      } else {
        global_start_indices_.push_back(size_);
      }
    }
  } // BlockMapper(...)

  BlockMapper(const ThisType& other) = default;
  BlockMapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  size_t numBlocks() const
  {
    return num_blocks_;
  }

  size_t localSize(const size_t block) const
  {
    assert(block < num_blocks_);
    return backend()[block]->mapper().size();
  }

  size_t mapToGlobal(const size_t block, const size_t localIndex) const
  {
    assert(block < num_blocks_);
    return global_start_indices_[block] + localIndex;
  }

  const BackendType& backend() const
  {
    return *local_spaces_;
  }

  size_t size() const
  {
    return size_;
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  template <int cd, class GridImp, template <int, int, class> class EntityImp>
  size_t numDofs(const Entity<cd, EntityType::dimension, GridImp, EntityImp>& entity) const
  {
    return Compute<LocalSpaceType, EntityType>::numDofs(*this, entity);
  }

  using BaseType::globalIndices;

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    Compute<LocalSpaceType, EntityType>::globalIndices(*this, entity, ret);
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    return Compute<LocalSpaceType, EntityType>::mapToGlobal(*this, entity, localIndex);
  }

private:
  template <class L, class E>
  friend class Compute;

  const std::shared_ptr<GridLayerType> global_grid_part_;
  const std::shared_ptr<const typename DdSubdomainsGridType::EntityToSubdomainMapType> entity_to_subdomain_map_;
  const std::shared_ptr<std::vector<std::shared_ptr<const LocalSpaceType>>> local_spaces_;
  size_t num_blocks_;
  size_t size_;
  size_t max_num_dofs_;
  std::vector<size_t> global_start_indices_;
}; // class BlockMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_MAPPER_BLOCK_HH
