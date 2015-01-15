// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_BLOCK_HH
#define DUNE_GDT_MAPPER_BLOCK_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/default.hh>
#endif

#include <dune/gdt/spaces/interface.hh>

#include "../../mapper/interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {

#if HAVE_DUNE_GRID_MULTISCALE


template <class LocalSpaceImp>
class Block;


namespace internal {


template <class LocalSpaceType>
class BlockTraits
{
  static_assert(std::is_base_of<SpaceInterface<typename LocalSpaceType::Traits, LocalSpaceType::dimDomain,
                                               LocalSpaceType::dimRange, LocalSpaceType::dimRangeCols>,
                                LocalSpaceType>::value,
                "LocalSpaceType has to be derived from SpaceInterface!");

public:
  typedef Block<LocalSpaceType> derived_type;
  typedef typename LocalSpaceType::EntityType EntityType;
  typedef typename LocalSpaceType::MapperType::BackendType BackendType;
}; // class BlockTraits


} // namespace internal


template <class LocalSpaceImp>
class Block : public MapperInterface<internal::BlockTraits<LocalSpaceImp>>
{
  typedef MapperInterface<internal::BlockTraits<LocalSpaceImp>> BaseType;

public:
  typedef internal::BlockTraits<LocalSpaceImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;
  typedef LocalSpaceImp LocalSpaceType;

  typedef grid::Multiscale::Default<typename LocalSpaceType::GridViewType::Grid> MsGridType;

private:
  template <class L, class E>
  class Compute
  {
    static_assert(AlwaysFalse<L>::value, "Not implemented for this kind of entity (only codim 0)!");
  };

  template <class L>
  class Compute<L, typename MsGridType::EntityType>
  {
    typedef typename MsGridType::EntityType Comdim0EntityType;

  public:
    static size_t numDofs(const MsGridType& ms_grid, const std::vector<std::shared_ptr<const L>>& local_spaces,
                          const Comdim0EntityType& entity)
    {
      const size_t block = find_block_of_(ms_grid, entity);
      return local_spaces[block]->mapper().numDofs(entity);
    }

    static void globalIndices(const MsGridType& ms_grid, const std::vector<std::shared_ptr<const L>>& local_spaces,
                              const std::vector<size_t>& global_start_indices, const Comdim0EntityType& entity,
                              Dune::DynamicVector<size_t>& ret)
    {
      const size_t block = find_block_of_(ms_grid, entity);
      local_spaces[block]->mapper().globalIndices(entity, ret);
      const size_t num_dofs = local_spaces[block]->mapper().numDofs(entity);
      assert(ret.size() >= num_dofs);
      for (size_t ii = 0; ii < num_dofs; ++ii)
        ret[ii] += global_start_indices[block];
    }

    static size_t mapToGlobal(const MsGridType& ms_grid, const std::vector<std::shared_ptr<const L>>& local_spaces,
                              const std::vector<size_t>& global_start_indices, const Comdim0EntityType& entity,
                              const size_t& localIndex)
    {
      const size_t block             = find_block_of_(ms_grid, entity);
      const size_t block_local_index = local_spaces[block]->mapper().mapToGlobal(entity, localIndex);
      return global_start_indices[block] + block_local_index;
    }

  private:
    static size_t find_block_of_(const MsGridType& ms_grid, const Comdim0EntityType& entity)
    {
      const auto global_entity_index = ms_grid.globalGridView().indexSet().index(entity);
      const auto result              = ms_grid.entityToSubdomainMap()->find(global_entity_index);
#ifndef NDEBUG
      if (result == ms_grid.entityToSubdomainMap()->end())
        DUNE_THROW(Stuff::Exceptions::internal_error,
                   "Entity " << global_entity_index
                             << " of the global grid view was not found in the multiscale grid!");
#endif // NDEBUG
      const size_t subdomain = result->second;
#ifndef NDEBUG
      if (subdomain >= ms_grid.size())
        DUNE_THROW(Stuff::Exceptions::internal_error,
                   "The multiscale grid is corrupted!\nIt reports Entity " << global_entity_index
                                                                           << " to be in subdomain "
                                                                           << subdomain
                                                                           << " while only having "
                                                                           << ms_grid.size()
                                                                           << " subdomains!");
#endif // NDEBUG
      return subdomain;
    } // ... find_block_of_(...)
  }; // class Compute< ..., EntityType >

public:
  Block(const std::shared_ptr<const MsGridType> ms_grid,
        const std::vector<std::shared_ptr<const LocalSpaceType>> local_spaces)
    : ms_grid_(ms_grid)
    , local_spaces_(local_spaces)
    , num_blocks_(local_spaces_.size())
    , size_(0)
    , max_num_dofs_(0)
  {
    if (local_spaces_.size() != ms_grid_->size())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                 "You have to provide a local space for each subdomain of the multiscale grid!\n"
                     << "  Size of the given multiscale grid: "
                     << ms_grid_->size()
                     << "\n"
                     << "  Number of local spaces given: "
                     << local_spaces_.size());
    for (size_t bb = 0; bb < num_blocks_; ++bb) {
      max_num_dofs_ = std::max(max_num_dofs_, local_spaces_[bb]->mapper().maxNumDofs());
      global_start_indices_.push_back(size_);
      size_ += local_spaces_[bb]->mapper().size();
    }
  } // Block(...)

  size_t numBlocks() const
  {
    return num_blocks_;
  }

  size_t localSize(const size_t block) const
  {
    assert(block < num_blocks_);
    return local_spaces_[block]->mapper().size();
  }

  size_t mapToGlobal(const size_t block, const size_t localIndex) const
  {
    assert(block < num_blocks_);
    return global_start_indices_[block] + localIndex;
  }

  const BackendType& backend() const
  {
    return local_spaces_[0]->mapper().backend();
  }

  size_t size() const
  {
    return size_;
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  size_t numDofs(const EntityType& entity) const
  {
    return Compute<LocalSpaceType, EntityType>::numDofs(*ms_grid_, local_spaces_, entity);
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    Compute<LocalSpaceType, EntityType>::globalIndices(*ms_grid_, local_spaces_, global_start_indices_, entity, ret);
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    return Compute<LocalSpaceType, EntityType>::mapToGlobal(
        *ms_grid_, local_spaces_, global_start_indices_, entity, localIndex);
  } // ... mapToGlobal(...)

private:
  std::shared_ptr<const MsGridType> ms_grid_;
  std::vector<std::shared_ptr<const LocalSpaceType>> local_spaces_;
  size_t num_blocks_;
  size_t size_;
  size_t max_num_dofs_;
  std::vector<size_t> global_start_indices_;
}; // class Block


#else // HAVE_DUNE_GRID_MULTISCALE


template <class LocalSpaceImp>
class Block
{
  static_assert(AlwaysFalse<LocalSpaceImp>::value, "You are missing dune-grid-multiscale!");
};


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_BLOCK_HH
