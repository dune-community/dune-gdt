// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BLOCK_HH
#define DUNE_GDT_SPACES_BLOCK_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>

#if HAVE_DUNE_GRID_MULTISCALE
#include <dune/grid/multiscale/default.hh>
#endif

#include "mapper/block.hh"

#include "../../spaces/interface.hh"

namespace Dune {
namespace GDT {

#if HAVE_DUNE_GRID_MULTISCALE


template <class LocalSpaceImp>
class BlockSpace;


namespace internal {


template <class LocalSpaceType>
class BlockSpaceTraits
{
  static_assert(std::is_base_of<SpaceInterface<typename LocalSpaceType::Traits, LocalSpaceType::dimDomain,
                                               LocalSpaceType::dimRange, LocalSpaceType::dimRangeCols>,
                                LocalSpaceType>::value,
                "LocalSpaceType has to be derived from SpaceInterface!");
  typedef grid::Multiscale::Default<typename LocalSpaceType::GridViewType::Grid> MsGridType;

public:
  typedef BlockSpace<LocalSpaceType> derived_type;
  static const int polOrder    = LocalSpaceType::polOrder;
  static const bool continuous = false;
  typedef typename LocalSpaceType::BackendType BackendType;
  typedef BlockMapper<LocalSpaceType> MapperType;
  typedef typename LocalSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename LocalSpaceType::CommunicatorType CommunicatorType;
  typedef typename MsGridType::GlobalGridViewType GridViewType;
  typedef typename LocalSpaceType::RangeFieldType RangeFieldType;

  static const XT::Grid::Backends part_view_type = LocalSpaceType::part_view_type;

  static const bool needs_grid_view = LocalSpaceType::needs_grid_view;
}; // class BlockSpaceTraits


} // namespace internal


/**
 * \todo This can be implemented easier by now. Since all local spaces hold a copy of their local grid view it should
 *       be enough to hold a copy of the global grid view in this space (if the ms_grid is not needed elsewhere)
 */
template <class LocalSpaceImp>
class BlockSpace : public SpaceInterface<internal::BlockSpaceTraits<LocalSpaceImp>, LocalSpaceImp::dimDomain,
                                         LocalSpaceImp::dimRange, LocalSpaceImp::dimRangeCols>
{
  typedef SpaceInterface<internal::BlockSpaceTraits<LocalSpaceImp>, LocalSpaceImp::dimDomain, LocalSpaceImp::dimRange,
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
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::CommunicatorType CommunicatorType;

  typedef grid::Multiscale::Default<typename GridViewType::Grid> MsGridType;

  BlockSpace(const std::shared_ptr<const MsGridType>& ms_grid,
             const std::vector<std::shared_ptr<const LocalSpaceType>>& local_spaces)
    : ms_grid_(ms_grid)
    , grid_view_(ms_grid_->globalGridView())
    , local_spaces_(local_spaces)
    , mapper_(std::make_shared<MapperType>(ms_grid_, local_spaces_))
  {
    if (local_spaces_.size() != ms_grid_->size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "You have to provide a local space for each subdomain of the multiscale grid!\n"
                     << "  Size of the given multiscale grid: "
                     << ms_grid_->size()
                     << "\n"
                     << "  Number of local spaces given: "
                     << local_spaces_.size());
  } // BlockSpace(...)

  BlockSpace(const ThisType& other) = default;

  BlockSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const std::shared_ptr<const MsGridType>& ms_grid() const
  {
    return ms_grid_;
  }

  const std::vector<std::shared_ptr<const LocalSpaceType>>& local_spaces() const
  {
    return local_spaces_;
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return local_spaces_[0]->backend();
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    const size_t block = find_block_of_(entity);
    return local_spaces_[block]->base_function_set(entity);
  }

  template <class ConstraintsType>
  void local_constraints(const EntityType& /*entity*/, ConstraintsType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
  }

  template <class G, class S, size_t d, size_t r, size_t rC>
  PatternType compute_pattern(const GridView<G>& /*local_grid_view*/,
                              const SpaceInterface<S, d, r, rC>& /*ansatz_space*/) const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
    return PatternType();
  }

  CommunicatorType& communicator() const
  {
    DUNE_THROW(NotImplemented, "I am not sure yet how to implement this!");
    return local_spaces_[0]->communicator();
  }

private:
  template <class EntityType>
  size_t find_block_of_(const EntityType& entity) const
  {
    const auto global_entity_index = ms_grid_->globalGridView().indexSet().index(entity);
    const auto result              = ms_grid_->entityToSubdomainMap()->find(global_entity_index);
#ifndef NDEBUG
    if (result == ms_grid_->entityToSubdomainMap()->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "Entity " << global_entity_index << " of the global grid view was not found in the multiscale grid!");
#endif // NDEBUG
    const size_t subdomain = result->second;
#ifndef NDEBUG
    if (subdomain >= ms_grid_->size())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "The multiscale grid is corrupted!\nIt reports Entity " << global_entity_index
                                                                         << " to be in subdomain "
                                                                         << subdomain
                                                                         << " while only having "
                                                                         << ms_grid_->size()
                                                                         << " subdomains!");
#endif // NDEBUG
    assert(subdomain < local_spaces_.size());
    return subdomain;
  } // ... find_block_of_(...)

  const std::shared_ptr<const MsGridType> ms_grid_;
  const GridViewType grid_view_;
  const std::vector<std::shared_ptr<const LocalSpaceType>> local_spaces_;
  const std::shared_ptr<const MapperType> mapper_;
}; // class Block


#else // HAVE_DUNE_GRID_MULTISCALE


template <class LocalSpaceImp>
class BlockSpace
{
  static_assert(Dune::AlwaysFalse<LocalSpaceImp>::value, "You are missing dune-grid-multiscale!");
};


#endif // HAVE_DUNE_GRID_MULTISCALE

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BLOCK_HH
