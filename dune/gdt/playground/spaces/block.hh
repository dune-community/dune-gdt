// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BLOCK_HH
#define DUNE_GDT_SPACES_BLOCK_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/playground/spaces/mapper/block.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/parallel.hh>

namespace Dune {
namespace GDT {


template <class LocalSpaceImp>
class BlockSpace;


namespace internal {


template <class LocalSpaceType>
class BlockSpaceTraits
{
  static_assert(is_space<LocalSpaceType>::value, "LocalSpaceType has to be derived from SpaceInterface!");
  typedef XT::Grid::DD::SubdomainGrid<XT::Grid::extract_grid_t<typename LocalSpaceType::GridLayerType>>
      DdSubdomainsGridType;

public:
  typedef BlockSpace<LocalSpaceType> derived_type;
  static const int polOrder = LocalSpaceType::polOrder;
  static const bool continuous = false;
  typedef std::vector<std::shared_ptr<const LocalSpaceType>> BackendType;
  typedef BlockMapper<LocalSpaceType> MapperType;
  typedef typename LocalSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename DdSubdomainsGridType::GlobalGridViewType GridLayerType;
  typedef typename LocalSpaceType::RangeFieldType RangeFieldType;

  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const constexpr Backends backend_type{Backends::gdt};
  using DofCommunicationChooserType = DofCommunicationChooser<GridLayerType, true>;
  using DofCommunicatorType = typename DofCommunicationChooserType::Type;

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
  typedef DofCommunicationChooser<GridLayerType> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;

  typedef XT::Grid::DD::SubdomainGrid<typename XT::Grid::extract_grid<GridLayerType>::type> DdSubdomainsGridType;
  static const constexpr Backends backend_type{Backends::gdt};

  BlockSpace(const DdSubdomainsGridType& grid, std::vector<std::shared_ptr<const LocalSpaceType>> spaces)
    : dd_grid_(grid)
    , entity_to_subdomain_map_(dd_grid_.entityToSubdomainMap())
    , global_grid_part_(new GridLayerType(dd_grid_.global_grid_view()))
    , local_spaces_(new std::vector<std::shared_ptr<const LocalSpaceType>>(spaces))
    , mapper_(new MapperType(dd_grid_, global_grid_part_, local_spaces_))
    , communicator_(Traits::DofCommunicationChooserType::create(*global_grid_part_))
    , communicator_prepared_(false)
  {
    if (local_spaces_->size() != dd_grid_.size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "You have to provide a local space (or a nullptr) for each subdomain of the DD subdomains grid!\n"
                     << "  Number of subdomains: "
                     << dd_grid_.size()
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
    if (backend()[block] == nullptr)
      DUNE_THROW(InvalidStateException, "You did not provide a local space for block " << block << "!");
    return backend()[block]->base_function_set(entity);
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

  typename Traits::DofCommunicatorType& dof_communicator() const
  {
    if (!communicator_prepared_)
      communicator_prepared_ = DofCommunicationChooserType::prepare(*this, *communicator_);
    return *communicator_;
    return *communicator_;
  }

  const DdSubdomainsGridType& dd_grid() const
  {
    return dd_grid_;
  }

  size_t num_blocks() const
  {
    return local_spaces_->size();
  }

  const LocalSpaceType& local_space(const size_t block) const
  {
    if (block >= num_blocks())
      DUNE_THROW(XT::Common::Exceptions::index_out_of_range,
                 "  num_blocks: " << num_blocks() << "\n  block: " << block);
    return *(local_spaces_->at(block));
  }

  static constexpr bool associates_data_with(int codim)
  {
    return LocalSpaceType::associates_data_with(codim);
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
    if (subdomain >= dd_grid_.size())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "The DD subdomains grid is corrupted!\nIt reports Entity " << global_entity_index
                                                                            << " to be in subdomain "
                                                                            << subdomain
                                                                            << " while only having "
                                                                            << dd_grid_.size()
                                                                            << " subdomains!");
#endif // NDEBUG
    return subdomain;
  } // ... find_block_of(...)

  const DdSubdomainsGridType& dd_grid_;
  const std::shared_ptr<const typename DdSubdomainsGridType::EntityToSubdomainMapType> entity_to_subdomain_map_;
  const std::shared_ptr<GridLayerType> global_grid_part_;
  const std::shared_ptr<std::vector<std::shared_ptr<const LocalSpaceType>>> local_spaces_;
  const std::shared_ptr<MapperType> mapper_;
  mutable std::shared_ptr<typename Traits::DofCommunicatorType> communicator_;
  mutable bool communicator_prepared_;
}; // class BlockSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BLOCK_HH
