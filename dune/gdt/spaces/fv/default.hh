// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_FV_SPACE_HH
#define DUNE_GDT_SPACES_FV_SPACE_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basis/finite-volume.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>
#include <dune/gdt/spaces/parallel.hh>
#include <dune/gdt/spaces/interface.hh>


namespace Dune {
namespace GDT {


// forward, to be used in the traits and to allow for specialization
template <class GL, class R, size_t r, size_t rC = 1>
class FvSpace
{
  static_assert(Dune::AlwaysFalse<GL>::value, "Untested for these dimensions!");
};


namespace internal {


/**
 *  \brief Traits class for FvSpace.
 */
template <class GL, class R, size_t r, size_t rC>
class FvSpaceTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");

public:
  typedef FvSpace<GL, R, r, rC> derived_type;
  static const int polOrder = 0;
  static const bool continuous = false;
  typedef GL GridLayerType;
  typedef typename GridLayerType::IndexSet BackendType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef R RangeFieldType;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const bool needs_grid_view = true;
  typedef DofCommunicationChooser<GridLayerType> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;
  static const constexpr Backends backend_type{GDT::Backends::gdt};
}; // class FvSpaceTraits


} // namespace internal


template <class GL, class R>
class FvSpace<GL, R, 1, 1> : public SpaceInterface<internal::FvSpaceTraits<GL, R, 1, 1>, GL::dimension, 1, 1>
{
  typedef FvSpace<GL, R, 1, 1> ThisType;
  typedef SpaceInterface<internal::FvSpaceTraits<GL, R, 1, 1>, GL::dimension, 1, 1> BaseType;

public:
  typedef typename internal::FvSpaceTraits<GL, R, 1, 1> Traits;
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::GlobalBasisType;

private:
  typedef typename Traits::DofCommunicationChooserType DofCommunicationChooserType;
  using MapperImplementation = FiniteVolumeMapper<GridLayerType, 1, 1>;
  using GlobalBasisImplementation = FiniteVolumeGlobalBasis<GridLayerType, R>;

public:
  using typename BaseType::DofCommunicatorType;

  FvSpace(GridLayerType grd_layr)
    : grid_layer_(grd_layr)
    , mapper_(grid_layer_)
    , basis_(grid_layer_)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
  {
  }

  FvSpace(const ThisType& other)
    : grid_layer_(other.grid_layer_)
    , mapper_(other.mapper_)
    , basis_(other.basis_)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
  {
  }

  FvSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    return grid_layer_.indexSet();
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  const GlobalBasisType& basis() const
  {
    return basis_;
  }

  DofCommunicatorType& dof_communicator() const
  {
    // no need to prepare the communicator
    return *communicator_;
  }

private:
  GridLayerType grid_layer_;
  const MapperImplementation mapper_;
  const GlobalBasisImplementation basis_;
  const std::unique_ptr<DofCommunicatorType> communicator_;
}; // class FvSpace< ..., 1, 1 >


template <class R, size_t r, size_t rC, class GL>
FvSpace<GL, R, r, rC> make_fv_space(const GL& grid_layer)
{
  return FvSpace<GL, R, r, rC>(grid_layer);
}

template <class R, size_t r, class GL>
FvSpace<GL, R, r, 1> make_fv_space(const GL& grid_layer)
{
  return FvSpace<GL, R, r, 1>(grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_SPACE_HH
