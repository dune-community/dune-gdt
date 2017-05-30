// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
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

#include <dune/gdt/spaces/basefunctionset/fv.hh>
#include <dune/gdt/spaces/mapper/fv.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits and to allow for specialization
template <class GridLayerImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FvSpace
{
  static_assert(Dune::AlwaysFalse<GridLayerImp>::value, "Untested for these dimensions!");
};


namespace internal {


/**
 *  \brief Traits class for FvSpace.
 */
template <class GridLayerImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FvSpaceTraits
{
  static_assert(XT::Grid::is_layer<GridLayerImp>::value, "");

public:
  typedef FvSpace<GridLayerImp, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  static const int polOrder = 0;
  static const bool continuous = false;
  typedef GridLayerImp GridLayerType;
  typedef typename GridLayerType::IndexSet BackendType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef RangeFieldImp RangeFieldType;
  typedef FvMapper<GridLayerType, rangeDim, rangeDimCols> MapperType;
  typedef BaseFunctionSet::FiniteVolume<EntityType,
                                        typename GridLayerType::ctype,
                                        GridLayerType::dimension,
                                        RangeFieldType,
                                        rangeDim,
                                        rangeDimCols>
      BaseFunctionSetType;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const bool needs_grid_view = true;
  typedef CommunicationChooser<GridLayerType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
  static const constexpr auto backend_type{GDT::Backends::gdt};
}; // class FvSpaceTraits


} // namespace internal


template <class GridLayerImp, class RangeFieldImp, size_t rangeDim>
class FvSpace<GridLayerImp, RangeFieldImp, rangeDim, 1>
    : public FvSpaceInterface<internal::FvSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1>,
                              GridLayerImp::dimension,
                              rangeDim,
                              1>
{
  typedef FvSpace<GridLayerImp, RangeFieldImp, rangeDim, 1> ThisType;
  typedef FvSpaceInterface<internal::FvSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1>,
                           GridLayerImp::dimension,
                           rangeDim,
                           1>
      BaseType;

public:
  typedef typename internal::FvSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1> Traits;
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;

public:
  using typename BaseType::CommunicatorType;

  FvSpace(GridLayerType grd_layr)
    : grid_layer_(grd_layr)
    , mapper_(grid_layer_)
    , communicator_(CommunicationChooserType::create(grid_layer_))
  {
  }

  FvSpace(const ThisType& other)
    : grid_layer_(other.grid_layer_)
    , mapper_(other.mapper_)
    , communicator_(CommunicationChooserType::create(grid_layer_))
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

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(entity);
  }

  CommunicatorType& communicator() const
  {
    // no need to prepare the communicator, since we are not pdelab based
    return *communicator_;
  }

private:
  GridLayerType grid_layer_;
  const MapperType mapper_;
  const std::unique_ptr<CommunicatorType> communicator_;
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
