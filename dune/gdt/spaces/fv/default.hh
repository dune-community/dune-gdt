// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_FV_SPACE_HH
#define DUNE_GDT_SPACES_FV_SPACE_HH

#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/spaces/basefunctionset/fv.hh>
#include <dune/gdt/spaces/mapper/fv.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FvSpace
{
  static_assert(Dune::AlwaysFalse<GridViewImp>::value, "Untested for these dimensions!");
};


namespace internal {


/**
 *  \brief Traits class for FvSpace.
 */
template <class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FvSpaceTraits
{
public:
  typedef FvSpace<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  static const int polOrder    = 0;
  static const bool continuous = false;
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::IndexSet BackendType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef RangeFieldImp RangeFieldType;
  typedef FvMapper<GridViewType, rangeDim, rangeDimCols> MapperType;
  typedef BaseFunctionSet::FiniteVolume<typename GridViewType::template Codim<0>::Entity, typename GridViewType::ctype,
                                        GridViewType::dimension, RangeFieldType, rangeDim,
                                        rangeDimCols> BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view                       = true;
  typedef CommunicationChooser<GridViewType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class FvSpaceTraits


} // namespace internal


template <class GridViewImp, class RangeFieldImp, size_t rangeDim>
class FvSpace<GridViewImp, RangeFieldImp, rangeDim, 1>
    : public FvSpaceInterface<internal::FvSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1>, GridViewImp::dimension,
                              rangeDim, 1>
{
  typedef FvSpace<GridViewImp, RangeFieldImp, rangeDim, 1> ThisType;
  typedef FvSpaceInterface<internal::FvSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1>, GridViewImp::dimension,
                           rangeDim, 1> BaseType;

public:
  typedef typename internal::FvSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1> Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;

public:
  using typename BaseType::CommunicatorType;

  FvSpace(GridViewType gv)
    : grid_view_(gv)
    , mapper_(grid_view_)
    , communicator_(CommunicationChooserType::create(grid_view_))
  {
  }

  FvSpace(const ThisType& other)
    : grid_view_(other.grid_view_)
    , mapper_(other.mapper_)
    , communicator_(CommunicationChooserType::create(grid_view_))
  {
  }

  FvSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return grid_view_.indexSet();
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
  const GridViewType grid_view_;
  const MapperType mapper_;
  const std::unique_ptr<CommunicatorType> communicator_;
}; // class FvSpace< ..., 1, 1 >


template <class R, size_t r, size_t rC, class GV>
FvSpace<GV, R, r, rC> make_fv_space(const GV& grid_view)
{
  return FvSpace<GV, R, r, rC>(grid_view);
}

template <class R, size_t r, class GV>
FvSpace<GV, R, r, 1> make_fv_space(const GV& grid_view)
{
  return FvSpace<GV, R, r, 1>(grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_SPACE_HH
