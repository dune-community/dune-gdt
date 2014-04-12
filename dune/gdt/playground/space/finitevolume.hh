// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_FINITEVOLUME_HH
#define DUNE_GDT_SPACE_FINITEVOLUME_HH

#include "../mapper/finitevolume.hh"
#include "../basefunctionset/finitevolume.hh"
#include "../../space/interface.hh"

namespace Dune {
namespace GDT {
namespace FiniteVolumeSpace {


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Default
{
  static_assert(Dune::AlwaysFalse<GridViewImp>::value, "Untested for these dimensions!");
};


/**
 *  \brief Traits class for ContinuousLagrangeSpace::FemWrapper.
 */
template <class GridViewImp, class RangeFieldImp, int rangeDim, int rangeDimCols>
class DefaultTraits
{
public:
  typedef Default<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  static const int polOrder = 0;
  typedef double BackendType;
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Mapper::FiniteVolume<GridViewType, dimRange, dimRangeCols> MapperType;
  typedef BaseFunctionSet::FiniteVolume<typename GridViewType::template Codim<0>::Entity, typename GridViewType::ctype,
                                        GridViewType::dimension, RangeFieldType, dimRange,
                                        dimRangeCols> BaseFunctionSetType;
  static const bool needs_grid_view = true;
}; // class DefaultTraits


template <class GridViewImp, class RangeFieldImp, int rangeDim>
class Default<GridViewImp, RangeFieldImp, rangeDim, 1>
    : public SpaceInterface<DefaultTraits<GridViewImp, RangeFieldImp, rangeDim, 1>>
{
  typedef SpaceInterface<DefaultTraits<GridViewImp, RangeFieldImp, rangeDim, 1>> BaseType;
  typedef Default<GridViewImp, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef DefaultTraits<GridViewImp, RangeFieldImp, rangeDim, 1> Traits;

  typedef typename Traits::GridViewType GridViewType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = Traits::dimRange;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  Default(const std::shared_ptr<const GridViewType>& gv)
    : grid_view_(gv)
    , mapper_(std::make_shared<MapperType>(*grid_view_))
    , backend_(1)
  {
  }

  Default(const ThisType& other)
    : grid_view_(other.grid_view_)
    , mapper_(other.mapper_)
  {
  }

  Default& operator=(const ThisType& other)
  {
    if (this != &other) {
      grid_view_ = other.grid_view_;
      mapper_    = other.mapper_;
    }
    return *this;
  }

  ~Default()
  {
  }

  const std::shared_ptr<const GridViewType>& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(entity);
  }

  using BaseType::compute_pattern;

  template <class G, class S>
  PatternType compute_pattern(const GridView<G>& local_grid_view, const SpaceInterface<S>& ansatz_space) const
  {
    return BaseType::compute_face_and_volume_pattern(local_grid_view, ansatz_space);
  }

private:
  std::shared_ptr<const GridViewType> grid_view_;
  std::shared_ptr<const MapperType> mapper_;
  const BackendType backend_;
}; // class Default< ..., 1, 1 >


} // namespace FiniteVolumeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_FINITEVOLUME_HH
