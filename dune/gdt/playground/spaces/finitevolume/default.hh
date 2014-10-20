// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_FINITEVOLUME_HH
#define DUNE_GDT_SPACE_FINITEVOLUME_HH

#include <dune/gdt/spaces/parallel.hh>

#include "../../mapper/finitevolume.hh"
#include "../../basefunctionset/finitevolume.hh"
#include "../../../spaces/interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace FiniteVolume {


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class Default
{
  static_assert(Dune::AlwaysFalse< GridViewImp >::value, "Untested for these dimensions!");
};


/**
 *  \brief Traits class for Spaces::ContinuousLagrange::FemBased.
 */
template< class GridViewImp, class RangeFieldImp, int rangeDim, int rangeDimCols >
class DefaultTraits
{
public:
  typedef Default< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  static const int  polOrder = 0;
  typedef GridViewImp                     GridViewType;
  typedef typename GridViewType::IndexSet BackendType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef RangeFieldImp     RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Mapper::FiniteVolume< GridViewType, dimRange, dimRangeCols > MapperType;
  typedef BaseFunctionSet::FiniteVolume< typename GridViewType::template Codim< 0 >::Entity
                                       , typename GridViewType::ctype, GridViewType::dimension
                                       , RangeFieldType, dimRange, dimRangeCols > BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool                        needs_grid_view = true;
  typedef CommunicationChooser< GridViewType >    CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class DefaultTraits


template< class GridViewImp, class RangeFieldImp, int rangeDim >
class Default< GridViewImp, RangeFieldImp, rangeDim, 1 >
  : public SpaceInterface< DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, 1 > >
{
  typedef SpaceInterface< DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, 1 > > BaseType;
  typedef Default< GridViewImp, RangeFieldImp, rangeDim, 1 >                         ThisType;
public:
  typedef DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, 1 > Traits;

  typedef typename Traits::GridViewType GridViewType;
  static const int                      polOrder = Traits::polOrder;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridViewType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
  typedef typename Traits::CommunicatorType         CommunicatorType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  Default(GridViewType gv)
    : grid_view_(gv)
    , mapper_(grid_view_)
    , communicator_(CommunicationChooserType::create(grid_view_))
  {}

  Default(const ThisType& other)
    : grid_view_(other.grid_view_)
    , mapper_(other.mapper_)
    , communicator_(CommunicationChooserType::create(grid_view_))
  {}

  Default(ThisType&& source) = default;

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

  using BaseType::compute_pattern;

  template< class G, class S >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S >& ansatz_space) const
  {
    return BaseType::compute_face_and_volume_pattern(local_grid_view, ansatz_space);
  }

  CommunicatorType& communicator() const
  {
    // no need to prepare the communicator, since we are not pdelab based
    return *communicator_;
  }

private:
  const GridViewType grid_view_;
  const MapperType mapper_;
  const std::unique_ptr< CommunicatorType > communicator_;
}; // class Default< ..., 1, 1 >


} // namespace FiniteVolume
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_FINITEVOLUME_HH
