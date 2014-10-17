// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH

#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>

#if HAVE_DUNE_FEM
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/fem/space/common/functionspace.hh>
#   include <dune/fem/space/lagrange/space.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_DUNE_FEM

#include <dune/gdt/spaces/parallel.hh>

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem.hh"

#include "base.hh"
#include "../constraints.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace ContinuousLagrange {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemBased
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "Untested for these dimensions!");
};


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols >
class FemBasedTraits
{
public:
  typedef FemBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int                            polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");
  static const unsigned int             dimDomain = GridPartType::dimension;
private:
  typedef typename GridPartType::ctype  DomainFieldType;
public:
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
  static const unsigned int             dimRangeCols = rangeDimCols;
private:
  typedef Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
public:
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > BackendType;
  typedef Mapper::FemDofWrapper< typename BackendType::BlockMapperType, 1 > MapperType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef BaseFunctionSet::FemWrapper
      < typename BackendType::ShapeFunctionSetType, EntityType, DomainFieldType, dimDomain,
        RangeFieldType, dimRange, dimRangeCols > BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::part;
  static const bool needs_grid_view = false;
  typedef CommunicationChooser< GridViewType, false > CommunicationChooserType;
  typedef typename CommunicationChooserType::Type     CommunicatorType;
}; // class SpaceWrappedFemContinuousLagrangeTraits


// untested for the vector-valued case, especially Spaces::ContinuousLagrangeBase
template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemBased< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
  : public Spaces::ContinuousLagrangeBase< FemBasedTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
                                      , GridPartImp::dimension, RangeFieldImp, 1, 1 >
{
  typedef Spaces::ContinuousLagrangeBase< FemBasedTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
                                     , GridPartImp::dimension, RangeFieldImp, 1, 1 >  BaseType;
  typedef FemBased< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >               ThisType;

public:
  typedef FemBasedTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >         Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::GridViewType GridViewType;
  static const int                      polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType              BackendType;
  typedef typename Traits::MapperType               MapperType;
  typedef typename Traits::BaseFunctionSetType      BaseFunctionSetType;
  typedef typename Traits::EntityType               EntityType;
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
  typedef typename Traits::CommunicatorType         CommunicatorType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  FemBased(GridPartType gridP)
    : gridPart_(gridP)
    , gridView_(gridPart_.gridView())
    , backend_(std::make_shared< BackendType >(const_cast< GridPartType& >(gridPart_)))
    , mapper_(std::make_shared< MapperType >(backend_->blockMapper()))
    , communicator_(CommunicationChooserType::create(gridView_))
  {}

  FemBased(const ThisType& other) = default;

  FemBased(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridPartType& grid_part() const
  {
    return gridPart_;
  }

  const GridViewType& grid_view() const
  {
    return gridView_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(*backend_, entity);
  }

  CommunicatorType& communicator() const
  {
    // no need to prepare the communicator, since we are not pdelab based
    return *communicator_;
  }

private:
  const GridPartType gridPart_;
  const GridViewType gridView_;
  const std::shared_ptr< const BackendType > backend_;
  const std::shared_ptr< const MapperType > mapper_;
  mutable std::shared_ptr< CommunicatorType > communicator_;
}; // class FemBased< ..., 1 >


#else // HAVE_DUNE_FEM


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemBased
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace ContinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH
