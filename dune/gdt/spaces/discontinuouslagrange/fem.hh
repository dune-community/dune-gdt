// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
#define DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH

#include "config.h"

#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
//# include <dune/fem/space/lagrange/space.hh>
#endif // HAVE_DUNE_FEM

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem.hh"

#include "../constraints.hh"
#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace DiscontinuousLagrange {

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
  typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > BackendType;
  typedef Mapper::FemDofWrapper< typename BackendType::BlockMapperType > MapperType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef BaseFunctionSet::FemWrapper
      < typename BackendType::ShapeFunctionSetType, EntityType, DomainFieldType, dimDomain,
        RangeFieldType, dimRange, dimRangeCols > BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::part;
  static const bool needs_grid_view = false;
  typedef double CommunicatorType;
}; // class SpaceWrappedFemDiscontinuousLagrangeTraits


// untested for the vector-valued case
template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemBased< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
    : public SpaceInterface< FemBasedTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > >
{
  typedef FemBased< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >               ThisType;
  typedef SpaceInterface< FemBasedTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > > BaseType;
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

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  FemBased(const std::shared_ptr< const GridPartType >& gridP)
    : gridPart_(gridP)
    , gridView_(std::make_shared< GridViewType >(gridPart_->gridView()))
    , backend_(std::make_shared< BackendType >(const_cast< GridPartType& >(*(gridPart_))))
    , mapper_(std::make_shared< MapperType >(backend_->blockMapper()))
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
    , communicator_(0.0)
  {}

  FemBased(const ThisType& other)
    : gridPart_(other.gridPart_)
    , gridView_(other.gridView_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
    , communicator_(0.0)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_ = other.gridPart_;
      gridView_ = other.gridView_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
      tmpMappedRows_.resize(mapper_->maxNumDofs());
      tmpMappedCols_.resize(mapper_->maxNumDofs());
    }
    return *this;
  }

  ~FemBased() {}

  using BaseType::compute_pattern;

  template< class G, class S >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S >& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(local_grid_view, ansatz_space);
  }

  const std::shared_ptr< const GridPartType >& grid_part() const
  {
    return gridPart_;
  }

  const std::shared_ptr< const GridViewType >& grid_view() const
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

  double& communicator() const
  {
    return communicator_;
  }

private:
  std::shared_ptr< const GridPartType > gridPart_;
  std::shared_ptr< const GridViewType > gridView_;
  std::shared_ptr< const BackendType > backend_;
  std::shared_ptr< const MapperType > mapper_;
  mutable Dune::DynamicVector< size_t > tmpMappedRows_;
  mutable Dune::DynamicVector< size_t > tmpMappedCols_;
  mutable double communicator_;
}; // class FemBased< ..., 1 >


#else // HAVE_DUNE_FEM


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemBased
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace DiscontinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
