// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CG_FEM_HH
#define DUNE_GDT_SPACES_CG_FEM_HH

#include <memory>

#include <dune/common/deprecated.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange/space.hh>
#endif // HAVE_DUNE_FEM

#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/spaces/parallel.hh>

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem.hh"

#include "interface.hh"
#include "../constraints.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace CG {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FemBased
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "Untested for these dimensions!");
};


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols>
class FemBasedTraits
{
public:
  typedef FemBased<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;

private:
  typedef Dune::Fem::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, rangeDim> FunctionSpaceType;

public:
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, polOrder> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::BlockMapperType, 1> MapperType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef BaseFunctionSet::FemWrapper<typename BackendType::ShapeFunctionSetType, EntityType, DomainFieldType,
                                      dimDomain, RangeFieldType, rangeDim, rangeDimCols> BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::part;
  static const bool needs_grid_view                       = false;
  typedef CommunicationChooser<GridViewType, false> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class SpaceWrappedFemContinuousLagrangeTraits


// untested for the vector-valued case, especially Spaces::CGInterface
template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemBased<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public Spaces::CGInterface<FemBasedTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>,
                                 GridPartImp::dimension, RangeFieldImp, 1, 1>
{
  typedef Spaces::CGInterface<FemBasedTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>, GridPartImp::dimension,
                              RangeFieldImp, 1, 1> BaseType;
  typedef FemBased<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;

public:
  typedef FemBasedTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  static const int polOrder              = Traits::polOrder;
  static const unsigned int dimDomain    = BaseType::dimDomain;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
  typedef typename Traits::CommunicatorType CommunicatorType;

  typedef typename GridPartType::ctype DomainFieldType;
  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;
  using typename BaseType::DomainType;
  using typename BaseType::BoundaryInfoType;

  explicit FemBased(GridPartType gridP)
    : gridPart_(new GridPartType(gridP))
    , gridView_(new GridViewType(gridPart_->gridView()))
    , backend_(new BackendType(*gridPart_))
    , mapper_(new MapperType(backend_->blockMapper()))
    , communicator_(CommunicationChooserType::create(*gridView_))
  {
  }

  FemBased(const ThisType& other) = default;
  explicit FemBased(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GridPartType& grid_part() const
  {
    return *gridPart_;
  }

  const GridViewType& grid_view() const
  {
    return *gridView_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    return BaseType::lagrange_points_order_1(entity);
  }

  std::set<size_t> local_dirichlet_DoFs(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    return BaseType::local_dirichlet_DoFs_order_1(entity, boundaryInfo);
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
  std::shared_ptr<GridPartType> gridPart_;
  const std::shared_ptr<const GridViewType> gridView_;
  const std::shared_ptr<const BackendType> backend_;
  const std::shared_ptr<const MapperType> mapper_;
  mutable std::shared_ptr<CommunicatorType> communicator_;
}; // class FemBased< ..., 1 >


#else // HAVE_DUNE_FEM


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FemBased
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM


} // namespace CG
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_FEM_HH
