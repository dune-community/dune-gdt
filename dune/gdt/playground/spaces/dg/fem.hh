// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
#define DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH

#include <memory>

#include <dune/common/unused.hh>
#include <dune/common/deprecated.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#endif

#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/spaces/parallel.hh>

#include "../../../mapper/fem.hh"
#include <dune/gdt/spaces/basefunctionset/dune-fem-wrapper.hh>

#include "../../../spaces/dg/interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace DG {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FemBased
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FemBasedTraits
{
public:
  typedef FemBased<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int polOrder    = polynomialOrder;
  static const bool continuous = false;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const size_t dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;

private:
  typedef Dune::Fem::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, rangeDim> FunctionSpaceType;

public:
  typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polOrder> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::BlockMapperType, BackendType::Traits::localBlockSize> MapperType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef BaseFunctionSet::DuneFemWrapper<typename BackendType::ShapeFunctionSetType, EntityType, DomainFieldType,
                                          dimDomain, RangeFieldType, rangeDim, rangeDimCols> BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::part;
  static const bool needs_grid_view                       = false;
  typedef CommunicationChooser<GridViewType, false> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class FemBasedTraits


} // namespace internal


// untested for the vector-valued case
template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemBased<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public DGInterface<internal::FemBasedTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>,
                         GridPartImp::dimension, 1, 1>
{
  typedef FemBased<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;
  typedef DGInterface<internal::FemBasedTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>,
                      GridPartImp::dimension, 1, 1> BaseType;

public:
  using typename BaseType::Traits;
  typedef typename Traits::GridPartType GridPartType;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::EntityType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;

public:
  using typename BaseType::CommunicatorType;

  FemBased(GridPartType gridP)
    : gridPart_(new GridPartType(gridP))
    , gridView_(new GridViewType(gridPart_->gridView()))
    , backend_(new BackendType(*gridPart_))
    , mapper_(new MapperType(backend_->blockMapper()))
    , communicator_(CommunicationChooserType::create(*gridView_))
  {
  }

  FemBased(const ThisType& other) = default;

  FemBased(ThisType&& source) = default;

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


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FemBased
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM


} // namespace DG
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
