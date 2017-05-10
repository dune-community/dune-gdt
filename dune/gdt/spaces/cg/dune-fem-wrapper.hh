// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_HH
#define DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_HH

#include <memory>

#include <dune/common/deprecated.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange/space.hh>
#endif

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"
#include "../mapper/dune-fem-wrapper.hh"
#include "../basefunctionset/dune-fem-wrapper.hh"
#include "../constraints.hh"

namespace Dune {
namespace GDT {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DuneFemCgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "Untested for these dimensions!");
};


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class DuneFemCgSpaceWrapperTraits
{
  static_assert(XT::Grid::is_layer<GridPartImp>::value && !XT::Grid::is_view<GridPartImp>::value, "");

public:
  typedef DuneFemCgSpaceWrapper<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridPartImp GridLayerType;
  static const int polOrder = polynomialOrder;
  static const bool continuous = true;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = GridLayerType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;

private:
  typedef Dune::Fem::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, rangeDim> FunctionSpaceType;

public:
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridLayerType, polOrder> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::BlockMapperType, BackendType::Traits::localBlockSize> MapperType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef BaseFunctionSet::DuneFemWrapper<typename BackendType::BasisFunctionSetType,
                                          EntityType,
                                          DomainFieldType,
                                          dimDomain,
                                          RangeFieldType,
                                          rangeDim,
                                          rangeDimCols>
      BaseFunctionSetType;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::part;
  static const bool needs_grid_view = false;
  typedef CommunicationChooser<GridLayerType, false> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class DuneFemCgSpaceWrapperTraits


// untested for the vector-valued case, especially CgSpaceInterface
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t r>
class DuneFemCgSpaceWrapper<GridPartImp, polynomialOrder, RangeFieldImp, r, 1>
    : public CgSpaceInterface<DuneFemCgSpaceWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, r, 1>,
                              GridPartImp::dimension,
                              r,
                              1>
{
  typedef CgSpaceInterface<DuneFemCgSpaceWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, r, 1>,
                           GridPartImp::dimension,
                           r,
                           1>
      BaseType;
  typedef DuneFemCgSpaceWrapper<GridPartImp, polynomialOrder, RangeFieldImp, r, 1> ThisType;

public:
  typedef DuneFemCgSpaceWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, r, 1> Traits;

  static const int polOrder = Traits::polOrder;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t dimRangeCols = BaseType::dimRangeCols;

  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;

public:
  typedef typename Traits::CommunicatorType CommunicatorType;

  using typename BaseType::DomainType;

  explicit DuneFemCgSpaceWrapper(GridLayerType grd_prt)
    : grid_part_(new GridLayerType(grd_prt))
    , backend_(new BackendType(*grid_part_))
    , mapper_(new MapperType(backend_->blockMapper()))
    , communicator_(CommunicationChooserType::create(*grid_part_))
  {
  }

  DuneFemCgSpaceWrapper(const ThisType& other) = default;
  explicit DuneFemCgSpaceWrapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& DUNE_DEPRECATED_MSG("Use grid_layer() instead (03.04.2017)!") grid_part() const
  {
    return *grid_part_;
  }

  const GridLayerType& grid_layer() const
  {
    return *grid_part_;
  }

  GridLayerType& grid_part()
  {
    return *grid_part_;
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

  std::set<size_t> local_dirichlet_DoFs(
      const EntityType& entity,
      const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundaryInfo) const
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
  std::shared_ptr<GridLayerType> grid_part_;
  const std::shared_ptr<const BackendType> backend_;
  const std::shared_ptr<const MapperType> mapper_;
  mutable std::shared_ptr<CommunicatorType> communicator_;
}; // class DuneFemCgSpaceWrapper< ..., 1 >


#else // HAVE_DUNE_FEM


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DuneFemCgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM


} // namespace GDT
} // namespace Dune


#include "dune-fem-wrapper.lib.hh"


#endif // DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_HH
