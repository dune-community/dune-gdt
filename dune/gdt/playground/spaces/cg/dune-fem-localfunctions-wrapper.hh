// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_CG_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_CG_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_FEM_LOCALFUNCTIONS
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/fem_localfunctions/localfunctions/transformations.hh>
#include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
#include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
#include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>
#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS

#include <dune/xt/common/type_traits.hh>

#include <dune/gdt/spaces/mapper/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/basefunctionset/dune-fem-localfunctions-wrapper.hh>
#include <dune/gdt/spaces/cg/interface.hh>

namespace Dune {
namespace GDT {

#if HAVE_DUNE_FEM_LOCALFUNCTIONS


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DuneFemLocalfunctionsCgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "Untested for these dimensions!");
};


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class DuneFemLocalfunctionsCgSpaceWrapperTraits
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int polOrder = polynomialOrder;
  static const bool continuous = true;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const size_t dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  typedef DuneFemLocalfunctionsCgSpaceWrapper<GridPartType, polOrder, RangeFieldType, rangeDim, rangeDimCols>
      derived_type;

private:
  typedef typename GridPartType::GridType GridType;
  static_assert(dimDomain == 1 || Dune::Capabilities::hasSingleGeometryType<GridType>::v,
                "This space is only implemented for fully simplicial grids!");
  static_assert(dimDomain == 1 || (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                                   == GenericGeometry::SimplexTopology<dimDomain>::type::id),
                "This space is only implemented for fully simplicial grids!");
  typedef DuneFemLocalfunctionsCgSpaceWrapperTraits<GridPartType, polOrder, RangeFieldType, rangeDim, rangeDimCols>
      ThisType;

public:
  typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dimDomain, DomainFieldType, RangeFieldType>
      FiniteElementType;

private:
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<GridPartType,
                                                      FiniteElementType,
                                                      Dune::FemLocalFunctions::NoTransformation,
                                                      Dune::FemLocalFunctions::SimpleStorage,
                                                      polOrder,
                                                      polOrder>
      BaseFunctionSetMapType;

public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace<BaseFunctionSetMapType> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::MapperType> MapperType;
  typedef BaseFunctionSet::DuneFemLocalfunctionsWrapper<BaseFunctionSetMapType,
                                                        DomainFieldType,
                                                        dimDomain,
                                                        RangeFieldType,
                                                        rangeDim,
                                                        rangeDimCols>
      BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType EntityType;
  static const XT::Grid::Backends part_view_type = XT::Grid::Backends::part;
  static const bool needs_grid_view = false;
  typedef double CommunicatorType;

private:
  template <class G, int p, class R, size_t r, size_t rC>
  friend class DuneFemLocalfunctionsCgSpaceWrapper;
}; // class DuneFemLocalfunctionsCgSpaceWrapperTraits


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class DuneFemLocalfunctionsCgSpaceWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public CgSpaceInterface<DuneFemLocalfunctionsCgSpaceWrapperTraits<GridPartImp,
                                                                        polynomialOrder,
                                                                        RangeFieldImp,
                                                                        1,
                                                                        1>,
                              GridPartImp::dimension,
                              RangeFieldImp,
                              1,
                              1>
{
  typedef CgSpaceInterface<DuneFemLocalfunctionsCgSpaceWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>,
                           GridPartImp::dimension,
                           RangeFieldImp,
                           1,
                           1>
      BaseType;
  typedef DuneFemLocalfunctionsCgSpaceWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;

public:
  typedef DuneFemLocalfunctionsCgSpaceWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::GridViewType GridViewType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype DomainFieldType;
  static const size_t dimDomain = BaseType::dimension;

private:
  static_assert(GridPartType::dimension == dimDomain, "Dimension of GridPart has to match dimDomain");

public:
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t dimRangeCols = BaseType::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

  typedef Dune::XT::LA::SparsityPatternDefault PatternType;
  using typename BaseType::DomainType;
  using typename BaseType::BoundaryInfoType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  DuneFemLocalfunctionsCgSpaceWrapper(GridPartType gridP)
    : gridPart_(gridP)
    , gridView_(gridPart_.gridView()))
    , baseFunctionSetMap_(gridPart_)
    , backend_(const_cast< GridPartType& >(gridPart_), baseFunctionSetMap_)
    , mapper_(backend_->mapper())
    , tmp_global_indices_(mapper_->maxNumDofs())
    , communicator_(0.0)
  {
  }

  DuneFemLocalfunctionsCgSpaceWrapper(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_ = other.gridPart_;
      gridView_ = other.gridView_;
      baseFunctionSetMap_ = other.baseFunctionSetMap_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
      tmp_global_indices_.resize(mapper_->maxNumDofs());
    }
    return *this;
  }

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
    return backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return mapper_;
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
    return BaseFunctionSetType(baseFunctionSetMap_, entity);
  }

  double& communicator() const
  {
    return communicator_;
  }

private:
  const GridPartType gridPart_;
  const GridViewType gridView_;
  BaseFunctionSetMapType baseFunctionSetMap_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable Dune::DynamicVector<size_t> tmp_global_indices_;
  mutable double communicator_;
}; // class DuneFemLocalfunctionsCgSpaceWrapper< ..., 1, 1 >


#else // HAVE_DUNE_FEM_LOCALFUNCTIONS


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DuneFemLocalfunctionsCgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "You are missing dune-fem-localfunctions!");
};


#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_CG_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH
