// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH

#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>

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

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem-localfunctions.hh"
#include "../continuouslagrange.hh"

namespace Dune {
namespace GDT {
namespace ContinuousLagrangeSpace {

#if HAVE_DUNE_FEM_LOCALFUNCTIONS


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FemLocalfunctionsWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "Untested for these dimensions!");
};


/**
 *  \brief Traits class for ContinuousLagrangeSpace::FemLocalfunctionsWrapper.
 */
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FemLocalfunctionsWrapperTraits
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;

public:
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef FemLocalfunctionsWrapper<GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols> derived_type;

private:
  typedef typename GridPartType::GridType GridType;
  static_assert(dimDomain == 1 || Dune::Capabilities::hasSingleGeometryType<GridType>::v,
                "This space is only implemented for fully simplicial grids!");
  static_assert(dimDomain == 1 || (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                                   == GenericGeometry::SimplexTopology<dimDomain>::type::id),
                "This space is only implemented for fully simplicial grids!");
  typedef FemLocalfunctionsWrapperTraits<GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols> ThisType;

public:
  typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dimDomain, DomainFieldType, RangeFieldType>
      FiniteElementType;

private:
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<GridPartType, FiniteElementType,
                                                      Dune::FemLocalFunctions::NoTransformation,
                                                      Dune::FemLocalFunctions::SimpleStorage, polOrder,
                                                      polOrder> BaseFunctionSetMapType;

public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace<BaseFunctionSetMapType> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::MapperType> MapperType;
  typedef BaseFunctionSet::FemLocalfunctionsWrapper<BaseFunctionSetMapType, DomainFieldType, dimDomain, RangeFieldType,
                                                    dimRange, dimRangeCols> BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType EntityType;
  static const bool needs_grid_view = false;

private:
  template <class G, int p, class R, int r, int rC>
  friend class FemLocalfunctionsWrapper;
}; // class FemLocalfunctionsWrapperTraits


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemLocalfunctionsWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public ContinuousLagrangeSpaceBase<FemLocalfunctionsWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1,
                                                                        1>,
                                         GridPartImp::dimension, RangeFieldImp, 1, 1>
{
  typedef ContinuousLagrangeSpaceBase<FemLocalfunctionsWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>,
                                      GridPartImp::dimension, RangeFieldImp, 1, 1> BaseType;
  typedef FemLocalfunctionsWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;

public:
  typedef FemLocalfunctionsWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::GridViewType GridViewType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = Traits::dimRange;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  FemLocalfunctionsWrapper(std::shared_ptr<const GridPartType> gridP)
    : gridPart_(gridP)
    , gridView_(std::make_shared<GridViewType>(gridPart_->gridView()))
    , baseFunctionSetMap_(new BaseFunctionSetMapType(*gridPart_))
    , backend_(new BackendType(const_cast<GridPartType&>(*gridPart_), *baseFunctionSetMap_))
    , mapper_(new MapperType(backend_->mapper()))
    , tmp_global_indices_(mapper_->maxNumDofs())
  {
  }

  FemLocalfunctionsWrapper(const ThisType& other)
    : gridPart_(other.gridPart_)
    , gridView_(other.gridView_)
    , baseFunctionSetMap_(other.baseFunctionSetMap_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , tmp_global_indices_(mapper_->maxNumDofs())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_           = other.gridPart_;
      gridView_           = other.gridView_;
      baseFunctionSetMap_ = other.baseFunctionSetMap_;
      backend_            = other.backend_;
      mapper_ = other.mapper_;
      tmp_global_indices_.resize(mapper_->maxNumDofs());
    }
    return *this;
  }

  const std::shared_ptr<const GridPartType>& grid_part() const
  {
    return gridPart_;
  }

  const std::shared_ptr<const GridViewType>& grid_view() const
  {
    return gridView_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(*baseFunctionSetMap_, entity);
  }

private:
  std::shared_ptr<const GridPartType> gridPart_;
  std::shared_ptr<const GridViewType> gridView_;
  std::shared_ptr<BaseFunctionSetMapType> baseFunctionSetMap_;
  std::shared_ptr<const BackendType> backend_;
  std::shared_ptr<const MapperType> mapper_;
  mutable Dune::DynamicVector<size_t> tmp_global_indices_;
}; // class FemLocalfunctionsWrapper< ..., 1, 1 >


#else // HAVE_DUNE_FEM_LOCALFUNCTIONS


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FemLocalfunctionsWrapper
{
  static_assert(Dune::AlwaysFalse<GridPartImp>::value, "You are missing dune-fem-localfunctions!");
};


#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS

} // namespace ContinuousLagrangeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH
