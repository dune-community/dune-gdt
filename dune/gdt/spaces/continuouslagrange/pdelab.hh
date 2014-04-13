// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH

#include <memory>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include "../../mapper/pdelab.hh"
#include "../../basefunctionset/pdelab.hh"

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace ContinuousLagrange {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class PdelabBased
{
  static_assert((Dune::AlwaysFalse<GridViewImp>::value), "Untested for this combination of dimensions!");
};


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class PdelabBasedTraits
{
public:
  typedef PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridViewImp GridViewType;
  static const int polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridViewType::ctype DomainFieldType;

public:
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;

private:
  template <class G, bool single_geom, bool is_simplex, bool is_cube>
  struct FeMap
  {
    static_assert(Dune::AlwaysFalse<G>::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template <class G>
  struct FeMap<G, true, true, false>
  {
    typedef PDELab::PkLocalFiniteElementMap<GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  template <class G>
  struct FeMap<G, true, false, true>
  {
    typedef PDELab::QkLocalFiniteElementMap<GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType<GridType>::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                                   == GenericGeometry::SimplexTopology<dimDomain>::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                              == GenericGeometry::CubeTopology<dimDomain>::type::id);
  typedef typename FeMap<GridType, single_geom_, simplicial_, cubic_>::Type FEMapType;

public:
  typedef PDELab::GridFunctionSpace<GridViewType, FEMapType> BackendType;
  typedef Mapper::SimplePdelabWrapper<BackendType> MapperType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef BaseFunctionSet::PdelabWrapper<BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                         dimRangeCols> BaseFunctionSetType;
  static const bool needs_grid_view = true;

private:
  friend class PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class SpaceWrappedFemContinuousLagrangeTraits


template <class GridViewImp, int polynomialOrder, class RangeFieldImp>
class PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public Spaces::ContinuousLagrangeBase<PdelabBasedTraits<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>,
                                            GridViewImp::dimension, RangeFieldImp, 1, 1>
{
  typedef Spaces::ContinuousLagrangeBase<PdelabBasedTraits<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>,
                                         GridViewImp::dimension, RangeFieldImp, 1, 1> BaseType;
  typedef PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;

public:
  typedef PdelabBasedTraits<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridViewType GridViewType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = Traits::dimRange;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::PatternType PatternType;
  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  PdelabBased(const std::shared_ptr<const GridViewType>& gV)
    : gridView_(gV)
    , fe_map_(std::make_shared<FEMapType>(*(gridView_)))
    , backend_(std::make_shared<BackendType>(const_cast<GridViewType&>(*gridView_), *fe_map_))
    , mapper_(std::make_shared<MapperType>(*backend_))
  {
  }

  PdelabBased(const ThisType& other)
    : gridView_(other.gridView_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridView_ = other.gridView_;
      fe_map_   = other.fe_map_;
      backend_  = other.backend_;
      mapper_   = other.mapper_;
    }
    return *this;
  }

  ~PdelabBased()
  {
  }

  const std::shared_ptr<const GridViewType>& grid_view() const
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

private:
  std::shared_ptr<const GridViewType> gridView_;
  std::shared_ptr<const FEMapType> fe_map_;
  std::shared_ptr<const BackendType> backend_;
  std::shared_ptr<const MapperType> mapper_;
}; // class PdelabBased< ..., 1 >


#else // HAVE_DUNE_PDELAB


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class PdelabBased
{
  static_assert((Dune::AlwaysFalse<GridViewImp>::value), "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace ContinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH
