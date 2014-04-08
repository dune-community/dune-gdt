// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH
#define DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH

#include <type_traits>

#include <dune/common/static_assert.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_PDELAB
# include<dune/pdelab/finiteelementmap/raviartthomasfem.hh>
# include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/gdt/basefunctionset/pdelab.hh>
#include <dune/gdt/mapper/pdelab.hh>

#include "../../../space/interface.hh"

namespace Dune {
namespace GDT {
namespace RaviartThomasSpace {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert(AlwaysFalse< GridViewImp >::value, "Untested for these dimensions or polynomial order!");
}; // class PdelabBased


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols >
class PdelabBasedTraits
{
public:
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridViewImp GridViewType;
  static const int    polOrder = polynomialOrder;
  static_assert(polOrder == 0, "Untested!");
  static_assert(rangeDim == GridViewType::dimension, "Untested!");
  static_assert(rangeDimCols == 1, "Untested!");
private:
  typedef typename GridViewType::ctype  DomainFieldType;
public:
  static const unsigned int             dimDomain = GridViewType::dimension;
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
  static const unsigned int             dimRangeCols = rangeDimCols;
private:
  template< class G, bool single_geom, bool is_simplex, bool is_cube >
  struct FeMap
  {
    static_assert(AlwaysFalse< G >::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template< class G >
  struct FeMap< G, true, true, false >
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap< GridViewType, DomainFieldType,
                                                        RangeFieldType, polOrder, Dune::GeometryType::simplex > Type;
  };
  template< class G >
  struct FeMap< G, true, false, true >
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap< GridViewType, DomainFieldType,
                                                        RangeFieldType, polOrder, Dune::GeometryType::cube > Type;
  };
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                                   == GenericGeometry::SimplexTopology< dimDomain >::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                              == GenericGeometry::CubeTopology< dimDomain >::type::id);
  typedef typename FeMap< GridType, single_geom_, simplicial_, cubic_ >::Type FEMapType;
public:
  typedef PDELab::GridFunctionSpace< GridViewType, FEMapType > BackendType;
  typedef Mapper::SimplePdelabWrapper< BackendType > MapperType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef BaseFunctionSet::PiolaTransformedPdelabWrapper
      < BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >
    BaseFunctionSetType;
  static const bool needs_grid_view = true;
private:
  friend class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >;
}; // class PdelabBasedTraits


template< class GridViewImp, class RangeFieldImp, int rangeDim >
class PdelabBased< GridViewImp, 0, RangeFieldImp, rangeDim, 1 >
  : public SpaceInterface< PdelabBasedTraits< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > >
{
  typedef PdelabBased< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > ThisType;
public:
  typedef PdelabBasedTraits< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > Traits;

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

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  PdelabBased(const std::shared_ptr< const GridViewType >& gV)
    : gridView_(gV)
    , fe_map_(std::make_shared< FEMapType >(*(gridView_)))
    , backend_(std::make_shared< BackendType >(const_cast< GridViewType& >(*gridView_), *fe_map_))
    , mapper_(std::make_shared< MapperType >(*backend_))
  {}

  PdelabBased(const ThisType& other)
    : gridView_(other.gridView_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridView_ = other.gridView_;
      fe_map_ = other.fe_map_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
    }
    return *this;
  }

  ~PdelabBased() {}

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

private:
  std::shared_ptr< const GridViewType > gridView_;
  std::shared_ptr< const FEMapType > fe_map_;
  std::shared_ptr< const BackendType > backend_;
  std::shared_ptr< const MapperType > mapper_;
}; // class PdelabBased< ..., 0, ..., 1 >


#else // HAVE_DUNE_PDELAB


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert(AlwaysFalse< GridViewImp >::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace RaviartThomasSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH
