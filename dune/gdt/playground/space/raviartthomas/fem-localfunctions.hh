// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_RAVIARTTHOMAS_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_SPACE_RAVIARTTHOMAS_FEM_LOCALFUNCTIONS_HH

#include <type_traits>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_FEM_LOCALFUNCTIONS
# include <dune/localfunctions/raviartthomas.hh>

# include <dune/fem_localfunctions/localfunctions/transformations.hh>
# include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
# include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
# include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>
#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS

#include <dune/stuff/common/color.hh>

#include "../../../mapper/fem.hh"
#include "../../../basefunctionset/fem-localfunctions.hh"
#include "../../../space/constraints.hh"
#include "../../../space/interface.hh"

namespace Dune {
namespace GDT {
namespace RaviartThomasSpace {

#if HAVE_DUNE_FEM_LOCALFUNCTIONS


// forward, to be used in the traits and to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapper
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "Untested for these dimensions!");
};


// forward, to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols >
class FemLocalfunctionsWrapperTraits
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "Untested for these dimensions!");
};


template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 >
{
  static_assert(GridPartImp::dimension == 2, "Only implemented for dimDomain 2!");
  static_assert(polynomialOrder == 0, "Wrong polOrder given!");
public:
  typedef GridPartImp                   GridPartType;
  typedef typename GridPartType::GridViewType GridViewType;
  static const int                      polOrder = polynomialOrder;
private:
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
public:
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = 2;
  static const unsigned int             dimRangeCols = 1;
  typedef FemLocalfunctionsWrapper< GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols > derived_type;
private:
  typedef typename GridPartType::GridType GridType;
  static_assert(Capabilities::hasSingleGeometryType< GridType >::v,
                "This space is only implemented for fully simplicial grids!");
  static_assert(Capabilities::hasSingleGeometryType< GridType >::topologyId
                  == GenericGeometry::SimplexTopology< dimDomain >::type::id,
                "This space is only implemented for fully simplicial grids!");
  typedef Dune::RaviartThomasSimplexLocalFiniteElement< dimDomain,
                                                        DomainFieldType,
                                                        RangeFieldType > FiniteElementType;
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<  GridPartType,
                                                        FiniteElementType,
                                                        Dune::FemLocalFunctions::PiolaTransformation,
                                                        Dune::FemLocalFunctions::SimpleStorage,
                                                        polOrder,
                                                        polOrder,
                                                        true >  BaseFunctionSetMapType;
public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace< BaseFunctionSetMapType >  BackendType;
  typedef Mapper::FemDofWrapper< typename BackendType::MapperType >                 MapperType;
  typedef BaseFunctionSet::FemLocalfunctionsWrapper< BaseFunctionSetMapType,
              DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >  BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType                                  EntityType;
  static const bool needs_grid_view = false;
private:
  template< class G, int p, class R, int r, int rC >
  friend class FemLocalfunctionsWrapper;
}; // class FemLocalfunctionsWrapperTraits< ..., 2, 1 >


template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapper< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 >
  : public SpaceInterface< FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 > >
{
  typedef FemLocalfunctionsWrapper< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 > ThisType;
  typedef SpaceInterface< FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 > > BaseType;
public:
  typedef FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 2, 1 > Traits;

  typedef typename Traits::GridPartType   GridPartType;
  typedef typename Traits::GridViewType   GridViewType;
  typedef typename GridPartType::ctype    DomainFieldType;
  static const int                        polOrder = Traits::polOrder;
  static const unsigned int               dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  FemLocalfunctionsWrapper(const std::shared_ptr< const GridPartType >& gridP)
    : gridPart_(gridP)
    , grid_view_(new GridViewType(gridPart_->gridView()))
    , baseFunctionSetMap_(new BaseFunctionSetMapType(*gridPart_))
    , backend_(new BackendType(const_cast< GridPartType& >(*gridPart_), *baseFunctionSetMap_))
    , mapper_(new MapperType(backend_->mapper()))
  {}

  FemLocalfunctionsWrapper(const ThisType& other)
    : gridPart_(other.gridPart_)
    , grid_view_(other.grid_view_)
    , baseFunctionSetMap_(other.baseFunctionSetMap_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_ = other.gridPart_;
      grid_view_ = other.grid_view_;
      baseFunctionSetMap_ = other.baseFunctionSetMap_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
    }
    return *this;
  }

  const std::shared_ptr< const GridPartType >& grid_part() const
  {
    return gridPart_;
  }

  const std::shared_ptr< const GridViewType >& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  bool continuous() const
  {
    return false;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(*baseFunctionSetMap_, entity);
  }

  template< class R >
  void local_constraints(const EntityType& /*entity*/, Constraints::LocalDefault< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for arbitrary constraints!");
  }

  using BaseType::compute_pattern;

  template< class G, class S >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S >& ansatz_space) const
  {
    return BaseType::compute_face_and_volume_pattern(local_grid_view, ansatz_space);
  }

private:
  std::shared_ptr< const GridPartType > gridPart_;
  std::shared_ptr< const GridViewType > grid_view_;
  std::shared_ptr< BaseFunctionSetMapType  >baseFunctionSetMap_;
  std::shared_ptr< const BackendType > backend_;
  std::shared_ptr< const MapperType > mapper_;
}; // class FemLocalfunctionsWrapper< ..., 2, 1 >


#else // HAVE_DUNE_FEM_LOCALFUNCTIONS


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapper
{
  static_assert(Dune::AlwaysFalse< GridPartImp >::value, "You are missing dune-fem-localfunctions!");
};


#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS

} // namespace RaviartThomasSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_RAVIARTTHOMAS_FEM_LOCALFUNCTIONS_HH
