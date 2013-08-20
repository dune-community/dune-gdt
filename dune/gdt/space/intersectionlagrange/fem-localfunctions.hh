// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_INTERSECTION_LAGRANGE_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_SPACE_INTERSECTION_LAGRANGE_FEM_LOCALFUNCTIONS_HH

#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/fem/space/common/allgeomtypes.hh>

//#include <dune/fem_localfunctions/localfunctions/transformations.hh>
#include <dune/fem_localfunctions/localfunctions/interfacelocalfiniteelement.hh>
//#include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
//#include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
#include <dune/fem_localfunctions/basefunctionsetmap/intersectionbasefunctionsetmap.hh>
//#include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>

#include <dune/stuff/common/color.hh>

//#include "../../mapper/fem.hh"
//#include "../../basefunctionset/fem-localfunctions.hh"
//#include "../constraints.hh"
#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace IntersectionLagrangeSpace {


// forward, to be used in the traits and to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapper;


// forward, to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapperTraits;


template< class Codim0GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapperTraits< Codim0GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
{
public:
  typedef typename FemLocalFunctions::IntersectionGridPart< Codim0GridPartImp > GridPartType;
//  typedef GridPartImp                   GridPartType;
  static const int                      polOrder = polynomialOrder;
  dune_static_assert((polOrder >= 0), "ERROR: wrong polOrder given!");
private:
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
public:
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = 1;
  static const unsigned int             dimRangeCols = 1;
  typedef FemLocalfunctionsWrapper< GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols > derived_type;
private:
  typedef FemLocalFunctions::IntersectionLocalFiniteElement<  LagrangeLocalFiniteElement< Dune::EquidistantPointSet,
                                                              dimDomain,
                                                              DomainFieldType,
                                                              RangeFieldType > > FiniteElementType;
  typedef FemLocalFunctions::IntersectionBaseFunctionSetMap<  Codim0GridPartImp,
                                                              FiniteElementType,
                                                              Dune::FemLocalFunctions::NoTransformation,
                                                              Dune::FemLocalFunctions::SimpleStorage,
                                                              polOrder,
                                                              polOrder >  BaseFunctionSetMapType;
public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace< BaseFunctionSetMapType >  BackendType;
  typedef Mapper::FemDofWrapper< typename BackendType::MapperType >                 MapperType;
  typedef BaseFunctionSet::FemLocalfunctionsWrapper< BaseFunctionSetMapType,
              DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >  BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType                                  EntityType;
private:
  template< class G, int p, class R, int r, int rC >
  friend class FemLocalfunctionsWrapper;
}; // class FemLocalfunctionsWrapperTraits< ..., 1, 1 >


template< class Codim0GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapper< Codim0GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
//  : public SpaceInterface< FemLocalfunctionsWrapperTraits< Codim0GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > >
{
public:
  typedef FemLocalfunctionsWrapperTraits< Codim0GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > Traits;

  typedef typename Traits::GridPartType   GridPartType;
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
private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  FemLocalfunctionsWrapper(const GridPartType& gridP)
    : gridPart_(assertGridPart(gridP))
    , baseFunctionSetMap_(gridPart_)
    , backend_(const_cast< GridPartType& >(gridPart_), baseFunctionSetMap_)
    , mapper_(backend_.mapper())
  {}

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  bool continuous() const
  {
    return false;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(baseFunctionSetMap_, entity);
  }

  template< class R >
  void localConstraints(const EntityType& /*entity*/,
                        Constraints::LocalDefault< R >& /*ret*/) const
  {
    dune_static_assert(Dune::AlwaysFalse< R >::value, "ERROR: not implemented for arbitrary constraints!");
  }

private:
  static const GridPartType& assertGridPart(const GridPartType& gP)
  {
    // check
    typedef typename Dune::Fem::AllGeomTypes< typename GridPartType::BaseType::IndexSetType,
                                              typename GridPartType::BaseType::GridType > AllGeometryTypes;
    const AllGeometryTypes allGeometryTypes(gP.baseGridPart().indexSet());
    const std::vector< Dune::GeometryType >& geometryTypes = allGeometryTypes.geomTypes(1);
    if (!(geometryTypes.size() == 1 && geometryTypes[0].isSimplex()))
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " this space is only implemented for simplicial intersections!");
    return gP;
  } // ... assertGridPart(...)

  const GridPartType& gridPart_;
  BaseFunctionSetMapType baseFunctionSetMap_;
  const BackendType backend_;
  const MapperType mapper_;
}; // class FemLocalfunctionsWrapper< ..., 1, 1 >


} // namespace IntersectionLagrangeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_INTERSECTION_LAGRANGE_FEM_LOCALFUNCTIONS_HH
