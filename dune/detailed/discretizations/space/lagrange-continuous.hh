#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/fem_localfunctions/localfunctions/transformations.hh>
#include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
#include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
#include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>

#include "../mapper/fem.hh"
#include "../basefunctionset/fem-localfunctions.hh"
#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, class RangeFieldImp, unsigned int rangeDim, int polynomialOrder>
class ContinuousLagrangeSpace;


// forward, to allow for specialization
template <class GridPartImp, class RangeFieldImp, unsigned int rangeDim, int polynomialOrder>
class ContinuousLagrangeSpaceTraits;


/**
 *  \brief Traits class for ContinuousLagrangeSpace for dimRange 1.
 */
template <class GridPartImp, class RangeFieldImp, int polynomialOrder>
class ContinuousLagrangeSpaceTraits<GridPartImp, RangeFieldImp, 1, polynomialOrder>
{
public:
  typedef GridPartImp GridPartType;

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = 1;
  static const int polOrder          = polynomialOrder;
  typedef ContinuousLagrangeSpace<GridPartType, RangeFieldType, dimRange, polOrder> derived_type;

private:
  typedef ContinuousLagrangeSpaceTraits<GridPartType, RangeFieldType, dimRange, polOrder> ThisType;
  typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dimDomain, DomainFieldType, RangeFieldType>
      FiniteElementType;
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<GridPartType, FiniteElementType,
                                                      Dune::FemLocalFunctions::NoTransformation,
                                                      Dune::FemLocalFunctions::SimpleStorage, polOrder,
                                                      polOrder> BaseFunctionSetMapType;

public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace<BaseFunctionSetMapType> BackendType;
  typedef Mapper::FemWrapper<ThisType> MapperType;
  typedef BaseFunctionSet::FemLocalfunctionsWrapper<ThisType> BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType EntityType;

private:
  friend class ContinuousLagrangeSpace<GridPartType, RangeFieldType, dimRange, polOrder>;
};


template <class GridPartImp, class RangeFieldImp, unsigned int rangeDim, int polynomialOrder>
class ContinuousLagrangeSpace
    : public SpaceInterface<ContinuousLagrangeSpaceTraits<GridPartImp, RangeFieldImp, rangeDim, polynomialOrder>>
{
public:
  typedef ContinuousLagrangeSpace<GridPartImp, RangeFieldImp, rangeDim, polynomialOrder> ThisType;
  typedef ContinuousLagrangeSpaceTraits<GridPartImp, RangeFieldImp, rangeDim, polynomialOrder> Traits;
  typedef SpaceInterface<Traits> InterfaceType;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = Traits::dimRange;
  static const int polOrder          = Traits::polOrder;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  ContinuousLagrangeSpace(const GridPartType& gridPart)
    : gridPart_(gridPart)
    , baseFunctionSetMap_(gridPart_)
    , backend_(const_cast<GridPartType&>(gridPart_), baseFunctionSetMap_)
    , mapper_(*this)
  {
  }

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
    return true;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(*this, entity);
  }

private:
  const GridPartType& gridPart_;
  BaseFunctionSetMapType baseFunctionSetMap_;
  const BackendType backend_;
  const MapperType mapper_;
}; // class ContinuousLagrangeSpace


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH
