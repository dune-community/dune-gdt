#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrangespace.hh>

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem.hh"

#include "../interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace ContinuousLagrangeSpace {


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class FemWrapper;


// forward, to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class FemWrapperTraits;


/**
 *  \brief Traits class for ContinuousLagrangeSpace for dimRange 1x1.
 */
template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
{
public:
  typedef FemWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> derived_type;
  typedef GridPartImp GridPartType;
  static const int polOrder = polynomialOrder;
  dune_static_assert((polOrder >= 1), "ERROR: wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = 1;
  static const unsigned int dimRangeCols = 1;

private:
  typedef Dune::Fem::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, dimRangeRows> FunctionSpaceType;

public:
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, polOrder> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::MapperType> MapperType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef BaseFunctionSet::FemWrapper<typename BackendType::BaseFunctionSetType, EntityType> BaseFunctionSetType;
}; // class SpaceWrappedFemContinuousLagrangeTraits< ..., 1, 1 >


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public SpaceInterface<FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>>
{
public:
  typedef FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const int polOrder           = Traits::polOrder;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = Traits::dimRangeRows;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

public:
  FemWrapper(const GridPartType& gridP)
    : gridPart_(gridP)
    , backend_(const_cast<GridPartType&>(gridPart_))
    , mapper_(backend_.mapper())
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
    return BaseFunctionSetType(backend_, entity);
  }

private:
  const GridPartType& gridPart_;
  const BackendType backend_;
  const MapperType mapper_;
}; // class FemWrapper< ..., 1, 1 >


} // namespace ContinuousLagrangeSpace
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH
