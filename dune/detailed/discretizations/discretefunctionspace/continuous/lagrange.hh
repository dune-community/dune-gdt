#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrangespace.hh>

#include <dune/stuff/la/container/pattern.hh>

#include <dune/detailed/discretizations/basefunctionset/continuous/lagrange.hh>
#include <dune/detailed/discretizations/mapper/continuous/lagrange.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace DiscreteFunctionSpace {

namespace Continuous {

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
class Lagrange
{
public:

  typedef FunctionSpaceImp FunctionSpaceType;

  typedef GridPartImp GridPartType;

  typedef typename GridPartType::GridViewType GridViewType;

  static const int polynomialOrder = polOrder;

  typedef Lagrange< FunctionSpaceType, GridPartType, polynomialOrder > ThisType;

  typedef Dune::Detailed::Discretizations
    ::Mapper::Continuous::Lagrange< FunctionSpaceType, GridPartType, polynomialOrder > MapperType;

  typedef Dune::Detailed::Discretizations
    ::BaseFunctionSet::Continuous::Lagrange< ThisType > BaseFunctionSetType;

//  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

//  typedef typename FunctionSpaceType::DomainType DomainType;

//  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

//  typedef typename FunctionSpaceType::RangeType RangeType;

//  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

//  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  typedef Dune::Stuff::LA::Container::SparsityPatternDefault PatternType;

  Lagrange(const GridPartType& gridPart)
    : gridPart_(gridPart)
    , gridView_(gridPart_.gridView())
    , mapper_(gridPart_)
    , baseFunctionSet_(*this)
  {}

  int order() const
  {
    return polynomialOrder;
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const GridViewType& gridView() const
  {
    return gridView_;
  }

  const MapperType& map() const
  {
    return mapper_;
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return baseFunctionSet_;
  }

  template< class LocalGridPartType, class OtherDiscreteFunctionSpaceType >
  std::shared_ptr< PatternType > computeLocalPattern(const LocalGridPartType& localGridPart,
                                                      const OtherDiscreteFunctionSpaceType& other) const
  {
    typedef typename PatternType::size_type size_type;
    std::shared_ptr< PatternType > ret(new PatternType(mapper_.size()));
    PatternType& pattern = *ret;
    // walk the grid part
    for (typename LocalGridPartType::template Codim< 0 >::IteratorType entityIt = localGridPart.template begin< 0 >();
         entityIt != localGridPart.template end< 0 >();
         ++entityIt) {
      const typename LocalGridPartType::template Codim< 0 >::EntityType& entity = *entityIt;
      for (unsigned int i = 0; i < baseFunctionSet().local(entity).size(); ++i) {
        const size_type globalI = map().toGlobal(entity, i);
        typename PatternType::SetType& columns = pattern.set(globalI);
        for (unsigned int j = 0; j < other.baseFunctionSet().local(entity).size(); ++j) {
          const size_type globalJ = other.map().toGlobal(entity, j);
          columns.insert(globalJ);
        }
      }
    } // walk the grid part
    return ret;
  } // computeLocalPattern()

  template< class LocalGridPartType >
  std::shared_ptr< PatternType > computeLocalPattern(const LocalGridPartType& localGridPart) const
  {
    return computeLocalPattern(localGridPart, *this);
  }

  template< class CouplingGridPartType, class OutsideDiscreteFunctionSpaceType >
  std::shared_ptr< PatternType > computeCouplingPattern(const CouplingGridPartType& couplingGridPart,
                                                         const OutsideDiscreteFunctionSpaceType& outerSpace) const
  {
    typedef typename PatternType::size_type size_type;
    std::shared_ptr< PatternType > ret(new PatternType(mapper_.size()));
    PatternType& pattern = *ret;
    // walk the coupling grid part
    for (typename CouplingGridPartType::template Codim< 0 >::IteratorType entityIt = couplingGridPart.template begin< 0 >();
         entityIt != couplingGridPart.template end< 0 >();
         ++entityIt) {
      // get the inside entity and basefunctionset
      const typename CouplingGridPartType::template Codim< 0 >::EntityType& insideEntity = *entityIt;
      const typename BaseFunctionSetType::LocalBaseFunctionSetType
          ansatzBaseFunctionSet = baseFunctionSet().local(insideEntity);
      // walk the neighbors
      for (typename CouplingGridPartType::IntersectionIteratorType intersectionIt = couplingGridPart.ibegin(insideEntity);
           intersectionIt != couplingGridPart.iend(insideEntity);
           ++intersectionIt) {
        // get the outside neighboring entity and basefunctionset (of the other space)
        const typename CouplingGridPartType::IntersectionIteratorType::Intersection& intersection = *intersectionIt;
        assert(intersection.neighbor() && !intersection.boundary());
        const typename CouplingGridPartType::IntersectionIteratorType::Intersection::EntityPointer outsideNeighborPtr = intersection.outside();
        const typename CouplingGridPartType::template Codim< 0 >::EntityType& outsideNeighbor = *outsideNeighborPtr;
        const typename BaseFunctionSetType::LocalBaseFunctionSetType
            testBaseFunctionSet = outerSpace.baseFunctionSet().local(outsideNeighbor);
        // compute pattern
        for (unsigned int i = 0; i < ansatzBaseFunctionSet.size(); ++i) {
          const size_type globalI = map().toGlobal(insideEntity, i);
          typename PatternType::SetType& columns = pattern.set(globalI);
          for (unsigned int j = 0; j < testBaseFunctionSet.size(); ++j ) {
            const size_type globalJ = outerSpace.map().toGlobal(outsideNeighbor, j);
            columns.insert(globalJ);
          }
        } // compute pattern
      } // walk the neighbors
    } // walk the coupling grid part
    return ret;
  } // computeCouplingPattern()

  template< class OtherDiscreteFunctionSpaceType >
  std::shared_ptr< PatternType > computePattern(const OtherDiscreteFunctionSpaceType& other) const
  {
    return computeLocalPattern(gridPart_, other);
  }

  std::shared_ptr< PatternType > computePattern() const
  {
    return computePattern(*this);
  }

protected:
  Lagrange(const ThisType&);
  ThisType& operator=(const ThisType&);

  const GridPartType& gridPart_;
  const GridViewType gridView_;
  const MapperType mapper_;
  const BaseFunctionSetType baseFunctionSet_;
}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
