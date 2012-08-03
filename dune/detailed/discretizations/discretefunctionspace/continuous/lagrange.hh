#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

// system
#include <map>
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/gridpart/gridpartview.hh>

// dune-detailed-discretizations includes
#include <dune/detailed/discretizations/basefunctionset/continuous/lagrange.hh>
#include <dune/detailed/discretizations/mapper/continuous/lagrange.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace DiscreteFunctionSpace
{

namespace Continuous
{

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
class Lagrange
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef GridPartImp
    GridPartType;

  typedef Dune::GridPartView< GridPartType > GridViewType;

  enum{ polynomialOrder = polOrder };

  typedef Lagrange< FunctionSpaceType, GridPartType, polynomialOrder >
    ThisType;

  typedef Dune::Detailed::Discretizations::Mapper::Continuous::Lagrange< FunctionSpaceType, GridPartType, polynomialOrder >
    MapperType;

  typedef Dune::Detailed::Discretizations::BaseFunctionSet::Continuous::Lagrange< ThisType >
    BaseFunctionSetType;

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType
    HessianRangeType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  typedef std::map< unsigned int, std::set< unsigned int > > PatternType;

  /**
      @name Convenience typedefs
      @{
   **/
  typedef typename GridPartType::template Codim< 0 >::IteratorType
    IteratorType;

  typedef typename IteratorType::Entity
    EntityType;
  /**
      @}
   **/

  Lagrange( const GridPartType& gridPart )
    : gridPart_( gridPart ),
      gridView_(gridPart_),
      mapper_( gridPart_ ),
      baseFunctionSet_( *this )
  {}

private:
  //! copy constructor
  Lagrange( const ThisType& other );

public:
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

  int order() const
  {
    return polynomialOrder;
  }

  bool continuous() const
  {
    if( order() > 0 )
      return false;
    else
      return true;
  }

  /**
      @name Convenience methods
      @{
   **/
  IteratorType begin() const
  {
    return gridPart_.template begin< 0 >();
  }

  IteratorType end() const
  {
    return gridPart_.template end< 0 >();
  }
  /**
      @}
   **/

  template< class OtherDiscreteFunctionSpaceType>
  PatternType computePattern(const OtherDiscreteFunctionSpaceType& other) const
  {
    std::map< unsigned int, std::set< unsigned int > > ret;
    // generate sparsity pattern
    for (IteratorType it = begin(); it != end(); ++it) {
      const EntityType& entity = *it;
      for(unsigned int i = 0; i < baseFunctionSet().local(entity).size(); ++i) {
        const unsigned int globalI = map().toGlobal(entity, i);
        std::map< unsigned int, std::set< unsigned int > >::iterator result = ret.find(globalI);
        if (result == ret.end())
          ret.insert(std::pair< unsigned int, std::set< unsigned int > >(globalI, std::set< unsigned int >()));
        result = ret.find(globalI);
        assert(result != ret.end());
        for(unsigned int j = 0; j < other.baseFunctionSet().local(entity).size(); ++j) {
          const unsigned int globalJ = other.map().toGlobal(entity, j);
          result->second.insert(globalJ);
        }
      }
    } // generate sparsity pattern
    return ret;
  } // PatternType computePattern(const OtherDiscreteFunctionSpaceType& other) const

  PatternType computePattern() const
  {
    return computePattern(*this);
  }

protected:

  //! assignment operator
  ThisType& operator=( const ThisType& );

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
