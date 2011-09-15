#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/basefunctionset/continuous/lagrange.hh>
#include <dune/detailed-discretizations/mapper/continuous/lagrange.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace DiscreteFunctionSpace {

namespace Continuous {

template <class FunctionSpaceImp, class GridPartImp, int polOrder>
class Lagrange
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef GridPartImp GridPartType;

  enum
  {
    polynomialOrder = polOrder
  };

  typedef Lagrange<FunctionSpaceType, GridPartType, polynomialOrder> ThisType;

  typedef Dune::DetailedDiscretizations::Mapper::Continuous::Lagrange<FunctionSpaceType, GridPartType, polynomialOrder>
      MapperType;

  typedef Dune::DetailedDiscretizations::BaseFunctionSet::Continuous::Lagrange<ThisType> BaseFunctionSetType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  /**
      @name Convenience typedefs
      @{
   **/
  typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;

  typedef typename IteratorType::Entity EntityType;
  /**
      @}
   **/

  Lagrange(const GridPartType& gridPart)
    : gridPart_(gridPart)
    , mapper_(gridPart_)
    , baseFunctionSet_(*this)
  {
  }

private:
  //! copy constructor
  Lagrange(const ThisType& other)
    : gridPart_(other.gridPart())
    , mapper_(gridPart_)
    , baseFunctionSet_(*this)
  {
  }

public:
  const GridPartType& gridPart() const
  {
    return gridPart_;
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
    if (order() > 0)
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
    return gridPart_.template begin<0>();
  }

  IteratorType end() const
  {
    return gridPart_.template end<0>();
  }
  /**
      @}
   **/

protected:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const GridPartType& gridPart_;
  const MapperType mapper_;
  const BaseFunctionSetType baseFunctionSet_;

}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
