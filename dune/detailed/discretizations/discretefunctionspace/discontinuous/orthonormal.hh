#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_DISCONTINUOUS_ORTHONORMAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_DISCONTINUOUS_ORTHONORMAL_HH

// dune-detailed-discretizations includes
#include <dune/detailed/discretizations/basefunctionset/orthonormal.hh>
#include <dune/detailed/discretizations/mapper/discontinuous.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace DiscreteFunctionSpace
{

namespace Discontinuous
{

template< class FunctionSpaceImp, class GridViewImp, int polOrder >
class Orthonormal
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef GridViewImp
    GridViewType;

  enum{ polynomialOrder = polOrder };

  typedef Orthonormal< FunctionSpaceType, GridViewType, polynomialOrder >
    ThisType;

  typedef Dune::Detailed::Discretizations::Mapper::Discontinuous< ThisType >
    MapperType;

  typedef Dune::Detailed::Discretizations::BaseFunctionSet::Orthonormal< ThisType >
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

  /**
      @name Convenience typedefs
      @{
   **/
  typedef typename GridViewType::template Codim< 0 >::Iterator
    GridElementIteratorType;

  typedef typename GridElementIteratorType::Entity
    GridElementType;
  /**
      @}
   **/

  Orthonormal( const GridViewType& gridView )
    : gridView_( gridView ),
      baseFunctionSet_( *this ),
      mapper_( *this )
  {
  }

public:
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
    return false;
  }

  /**
      @name Convenience methods
      @{
   **/
  GridElementIteratorType gridElementBegin() const
  {
    return gridView_.template begin< 0 >();
  }

  GridElementIteratorType gridElementEnd() const
  {
    return gridView_.template end< 0 >();
  }
  /**
      @}
   **/

protected:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  //! copy constructor
  Orthonormal( const ThisType& other );

  const GridViewType& gridView_;
  const BaseFunctionSetType baseFunctionSet_;
  const MapperType mapper_;

}; // end class Orthonormal

} // end namespace Discontinuous

} // end namespace DiscreteFunctionSpace

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_DISCONTINUOUS_ORTHONORMAL_HH
