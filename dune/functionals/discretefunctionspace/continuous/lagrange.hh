#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

//dune-functionals includes
#include <dune/functionals/basefunctionset/local/lagrange.hh>
#include <dune/functionals/mapper/lagrange.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctionSpace
{

namespace Continuous
{

template< class HostSpaceImp >
class LagrangeFemAdapter;

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
class Lagrange
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef GridPartImp
    GridPartType;

  enum{ polynomialOrder = polOrder };

  typedef Lagrange< FunctionSpaceType, GridPartType, polynomialOrder >
    ThisType;

  typedef Dune::Functionals::Mapper::Lagrange< FunctionSpaceType, GridPartType, polynomialOrder >
    MapperType;

  typedef Dune::Functionals::BaseFunctionSet::Local::Lagrange< ThisType >
    LocalBaseFunctionSetType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

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
      @name Convenience
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
      mapper_( gridPart_ )
  {
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const MapperType& map() const
  {
    return mapper_;
  }

  int order() const
  {
    return polynomialOrder;
  }

  template< class EntityType >
  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  {
    return LocalBaseFunctionSetType( *this, entity );
  }

  bool continuous() const
  {
    return true;
  }

protected:

  template< class >
  friend class Dune::Functionals::DiscreteFunctionSpace::Continuous::LagrangeFemAdapter;

  const GridPartType& gridPart_;
  const MapperType mapper_;

}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
