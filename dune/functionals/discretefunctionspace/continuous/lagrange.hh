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

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
class Lagrange
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef GridPartImp
    GridPartType;

  typedef Lagrange< FunctionSpaceType, GridPartType, polOrder >
    ThisType;

//  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
//    HostSpaceType;

  typedef Dune::Functionals::Mapper::Lagrange< FunctionSpaceType, GridPartType, polOrder >
    MapperType;

  typedef Dune::Functionals::BaseFunctionSet::Local::Lagrange< ThisType >
    LocalBaseFunctionSetType;

//  typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity
//    EntityType;

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

//  /**
//    \defgroup dune-fem related
//    \{
//    **/
//  typedef typename HostSpaceType::BaseFunctionSetType
//    BaseFunctionSetType;

//  typedef typename HostSpaceType::IteratorType
//    IteratorType;
//  /**
//    \}
//    **/

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  Lagrange( GridPartType& gridPart )
    : gridPart_( gridPart )/*,
      hostSpace_( gridPart_ )*/,
      mapper_( gridPart_ )/*,
      numMaxLocalDoFs_( -1 )*/
  {
//    // in the simple case, there should be the same number of dofs on each entity
//    const IteratorType entityIterator = begin();
//    const EntityType& entity = *entityIterator;
//    numMaxLocalDoFs_ = hostSpace_.baseFunctionSet( entity ).numBaseFunctions();
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

//  const HostSpaceType& hostSpace() const
//  {
//    return hostSpace_;
//  }

  const MapperType& map() const
  {
    return mapper_;
  }

  template< class EntityType >
  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  {
    return LocalBaseFunctionSetType( *this, entity );
  }

//  const unsigned int size() const
//  {
//    return hostSpace_.size();
//  }

//  const int numMaxLocalDoFs() const
//  {
//    return numMaxLocalDoFs_;
//  }

//  const int order() const
//  {
//    return hostSpace_.order();
//  }

//  /**
//    \defgroup dune-fem related
//    \{
//    **/
//  IteratorType begin() const
//  {
//    return hostSpace_.begin();
//  }

//  const IteratorType end() const
//  {
//    return hostSpace_.end();
//  }

//  const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
//  {
//    return hostSpace_.baseFunctionSet( entity );
//  }

//  int mapToGlobal( const EntityType& entity, const int localDof) const
//  {
//    return hostSpace_.mapToGlobal( entity, localDof);
//  }
//  /**
//    \}
//    **/

private:

  const GridPartType& gridPart_;
//  const HostSpaceType hostSpace_;
  const MapperType mapper_;
//  unsigned int numMaxLocalDoFs_;

}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
