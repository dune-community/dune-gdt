#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

//dune-functionals includes
#include <dune/functionals/common/localbasefunctionset.hh>

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

  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
    HostSpaceType;

  typedef Dune::Functionals::Common::LocalBaseFunctionSet< ThisType >
    LocalBaseFunctionSetType;

  typedef typename HostSpaceType::EntityType
    EntityType;

  typedef typename HostSpaceType::DomainType
    DomainType;

  typedef typename HostSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename HostSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename HostSpaceType::RangeType
    RangeType;

  typedef typename HostSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef typename HostSpaceType::HessianRangeType
    HessianRangeType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename HostSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename HostSpaceType::IteratorType
    IteratorType;
  /**
    \}
    **/

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

  Lagrange( GridPartType& gridPart )
    : gridPart_( gridPart ),
      hostSpace_( gridPart ),
      numMaxLocalDoFs_( -1 )
  {
    // in the simple case, there should be the same number of dofs on each entity
    const IteratorType entityIterator = begin();
    const EntityType& entity = *entityIterator;
    numMaxLocalDoFs_ = hostSpace_.baseFunctionSet( entity ).numBaseFunctions();
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const HostSpaceType& hostSpace() const
  {
    return hostSpace_;
  }

  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  {
    return LocalBaseFunctionSetType( *this, entity );
  }

  const unsigned int size() const
  {
    return hostSpace_.size();
  }

  const int numMaxLocalDoFs() const
  {
    return numMaxLocalDoFs_;
  }

  const int order() const
  {
    return hostSpace_.order();
  }

  /**
    \defgroup dune-fem related
    \{
    **/
  IteratorType begin() const
  {
    return hostSpace_.begin();
  }

  const IteratorType end() const
  {
    return hostSpace_.end();
  }

  const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
  {
    return hostSpace_.baseFunctionSet( entity );
  }

  const int mapToGlobal( const EntityType& entity, const int localDof) const
  {
    return hostSpace_.mapToGlobal( entity, localDof);
  }
  /**
    \}
    **/

private:

  const GridPartType& gridPart_;
  const HostSpaceType hostSpace_;
  unsigned int numMaxLocalDoFs_;

}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
