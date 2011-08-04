#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH

// dune-fem-tools includes, should be removed in the end!
//#include <dune/fem-tools/grid/entity.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunction
{

template< class DiscreteFunctionImp >
class LocalConst
{
public:

  typedef DiscreteFunctionImp
    DiscreteFunctionType;

  typedef LocalConst< DiscreteFunctionType >
    ThisType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType::template Codim< 0 >::IteratorType::Entity
    EntityType;

  typedef typename DiscreteFunctionType::DomainType
    DomainType;

  typedef typename DiscreteFunctionType::RangeFieldType
    RangeFieldType;

  typedef typename DiscreteFunctionType::RangeType
    RangeType;

  typedef typename DiscreteFunctionType::JacobianRangeType
    JacobianRangeType;

private:
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType
    LocalBaseFunctionSetType;

public:
  LocalConst( const DiscreteFunctionType& discreteFunction, const EntityType& entity )
    : discreteFunction_( discreteFunction ),
      entity_( entity ),
      localBaseFunctionSet_( discreteFunction_.space().baseFunctionSet().local( entity_ ) ),
      size_( localBaseFunctionSet_.size() ),
      order_( localBaseFunctionSet_.order() )
  {
//std::cout << "DiscreteFunction::LocalConst::LocalConst()" << std::endl;
//Dune::FemTools::Grid::Entity::print( entity_, std::cout, "  " );
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const RangeFieldType& operator[]( const unsigned int localDofNumber ) const
  {
    const unsigned int globalDofNumber = discreteFunction_.space().map().toGlobal( entity_, localDofNumber );
    return discreteFunction_[globalDofNumber];
  }

  int order() const
  {
    return order_;
  }

  unsigned int size() const
  {
    return size_;
  }

  void evaluate( const DomainType& x, RangeType& ret) const
  {
//std::cout << "DiscreteFunction::LocalConst::evaluate()" << std::endl;
//Dune::FemTools::Grid::Entity::print( entity_, std::cout, "  " );
    std::vector< RangeType > baseFunctionValues( size(), RangeType( 0.0 ) );
    localBaseFunctionSet_.evaluate( x, baseFunctionValues );
    ret = 0.0;
    for( unsigned int i = 0; i < size(); ++i )
    {
      baseFunctionValues[i] *= operator[](i);
      ret += baseFunctionValues[i];
    }
  }

  void jacobian( const DomainType& x, JacobianRangeType& ret ) const
  {
std::cout << "DiscreteFunction::LocalConst::jacobian()" << std::endl;
Dune::FemTools::Grid::Entity::print( entity_, std::cout, "  " );
    std::vector< JacobianRangeType > baseFunctionJacobianValues( size(), JacobianRangeType( 0.0 ) );
    localBaseFunctionSet_.jacobian( x, baseFunctionJacobianValues );
    ret = 0.0;
    for( unsigned int i = 0; i < size(); ++i )
    {
      baseFunctionJacobianValues[i] *= operator[](i);
      ret += baseFunctionJacobianValues[i];
    }
  }

private:

  const DiscreteFunctionType& discreteFunction_;
  const EntityType& entity_;
  const LocalBaseFunctionSetType& localBaseFunctionSet_;
  const unsigned int size_;
  const int order_;

}; // end class LocalConst

template< class DiscreteFunctionImp >
class Local
  : public LocalConst< DiscreteFunctionImp >
{
public:

  typedef DiscreteFunctionImp
    DiscreteFunctionType;

  typedef Local< DiscreteFunctionType >
    ThisType;

  typedef LocalConst< DiscreteFunctionType >
    BaseType;

  typedef typename BaseType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename BaseType::EntityType
    EntityType;

  typedef typename BaseType::DomainType
    DomainType;

  typedef typename BaseType::RangeFieldType
    RangeFieldType;

  typedef typename BaseType::RangeType
    RangeType;

  typedef typename BaseType::JacobianRangeType
    JacobianRangeType;

  Local( DiscreteFunctionType& discreteFunction, const EntityType& entity )
    : BaseType( const_cast< const DiscreteFunctionType& >( discreteFunction ), entity ),
      discreteFunction_( discreteFunction )
  {
  }

  using BaseType::entity;

  using BaseType::operator[];

  RangeFieldType& operator[]( const unsigned int localDofNumber )
  {
    const unsigned int globalDofNumber = discreteFunction_.space().map().toGlobal( entity(), localDofNumber );
    return discreteFunction_[globalDofNumber];
  }

  using BaseType::order;

  using BaseType::size;

  using BaseType::evaluate;

  using BaseType::jacobian;

private:

  DiscreteFunctionType& discreteFunction_;
}; // end class Local

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
