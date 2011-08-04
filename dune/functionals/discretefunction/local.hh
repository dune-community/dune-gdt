#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunction
{

template< class DiscreteFunctionImp >
class Local
{
public:

  typedef DiscreteFunctionImp
    DiscreteFunctionType;

  typedef Local< DiscreteFunctionType >
    ThisType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

//  typedef typename DiscreteFunctionSpaceType::LocalBaseFunctionSetType
//    LocalBaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType::template Codim< 0 >::IteratorType::Entity
    EntityType;

//  typedef typename DiscreteFunctionType::StorageType
//    StorageType;

  typedef typename DiscreteFunctionType::DomainType
    DomainType;

  typedef typename DiscreteFunctionType::RangeFieldType
    RangeFieldType;

  typedef typename DiscreteFunctionType::RangeType
    RangeType;

  typedef typename DiscreteFunctionType::JacobianRangeType
    JacobianRangeType;

  Local( DiscreteFunctionType& discreteFunction, const EntityType& entity )
    : discreteFunction_( discreteFunction ),
      entity_( entity ),
      size_( discreteFunction_.space().localBaseFunctionSet( entity ).size() ),
      order_( discreteFunction_.space().localBaseFunctionSet( entity ).order() )
  {
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  RangeFieldType& operator[] ( const unsigned int localDofNumber )
  {
    const unsigned int globalDofNumber = discreteFunction_.space().map().toGlobal( entity_, localDofNumber );
    return discreteFunction_[globalDofNumber];
  }

  const RangeFieldType& operator[] ( const unsigned int localDofNumber ) const
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
    std::vector< RangeType > baseFunctionValues( 0.0 );
    discreteFunction_.space().localBaseFunctionSet( entity_ ).evaluate( x, baseFunctionValues );
    ret = 0.0;
    for( unsigned int i = 0; i < size(); ++i )
    {
      baseFunctionValues[i] *= operator[](i);
      ret += baseFunctionValues[i];
    }
  }

  void jacobian( const DomainType& x, JacobianRangeType& ret ) const
  {
    std::vector< JacobianRangeType > baseFunctionJacobianValues( 0.0 );
    discreteFunction_.space().localBaseFunctionSet( entity_ ).jacobian( x, baseFunctionJacobianValues );
    ret = 0.0;
    for( unsigned int i = 0; i < size(); ++i )
    {
      baseFunctionJacobianValues[i] *= operator[](i);
      ret += baseFunctionJacobianValues[i];
    }
  }

private:

  DiscreteFunctionType& discreteFunction_;
  const EntityType& entity_;
  const unsigned int size_;
  const int order_;

}; // end class Local

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

//  typedef typename DiscreteFunctionSpaceType::LocalBaseFunctionSetType
//    LocalBaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType::template Codim< 0 >::IteratorType::Entity
    EntityType;

//  typedef typename DiscreteFunctionType::StorageType
//    StorageType;

  typedef typename DiscreteFunctionType::DomainType
    DomainType;

  typedef typename DiscreteFunctionType::RangeFieldType
    RangeFieldType;

  typedef typename DiscreteFunctionType::RangeType
    RangeType;

  typedef typename DiscreteFunctionType::JacobianRangeType
    JacobianRangeType;

  LocalConst( const DiscreteFunctionType& discreteFunction, const EntityType& entity )
    : discreteFunction_( discreteFunction ),
      entity_( entity ),
      size_( discreteFunction_.space().localBaseFunctionSet( entity ).size() ),
      order_( discreteFunction_.space().localBaseFunctionSet( entity ).order() )
  {
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const RangeFieldType& operator[] ( const unsigned int localDofNumber ) const
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
    std::vector< RangeType > baseFunctionValues( 0.0 );
    discreteFunction_.space().localBaseFunctionSet( entity_ ).evaluate( x, baseFunctionValues );
    ret = 0.0;
    for( unsigned int i = 0; i < size(); ++i )
    {
      baseFunctionValues[i] *= operator[](i);
      ret += baseFunctionValues[i];
    }
  }

  void jacobian( const DomainType& x, JacobianRangeType& ret ) const
  {
    std::vector< JacobianRangeType > baseFunctionJacobianValues( 0.0 );
    discreteFunction_.space().localBaseFunctionSet( entity_ ).jacobian( x, baseFunctionJacobianValues );
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
  const unsigned int size_;
  const int order_;

}; // end class LocalConst

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
