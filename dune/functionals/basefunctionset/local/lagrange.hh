#ifndef DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
#define DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

namespace Dune
{

namespace Functionals
{

namespace BaseFunctionSet
{

namespace Local
{

template< class DiscreteFunctionSpaceImp >
class Lagrange
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FucntionSpaceType
    FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  static const unsigned int order = DiscreteFunctionSpaceType::order;

  typedef Lagrange< DiscreteFunctionSpaceType >
    ThisType;

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

  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;

private:

  typedef typename GridPartType::GridType
    GridType;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                      GridPartType,
                                                      order >
    LagrangeDiscreteFunctionSpaceTraitsType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetImp
    BaseFunctionSetImp;

  typedef std::map< const GeometryType, const BaseFunctionSetImp* >
    BaseFunctionMapType;

public:

  //! does, whatever the constructor of the fem LagrangeDiscreteFunctionSpace does
  Lagrange( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    : space_( space ),
      entity_( entity ),
      baseFunctionSet_()
  {
    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::IndexSetType
      IndexSetType;

    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSpaceType
      BaseFunctionSpaceType;

    enum{ dimension = GridType::dimension };

    typedef LagrangeBaseFunctionFactory< typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension, order >
      ScalarFactoryType;

    typedef BaseFunctionSetSingletonFactory < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
      BaseFunctionSetSingletonFactoryType;

    typedef SingletonList< GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
      BaseFunctionSetSingletonProviderType;

    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetType
      BaseFunctionSetType;

    const IndexSetType& indexSet = space_.gridPart().indexSet();

    AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );

    const std::vector< GeometryType >& geometryTypes = allGeometryTypes.geomTypes( 0 );

    for( unsigned int i = 0; i < geometryTypes.size(); ++i )
    {
      const GeometryType& geometryType = geometryTypes[ i ];

      if( baseFunctionSet_.find( geometryType ) == baseFunctionSet_.end() )
      {
        const BaseFunctionSetImp* baseFunctionSet = &( BaseFunctionSetSingletonProviderType::getObject( geometryType ) );
        assert( baseFunctionSet != NULL );

        baseFunctionSet_[ geometryType ] = baseFunctionSet;
      }
    }

    assert( baseFunctionSet_.find( entity.type() ) != baseFunctionSet_.end() );
    assert( baseFunctionSet_[ entity.type() ] != NULL );

    BaseFunctionSetType baseFunctionSet( baseFunctionSet_[ entity.type() ] );

    size_ = baseFunctionSet.numBaseFunctions();
  }

  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSpaceType
      BaseFunctionSpaceType;
    enum{ dimension = GridType::dimension };
    typedef LagrangeBaseFunctionFactory< typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension, order >
      ScalarFactoryType;
    typedef BaseFunctionSetSingletonFactory < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
      BaseFunctionSetSingletonFactoryType;
    typedef SingletonList< GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
      BaseFunctionSetSingletonProviderType;
    typedef typename BaseFunctionMapType::iterator
      BFIteratorType;
    BFIteratorType bfend = baseFunctionSet_.end();
    for( BFIteratorType it = baseFunctionSet_.begin(); it != bfend; ++it )
    {
      const BaseFunctionSetImp* baseFunctionSet = (*it).second;
      if( baseFunctionSet != NULL )
        BaseFunctionSetSingletonProviderType::removeObject( *baseFunctionSet );
    }
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  unsigned int size() const
  {
    return size_;
  }

  void evaluate( const DomainType& x, std::vector< RangeType >& ret)
  {
    // get the basefunctionset
    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetType
      BaseFunctionSetType;
    BaseFunctionSetType baseFunctionSet( baseFunctionSet_[ entity_.type() ] );

    // and evaluate
    for( unsigned int i = 0; i < size_; ++i )
    {
      baseFunctionSet.evaluate( i, x, ret[i] );
    }
  }

  void jacobian( const DomainType& x, std::vector< JacobianRangeType >& ret ) const
  {
    // some types we will need
    typedef typename EntityType::Geometry
      EntityGeometryType;
    typedef typename EntityGeometryType::Jacobian
      JacobianInverseTransposedType;
    typedef typename JacobianRangeType::row_type
      JacobianRowType;

    // get the basefunctionset
    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetType
      BaseFunctionSetType;
    BaseFunctionSetType baseFunctionSet( baseFunctionSet_[ entity_.type() ] );

    // geometry and jacobian inverse transposed
    const EntityGeometryType& entityGeometry = entity_.geometry();
    const JacobianInverseTransposedType& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed( x );

    // some tmp storage
    JacobianRangeType jacobianUntransposed( 0.0 );
    JacobianRowType jacobian( 0.0 );

    // evaluate
    for( unsigned int i = 0; i < size_; ++i )
    {
      // get untransposed jacobian
      baseFunctionSet.jacobian( i, x, jacobianUntransposed );

      // transpose for each dim of range
      for( unsigned int row = 0; row < ret.N(); ++row )
      {
        // transpose
        jacobianInverseTransposed.mv( jacobianUntransposed[0], jacobian );

        // return
        ret[i][row] = jacobian;
      }
    }
  }

private:

  const DiscreteFunctionSpaceType& space_;
  const EntityType& entity_;
  mutable BaseFunctionMapType baseFunctionSet_;
  unsigned int size_;

}; // end class LocalBaseFunctionSet

//template< class InducingDiscreteFunctionImp >
//class LocalBaseFunctionSetWrapper
//{
//public:

//  typedef InducingDiscreteFunctionImp
//    InducingDiscreteFunctionType;

//  typedef typename InducingDiscreteFunctionType::LocalFunctionType
//    LocalBaseFunctionType;

//  typedef typename InducingDiscreteFunctionType::DiscreteFunctionSpaceType
//    DiscreteFunctionSpaceType;

//  typedef typename DiscreteFunctionSpaceType::EntityType
//    EntityType;

//  typedef typename DiscreteFunctionSpaceType::DomainFieldType
//    DomainFieldType;

//  typedef typename DiscreteFunctionSpaceType::DomainType
//    DomainType;

//  typedef typename DiscreteFunctionSpaceType::RangeFieldType
//    RangeFieldType;

//  typedef typename DiscreteFunctionSpaceType::RangeType
//    RangeType;

//  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
//    JacobianRangeType;

//  typedef typename DiscreteFunctionSpaceType::HessianRangeType
//    HessianRangeType;

//private:

//  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
//    HostBaseFunctionSetType;

//public:

//  LocalBaseFunctionSetWrapper( const InducingDiscreteFunctionType& inducingDiscreteFunction, const EntityType& entity )
//    : inducingDiscreteFunction_( inducingDiscreteFunction ),
//      entity_( entity )
//  {
//  }

//  const DiscreteFunctionSpaceType& space() const
//  {
//    return inducingDiscreteFunction_.space();
//  }

//  const EntityType& entity() const
//  {
//    return entity_;
//  }

//  const LocalBaseFunctionType baseFunction( const int i ) const
//  {
//    assert( i < numBaseFunctions() );
//    return inducingDiscreteFunction_.localFunction( entity );
//  }

//  const int order() const
//  {
//    return inducingDiscreteFunction_.space().order();
//  }

//  const int numBaseFunctions() const
//  {
//    return 1;
//  }

//  void evaluate( const int i, const DomainType& x, RangeType& ret ) const
//  {
//    assert( i < numBaseFunctions() );
//    const LocalBaseFunctionType localFunction = inducingDiscreteFunction_.localFunction( entity );
//    localFunction.evaluate( x, ret );
//  }

//  /**
//    \brief      evaluates the jacobian of the ith ocal basefunction
//    \attention  the evalaution is already multiplied by entityGeometry.jacobianInverseTransposed( x )
//    **/
//  void jacobian( const int i, const DomainType& x, JacobianRangeType& ret ) const
//  {
//    assert( i < numBaseFunctions() );
//    const LocalBaseFunctionType localFunction = inducingDiscreteFunction_.localFunction( entity );
//    localFunction.jacobian( x, ret );
//  }

//private:

//  const InducingDiscreteFunctionType& inducingDiscreteFunction_;
//  const EntityType& entity_;

//}; // end class LocalBaseFunctionSetWrapper

} // end namespace Local

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
