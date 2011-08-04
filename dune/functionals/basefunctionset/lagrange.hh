#ifndef DUNE_FUNCTIONALS_BASEFUNCTIONSET_LAGRANGE_HH
#define DUNE_FUNCTIONALS_BASEFUNCTIONSET_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// dune-functionals includes
#include <dune/functionals/basefunctionset/local/lagrange.hh>
//#include <dune/functionals/discretefunctionspace/continuous/lagrangefemadapter.hh>

//// dune-fem-tools includes
//#include <dune/fem-tools/common/printing.hh>
//#include <dune/fem-tools/common/string.hh>
//#include <dune/fem-tools/grid/entity.hh>

namespace Dune
{

namespace Functionals
{

namespace BaseFunctionSet
{

template< class DiscreteFunctionSpaceImp >
class Lagrange
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  enum{ polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

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

  typedef Dune::Functionals::BaseFunctionSet::Local::Lagrange< ThisType >
    LocalBaseFunctionSetType;

private:

  typedef typename GridPartType::GridType
    GridType;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                      GridPartType,
                                                      polynomialOrder >
    LagrangeDiscreteFunctionSpaceTraitsType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetImp
    BaseFunctionSetImp;

  typedef std::map< const GeometryType, const BaseFunctionSetImp* >
    HostBaseFunctionMapType;

  typedef typename GridPartType::template Codim< 0 >::IteratorType
    EntityIteratorType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::IndexSetType
    IndexSetType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSpaceType
    BaseFunctionSpaceType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetType
    BaseFunctionSetType;

  enum{ dimension = GridType::dimension };

  typedef LagrangeBaseFunctionFactory< typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension, polynomialOrder >
    ScalarFactoryType;

  typedef BaseFunctionSetSingletonFactory < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
    BaseFunctionSetSingletonFactoryType;

  typedef SingletonList< GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
    BaseFunctionSetSingletonProviderType;

public:

  //! does, whatever the constructor of the fem LagrangeDiscreteFunctionSpace does
  Lagrange( const DiscreteFunctionSpaceType& space )
    : space_( space ),
      hostBaseFunctionSetMap_()
  {
    std::cout << "BaseFunctionSet::Lagrange::Lagrange(){" << std::endl;
    mapInspector( "  " );
    const IndexSetType& indexSet = space_.gridPart().indexSet();

    AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );

    const std::vector< GeometryType >& geometryTypes = allGeometryTypes.geomTypes( 0 );

    for( unsigned int i = 0; i < geometryTypes.size(); ++i )
    {
      const GeometryType& geometryType = geometryTypes[ i ];

      if( hostBaseFunctionSetMap_.find( geometryType ) == hostBaseFunctionSetMap_.end() )
      {
        const BaseFunctionSetImp* baseFunctionSet = &( BaseFunctionSetSingletonProviderType::getObject( geometryType ) );
        assert( baseFunctionSet != NULL );

        hostBaseFunctionSetMap_[ geometryType ] = baseFunctionSet;
      }
    }
    mapInspector( "  " );
    std::cout << "}" << std::endl;
  }

  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
    std::cout << "BaseFunctionSet::Lagrange::~Lagrange(){" << std::endl;
    mapInspector( "  " );
    typedef typename HostBaseFunctionMapType::iterator
      BFIteratorType;
    BFIteratorType bfend = hostBaseFunctionSetMap_.end();
    for( BFIteratorType it = hostBaseFunctionSetMap_.begin(); it != bfend; ++it )
    {
      const BaseFunctionSetImp* baseFunctionSet = (*it).second;
      if( baseFunctionSet != NULL )
        BaseFunctionSetSingletonProviderType::removeObject( *baseFunctionSet );
    }
    mapInspector( "  " );
    std::cout << "}" << std::endl;
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  template< class EntityType >
  LocalBaseFunctionSetType local( const EntityType& entity ) const
  {
    return LocalBaseFunctionSetType( *this, entity );
  }

private:

  void mapInspector( const std::string prefix = "" ) const
  {
    std::cout << prefix << "BaseFunctionSet::Lagrange::mapInspector()" << std::endl;
    if( hostBaseFunctionSetMap_.empty() )
      std::cout << prefix << "  map is empty!" << std::endl;
    else
    {
      const unsigned int size = hostBaseFunctionSetMap_.size();
      std::cout << prefix << "  map has " << size << " element";
      if( size > 1 )
        std::cout << "s";
      std::cout << "!" << std::endl;
    }
  }

  template< class EntityType >
  BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
  {
    std::cout << "BaseFunctionSet::Lagrange::baseFunctionSet(){" << std::endl;
    mapInspector( "  " );
    // get the basefunctionset
    assert( hostBaseFunctionSetMap_.find( entity.type() ) != hostBaseFunctionSetMap_.end() );
    assert( hostBaseFunctionSetMap_[ entity.type() ] != NULL );
    std::cout << "}" << std::endl;
    return BaseFunctionSetType( hostBaseFunctionSetMap_[ entity.type() ] );
  }

//  friend class Dune::Functionals::DiscreteFunctionSpace::Continuous::LagrangeFemAdapter< DiscreteFunctionSpaceType >;
  friend class Dune::Functionals::BaseFunctionSet::Local::Lagrange< ThisType >;

  const DiscreteFunctionSpaceType& space_;
  mutable HostBaseFunctionMapType hostBaseFunctionSetMap_;

}; // end class Lagrange

} // end namespace BaseFunctionSet

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_BASEFUNCTIONSET_LAGRANGE_HH
