#ifndef DUNE_FUNCTIONALS_MAPPER_LAGRANGE_HH
#define DUNE_FUNCTIONALS_MAPPER_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>


namespace Dune
{

namespace Functionals
{

namespace Mapper
{

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
class Lagrange
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef GridPartImp
    GridPartType;

  static const unsigned int order = polOrder;

  typedef Lagrange< FunctionSpaceType, GridPartType, order >
    ThisType;

private:

  typedef typename GridPartType::GridType
    GridType;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                      GridPartType,
                                                      order >
    LagrangeDiscreteFunctionSpaceTraitsType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::MapperType
    MapperType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BlockMapperType
    BlockMapperType;

  typedef LagrangePointSet< GridPartType, order >
    LagrangePointSetType;

  typedef std::map< const GeometryType, const LagrangePointSetType* >
    LagrangePointSetMapType;

public:

  //! does, whatever the constructor of the fem LagrangeDiscreteFunctionSpace does
  Lagrange( const GridPartType& gridPart )
    : gridPart_( gridPart ),
      lagrangePointSet_(),
      mapper_( 0 ),
      blockMapper_( 0 )
  {
    typedef typename LagrangeDiscreteFunctionSpaceTraitsType::IndexSetType
      IndexSetType;

    typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetMapType >
      MapperSingletonKeyType;

    typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetMapType >
      MapperSingletonKeyType;

    typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
      BlockMapperSingletonFactoryType;

    typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

    const IndexSetType& indexSet = gridPart_.indexSet();

    AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );

    const std::vector< GeometryType >& geometryTypes = allGeometryTypes.geomTypes( 0 );

    for( unsigned int i = 0; i < geometryTypes.size(); ++i )
    {
      const GeometryType& geometryType = geometryTypes[ i ];

      if( lagrangePointSet_.find( geometryType ) == lagrangePointSet_.end() )
      {
        const LagrangePointSetType* lagrangePointSet = new LagrangePointSetType( geometryType );
        assert( lagrangePointSet != NULL );
        lagrangePointSet_[ geometryType ] = lagrangePointSet;
      }
    }

    MapperSingletonKeyType key( gridPart_, lagrangePointSet_, order );

    blockMapper_ = &( BlockMapperProviderType::getObject( key ) );
    assert( blockMapper_ != 0 );

    mapper_ = new MapperType( *blockMapper_ );
    assert( mapper_ != 0 );
  }

  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
    typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetMapType >
      MapperSingletonKeyType;

    typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
      BlockMapperSingletonFactoryType;

    typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

    delete mapper_;

    BlockMapperProviderType::removeObject( *blockMapper_ );

    typedef typename LagrangePointSetMapType::iterator LPIteratorType;

    LPIteratorType lpend = lagrangePointSet_.end();
    for( LPIteratorType it = lagrangePointSet_.begin(); it != lpend; ++it )
    {
      const LagrangePointSetType *lagrangePointSet = (*it).second;
      if( lagrangePointSet != NULL )
        delete lagrangePointSet;
    }
  }

  template< class EntityType >
  unsigned int toGlobal( const EntityType& entity, const unsigned int localDofNumber )
  {
    return mapper_->mapToGlobal( entity, localDofNumber );
  }

private:

  const GridPartType& gridPart_;
  mutable LagrangePointSetMapType lagrangePointSet_;
  MapperType* mapper_;
  BlockMapperType* blockMapper_;

}; // end class Lagrange

} // end namespace Mapper

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_MAPPER_LAGRANGE_HH
