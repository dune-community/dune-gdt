#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

//// dune-functionals includes
//#include <dune/functionals/discretefunctionspace/continuous/lagrangefemadapter.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace Mapper
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

  enum{ polynomialOrder = polOrder };

  typedef Lagrange< FunctionSpaceType, GridPartType, polynomialOrder >
    ThisType;

  typedef LagrangePointSet< GridPartType, polynomialOrder >
    LagrangePointSetType;

private:

  typedef typename GridPartType::GridType
    GridType;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                      GridPartType,
                                                      polynomialOrder >
    LagrangeDiscreteFunctionSpaceTraitsType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::MapperType
    MapperType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BlockMapperType
    BlockMapperType;

  typedef std::map< const GeometryType, const LagrangePointSetType* >
    LagrangePointSetMapType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::IndexSetType
    IndexSetType;

  typedef LagrangeMapperSingletonKey< GridPartType, LagrangePointSetMapType >
    MapperSingletonKeyType;

  typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, BlockMapperType >
    BlockMapperSingletonFactoryType;

  typedef SingletonList< MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType >
    BlockMapperProviderType;

public:

  //! does, whatever the constructor of the fem LagrangeDiscreteFunctionSpace does
  Lagrange( const GridPartType& gridPart )
    : gridPart_( gridPart ),
      lagrangePointSet_(),
      mapper_( 0 ),
      blockMapper_( 0 )
  {
    const IndexSetType& indexSet = gridPart_.indexSet();

    AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );

    const std::vector< GeometryType >& geometryTypes = allGeometryTypes.geomTypes( 0 );

    for( unsigned int i = 0; i < geometryTypes.size(); ++i )
    {
      const GeometryType& geometryType = geometryTypes[ i ];

      if( lagrangePointSet_.find( geometryType ) == lagrangePointSet_.end() )
      {
        const LagrangePointSetType* lagrangePointSet = new LagrangePointSetType(geometryType);
        assert( lagrangePointSet != NULL );
        lagrangePointSet_[ geometryType ] = lagrangePointSet;
      }
    }

    MapperSingletonKeyType key( gridPart_, lagrangePointSet_, polynomialOrder );

    blockMapper_ = &( BlockMapperProviderType::getObject( key ) );
    assert( blockMapper_ != 0 );

    mapper_ = new MapperType( *blockMapper_ );
    assert( mapper_ != 0 );
  }

  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
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

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  template< class EntityType >
  unsigned int toGlobal( const EntityType& entity, const unsigned int localDofNumber ) const
  {
    return mapper_->mapToGlobal( entity, localDofNumber );
  }

  unsigned int size() const
  {
    return mapper_->size();
  }

  unsigned int maxLocalSize() const
  {
    return mapper_->maxNumDofs();
  }

  template< class EntityType >
  const LagrangePointSetType& lagrangePointSet( const EntityType& entity ) const
  {
    assert( lagrangePointSet_.find( entity.type() ) != lagrangePointSet_.end() );
    assert( lagrangePointSet_[entity.type()] != NULL );
    return *( lagrangePointSet_[entity.type()] );
  }

private:
  //! copy constructor
  Lagrange( const ThisType& );

  //! assignment operator
  ThisType& operator=( const ThisType& );

//  template< class >
//  friend class Dune::Functionals::DiscreteFunctionSpace::Continuous::LagrangeFemAdapter;

  const GridPartType& gridPart_;
  mutable LagrangePointSetMapType lagrangePointSet_;
  MapperType* mapper_;
  BlockMapperType* blockMapper_;

}; // end class Lagrange

} // end namespace Mapper

} // end namespace Continuous

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH
