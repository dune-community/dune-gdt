#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH

// system
#include <vector>

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Mapper {

namespace Continuous {

template <class FunctionSpaceImp, class GridPartImp, int polOrder>
class Lagrange
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef GridPartImp GridPartType;

  static const int polynomialOrder = polOrder;

  typedef Lagrange<FunctionSpaceType, GridPartType, polynomialOrder> ThisType;

  typedef LagrangePointSet<GridPartType, polynomialOrder> LagrangePointSetType;

private:
  typedef typename GridPartType::GridType GridType;

  static const int dimension = GridType::dimension;

  typedef LagrangeDiscreteFunctionSpaceTraits<FunctionSpaceType, GridPartType, polynomialOrder> HostTraits;

  typedef typename HostTraits::MapperType MapperType;

  typedef typename HostTraits::BlockMapperType BlockMapperType;

  typedef std::vector<const LagrangePointSetType*> LagrangePointSetContainerType;

  typedef typename GridPartType::IndexSetType IndexSetType;

  typedef LagrangeMapperSingletonKey<GridPartType, LagrangePointSetContainerType> MapperSingletonKeyType;

  typedef LagrangeMapperSingletonFactory<MapperSingletonKeyType, BlockMapperType> BlockMapperSingletonFactoryType;

  typedef SingletonList<MapperSingletonKeyType, BlockMapperType, BlockMapperSingletonFactoryType>
      BlockMapperProviderType;

public:
  typedef typename IndexSetType::IndexType IndexType;

  //! does, whatever the constructor of the dune-fem LagrangeDiscreteFunctionSpace does
  Lagrange(const GridPartType& gridPart)
    : gridPart_(gridPart)
    , lagrangePointSetContainer_(LocalGeometryTypeIndex::size(dimension), nullptr)
    , mapper_(0)
    , blockMapper_(0)
  {
    const IndexSetType& indexSet = gridPart.indexSet();
    const AllGeomTypes<IndexSetType, GridType> allGeometryTypes(indexSet);
    const std::vector<GeometryType>& geometryTypes = allGeometryTypes.geomTypes(0);
    for (unsigned int i = 0; i < geometryTypes.size(); ++i) {
      const GeometryType& gt                        = geometryTypes[i];
      const LagrangePointSetType*& lagrangePointSet = lagrangePointSetContainer_[LocalGeometryTypeIndex::index(gt)];
      if (!lagrangePointSet)
        lagrangePointSet = new LagrangePointSetType(gt, polynomialOrder);
      assert(lagrangePointSet);
    }
    MapperSingletonKeyType key(gridPart, lagrangePointSetContainer_, polynomialOrder);
    blockMapper_ = &BlockMapperProviderType::getObject(key);
    assert(blockMapper_ != 0);
    mapper_ = new MapperType(*blockMapper_);
    assert(mapper_ != 0);
  } // Lagrange(const GridPartType& gridPart)

  //! does, whatever the destructor of the dune-fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
    delete mapper_;
    BlockMapperProviderType::removeObject(*blockMapper_);
    typedef typename LagrangePointSetContainerType::const_iterator IteratorType;
    const IteratorType end = lagrangePointSetContainer_.end();
    for (IteratorType it = lagrangePointSetContainer_.begin(); it != end; ++it) {
      delete *it;
    }
  } // ~Lagrange()

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  template <class EntityType>
  IndexType toGlobal(const EntityType& entity, const IndexType localDofNumber) const
  {
    return mapper_->mapToGlobal(entity, localDofNumber);
  }

  IndexType size() const
  {
    return mapper_->size();
  }

  IndexType maxLocalSize() const
  {
    return mapper_->maxNumDofs();
  }

  template <class EntityType>
  const LagrangePointSetType& lagrangePointSet(const EntityType& entity) const
  {
    const LagrangePointSetType* lagrangePointSet =
        lagrangePointSetContainer_[LocalGeometryTypeIndex::index(entity.type())];
    assert(lagrangePointSet);
    return *lagrangePointSet;
  } // const LagrangePointSetType& lagrangePointSet(const EntityType& entity) const

private:
  Lagrange(const ThisType&);
  ThisType& operator=(const ThisType&);

  const GridPartType& gridPart_;
  mutable LagrangePointSetContainerType lagrangePointSetContainer_;
  MapperType* mapper_;
  BlockMapperType* blockMapper_;
}; // end class Lagrange

} // end namespace Mapper

} // end namespace Continuous

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_CONTINUOUS_LAGRANGE_HH
