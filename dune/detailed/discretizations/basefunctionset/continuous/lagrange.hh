#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// dune-detailed-discretizations includes
#include <dune/detailed/discretizations/basefunctionset/local/lagrange.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace BaseFunctionSet {

namespace Continuous {

template <class DiscreteFunctionSpaceImp>
class Lagrange
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;

  typedef Lagrange<DiscreteFunctionSpaceType> ThisType;

  typedef Dune::Detailed::Discretizations::BaseFunctionSet::Local::Lagrange<ThisType> LocalBaseFunctionSetType;

private:
  typedef typename GridPartType::GridType GridType;

  static const int dimension = GridType::dimension;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<FunctionSpaceType, GridPartType, polynomialOrder> HostTraits;

  typedef typename HostTraits::ShapeFunctionSetType ShapeFunctionSetType;

  typedef Fem::BaseSetLocalKeyStorage<ShapeFunctionSetType> ShapeSetStorageType;

  typedef typename GridPartType::IndexSetType IndexSetType;

  typedef typename HostTraits::BaseFunctionSpaceType BaseFunctionSpaceType;

  typedef LagrangeBaseFunctionFactory<typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension,
                                      polynomialOrder> ScalarFactoryType;

  typedef BaseFunctionSetSingletonFactory<GeometryType, ShapeFunctionSetType, ScalarFactoryType>
      BaseFunctionSetSingletonFactoryType;

  typedef SingletonList<GeometryType, ShapeFunctionSetType, BaseFunctionSetSingletonFactoryType>
      BaseFunctionSetSingletonProviderType;

  typedef typename HostTraits::BaseFunctionSetType BaseFunctionSetType;

public:
  //! does, whatever the constructor of the dune-fem LagrangeDiscreteFunctionSpace does
  Lagrange(const DiscreteFunctionSpaceType& space)
    : space_(space)
    , shapeFunctionSets_()
  {
    const IndexSetType& indexSet = space.gridPart().indexSet();
    const AllGeomTypes<IndexSetType, GridType> allGeometryTypes(indexSet);
    const std::vector<GeometryType>& geometryTypes = allGeometryTypes.geomTypes(0);
    for (unsigned int i = 0; i < geometryTypes.size(); ++i) {
      const GeometryType& gt = geometryTypes[i];
      shapeFunctionSets_.template insert<BaseFunctionSetSingletonProviderType>(gt);
    }
  } // Lagrange(const DiscreteFunctionSpaceType& space)

  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  template <class EntityType>
  const LocalBaseFunctionSetType local(const EntityType& entity) const
  {
    return LocalBaseFunctionSetType(*this, entity);
  }

private:
  Lagrange(const ThisType&);
  ThisType& operator=(const ThisType&);

  template <class EntityType>
  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(&shapeFunctionSets_[entity.type()]);
  }

  friend class Dune::Detailed::Discretizations::BaseFunctionSet::Local::Lagrange<ThisType>;

  const DiscreteFunctionSpaceType& space_;
  mutable ShapeSetStorageType shapeFunctionSets_;
}; // end class Lagrange

} // end namespace BaseFunctionSet

} // end namespace Continuous

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH
