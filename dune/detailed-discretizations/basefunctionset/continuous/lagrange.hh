#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/basefunctionset/local/lagrange.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace BaseFunctionSet {

namespace Continuous {

template <class DiscreteFunctionSpaceImp>
class Lagrange
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef Lagrange<DiscreteFunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef Dune::DetailedDiscretizations::BaseFunctionSet::Local::Lagrange<ThisType> LocalBaseFunctionSetType;

private:
  typedef typename GridPartType::GridType GridType;

  typedef Dune::LagrangeDiscreteFunctionSpaceTraits<FunctionSpaceType, GridPartType, polynomialOrder>
      LagrangeDiscreteFunctionSpaceTraitsType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetImp BaseFunctionSetImp;

  typedef std::map<const GeometryType, const BaseFunctionSetImp*> HostBaseFunctionMapType;

  class FieldVectorComp
  {
  public:
    bool operator()(const DomainType& lhs, const DomainType& rhs) const
    {
      DomainFieldType lhsNorm = lhs.two_norm();
      DomainFieldType rhsNorm = rhs.two_norm();
      if (lhsNorm < rhsNorm) {
        return true;
      } else {
        if (lhsNorm > rhsNorm) {
          return false;
        } else {
          return (lhs[0] < rhs[0]);
        }
      }
    }
  };

  //  typedef /*std::map< GeometryType,*/ std::map< DomainType, std::vector< RangeType >, FieldVectorComp > /*>*/
  //  EvaluationCacheMapType;

  //  typedef /*std::map< GeometryType,*/ std::map< DomainType, std::vector< JacobianRangeType >, FieldVectorComp >
  //  /*>*/ JacobianCacheMapType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::IndexSetType IndexSetType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSpaceType BaseFunctionSpaceType;

  typedef typename LagrangeDiscreteFunctionSpaceTraitsType::BaseFunctionSetType BaseFunctionSetType;

  enum
  {
    dimension = GridType::dimension
  };

  typedef LagrangeBaseFunctionFactory<typename BaseFunctionSpaceType::ScalarFunctionSpaceType, dimension,
                                      polynomialOrder> ScalarFactoryType;

  typedef BaseFunctionSetSingletonFactory<GeometryType, BaseFunctionSetImp, ScalarFactoryType>
      BaseFunctionSetSingletonFactoryType;

  typedef SingletonList<GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType>
      BaseFunctionSetSingletonProviderType;

public:
  //! does, whatever the constructor of the fem LagrangeDiscreteFunctionSpace does
  Lagrange(const DiscreteFunctionSpaceType& space)
    : space_(space)
    , hostBaseFunctionSetMap_() /*,
       evaluationCacheMap_(),
       jacobianCacheMap_()*/
  {
    const IndexSetType& indexSet = space_.gridPart().indexSet();

    const AllGeomTypes<IndexSetType, GridType> allGeometryTypes(indexSet);

    const std::vector<GeometryType>& geometryTypes = allGeometryTypes.geomTypes(0);

    for (unsigned int i = 0; i < geometryTypes.size(); ++i) {
      const GeometryType& geometryType = geometryTypes[i];

      if (hostBaseFunctionSetMap_.find(geometryType) == hostBaseFunctionSetMap_.end()) {
        const BaseFunctionSetImp* baseFunctionSet = &(BaseFunctionSetSingletonProviderType::getObject(geometryType));
        assert(baseFunctionSet != NULL);

        hostBaseFunctionSetMap_[geometryType] = baseFunctionSet;
      }
    }
  } // end constructor

private:
  //! copy constructor
  Lagrange(const ThisType& other);
  //    : space_( other.space() ),
  //      hostBaseFunctionSetMap_()
  //  {
  //    if( !other.hostBaseFunctionSetMap_.empty() )
  //    {
  //      for( unsigned int i = 0; i < other.hostBaseFunctionSetMap_.size(); ++i )
  //      {
  //        hostBaseFunctionSetMap_[i] = other.hostBaseFunctionSetMap_[i];
  //      }
  //    }
  //  }

public:
  //! does, whatever the destructor of the fem LagrangeDiscreteFunctionSpace does
  ~Lagrange()
  {
    typedef typename HostBaseFunctionMapType::iterator BFIteratorType;
    BFIteratorType bfend = hostBaseFunctionSetMap_.end();
    for (BFIteratorType it = hostBaseFunctionSetMap_.begin(); it != bfend; ++it) {
      const BaseFunctionSetImp* baseFunctionSet = (*it).second;
      if (baseFunctionSet != NULL)
        BaseFunctionSetSingletonProviderType::removeObject(*baseFunctionSet);
    }
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
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class EntityType>
  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    // get the basefunctionset
    assert(hostBaseFunctionSetMap_.find(entity.type()) != hostBaseFunctionSetMap_.end());
    assert(hostBaseFunctionSetMap_[entity.type()] != NULL);
    return BaseFunctionSetType(hostBaseFunctionSetMap_[entity.type()]);
  }

  friend class Dune::DetailedDiscretizations::BaseFunctionSet::Local::Lagrange<ThisType>;

  //  const EvaluationCacheMapType& evaluationCacheMap() const
  //  {
  //    return evaluationCacheMap_;
  //  }

  //  const JacobianCacheMapType& jacobianCacheMap() const
  //  {
  //    return jacobianCacheMap_;
  //  }

  const DiscreteFunctionSpaceType& space_;
  mutable HostBaseFunctionMapType hostBaseFunctionSetMap_;
  //  mutable EvaluationCacheMapType evaluationCacheMap_;
  //  mutable JacobianCacheMapType jacobianCacheMap_;
}; // end class Lagrange

} // end namespace BaseFunctionSet

} // end namespace Continuous

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONTINUOUS_BASEFUNCTIONSET_LAGRANGE_HH
