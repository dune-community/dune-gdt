#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class SpaceTraits >
class FemLocalfunctionsWrapper;


template< class SpaceTraits >
class FemLocalfunctionsWrapperTraits
{
public:
  typedef FemLocalfunctionsWrapper< SpaceTraits >                                         derived_type;
  typedef typename SpaceTraits::derived_type                                              SpaceType;
  typedef typename SpaceTraits::BackendType::BaseFunctionSetMapType::BaseFunctionSetType  BackendType;
  typedef typename SpaceTraits::MapperType::IndexType                                     IndexType;

  typedef typename BackendType::EntityType        EntityType;
  typedef typename BackendType::DomainType        DomainType;
  typedef typename BackendType::RangeType         RangeType;
  typedef typename BackendType::JacobianRangeType JacobianRangeType;
};


template< class SpaceTraits >
class FemLocalfunctionsWrapper
  : public Interface< FemLocalfunctionsWrapperTraits< SpaceTraits > >
{
public:
  typedef FemLocalfunctionsWrapper< SpaceTraits >       ThisType;
  typedef FemLocalfunctionsWrapperTraits< SpaceTraits > Traits;
  typedef Interface< Traits >                           InterfaceType;

  typedef typename Traits::SpaceType    SpaceType;
  typedef typename Traits::BackendType  BackendType;

  typedef typename Traits::IndexType          IndexType;
  typedef typename Traits::EntityType         EntityType;
  typedef typename Traits::DomainType         DomainType;
  typedef typename Traits::RangeType          RangeType;
  typedef typename Traits::JacobianRangeType  JacobianRangeType;

  FemLocalfunctionsWrapper(const SpaceType& _space, const EntityType& _entity)
    : space_(_space)
    , entity_(_entity)
    , backend_(space_.backend().baseFunctionSetMap().find(entity_))
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  IndexType size() const
  {
    return backend_.size();
  }

  unsigned int order() const
  {
    return space_.backend().baseFunctionSetMap().getOrder(entity_);
  }

  void evaluate(const DomainType& x, std::vector< RangeType >& ret) const
  {
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector< JacobianRangeType >& ret) const
  {
    backend_.jacobianAll(x, entity_.geometry().jacobianInverseTransposed(x), ret);
  }

private:
  const SpaceType& space_;
  const EntityType& entity_;
  const BackendType backend_;
}; // class FemLocalfunctionsWrapper


} // namespace BaseFunctionSet
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
