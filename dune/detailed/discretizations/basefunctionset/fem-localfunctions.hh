#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class BaseFunctionSetMapImp >
class FemLocalfunctionsWrapper;


template< class BaseFunctionSetMapImp >
class FemLocalfunctionsWrapperTraits
{
public:
  typedef FemLocalfunctionsWrapper< BaseFunctionSetMapImp >  derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType               BackendType;
  typedef typename BackendType::EntityType        EntityType;
  typedef typename BackendType::DomainFieldType   DomainFieldType;
  static const unsigned int                       dimDomain = BackendType::dimDomain;
  typedef typename BackendType::DomainType        DomainType;
  typedef typename BackendType::RangeFieldType    RangeFieldType;
  static const unsigned int                       dimRange = BackendType::dimRange;
  typedef typename BackendType::RangeType         RangeType;
  typedef typename BackendType::JacobianRangeType JacobianRangeType;
};


template< class BaseFunctionSetMapImp >
class FemLocalfunctionsWrapper
  : public BaseFunctionSetInterface< FemLocalfunctionsWrapperTraits< BaseFunctionSetMapImp > >
{
public:
  typedef FemLocalfunctionsWrapperTraits< BaseFunctionSetMapImp > Traits;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;
  typedef typename Traits::DomainFieldType  DomainFieldType;
  static const unsigned int                 dimDomain = Traits::dimDomain;
  typedef typename Traits::DomainType       DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  typedef typename Traits::RangeType      RangeType;
  typedef typename Traits::JacobianRangeType  JacobianRangeType;

  FemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& en)
    : baseFunctionSetMap_(baseFunctionSetMap)
    , entity_(en)
    , backend_(baseFunctionSetMap_.find(entity_))
  {}

  const EntityType& entity() const
  {
    return entity_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  size_t order() const
  {
    return baseFunctionSetMap_.getOrder(entity_);
  }

  void evaluate(const DomainType& x, std::vector< RangeType >& ret) const
  {
    assert(ret.size() >= size());
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector< JacobianRangeType >& ret) const
  {
    assert(ret.size() >= size());
    backend_.jacobianAll(x, entity_.geometry().jacobianInverseTransposed(x), ret);
  }

private:
  const BaseFunctionSetMapImp& baseFunctionSetMap_;
  const EntityType& entity_;
  const BackendType backend_;
}; // class FemLocalfunctionsWrapper


} // namespace BaseFunctionSet
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
