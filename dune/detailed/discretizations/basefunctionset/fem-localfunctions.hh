#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits and to allow for specialization
template <class BaseFunctionSetMapImp>
class BaseFunctionSetFemLocalfunctionsWrapper;


template <class BaseFunctionSetMapImp>
class BaseFunctionSetFemLocalfunctionsWrapperTraits
{
public:
  typedef BaseFunctionSetFemLocalfunctionsWrapper<BaseFunctionSetMapImp> derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType BackendType;
  typedef typename BackendType::EntityType EntityType;
  typedef typename BackendType::DomainType DomainType;
  typedef typename BackendType::RangeType RangeType;
  typedef typename BackendType::JacobianRangeType JacobianRangeType;
};


template <class BaseFunctionSetMapImp>
class BaseFunctionSetFemLocalfunctionsWrapper
    : public BaseFunctionSetInterface<BaseFunctionSetFemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp>>
{
public:
  typedef BaseFunctionSetFemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::JacobianRangeType JacobianRangeType;

  BaseFunctionSetFemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& en)
    : baseFunctionSetMap_(baseFunctionSetMap)
    , entity_(en)
    , backend_(baseFunctionSetMap_.find(entity_))
  {
  }

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

  void evaluate(const DomainType& x, std::vector<RangeType>& ret) const
  {
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    backend_.jacobianAll(x, entity_.geometry().jacobianInverseTransposed(x), ret);
  }

private:
  const BaseFunctionSetMapImp& baseFunctionSetMap_;
  const EntityType& entity_;
  const BackendType backend_;
}; // class BaseFunctionSetFemLocalfunctionsWrapper


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
