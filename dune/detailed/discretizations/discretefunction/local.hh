#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_LOCAL_HH

#include <vector>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace DiscreteFunction {

template <class DiscreteFunctionImp>
class LocalConst
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef LocalConst<DiscreteFunctionType> ThisType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType::template Codim<0>::EntityType EntityType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionType::size_type size_type;

private:
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalBaseFunctionSetType;

public:
  LocalConst(const DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : discreteFunction_(discreteFunction)
    , entity_(entity)
    , localBaseFunctionSet_(discreteFunction_.space().baseFunctionSet().local(entity_))
    , size_(localBaseFunctionSet_.size())
    , order_(localBaseFunctionSet_.order())
  {
  }

  LocalConst(const ThisType& other)
    : discreteFunction_(other.discreteFunction_)
    , entity_(other.entity_)
    , localBaseFunctionSet_(other.localBaseFunctionSet_)
    , size_(other.size_)
    , order_(other.order_)
  {
  }

  const DiscreteFunctionType& discreteFunction() const
  {
    return discreteFunction_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const RangeFieldType get(const size_type localDofNumber) const
  {
    assert(localDofNumber < size());
    const size_type globalDofNumber = discreteFunction_.space().map().toGlobal(entity_, localDofNumber);
    assert(globalDofNumber < discreteFunction_.vector()->size());
    return discreteFunction_.vector()->get(globalDofNumber);
  }

  int order() const
  {
    return order_;
  }

  size_type size() const
  {
    return size_;
  }

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    std::vector<RangeType> baseFunctionValues(size(), RangeType(0.0));
    localBaseFunctionSet_.evaluate(x, baseFunctionValues);
    ret = 0.0;
    for (size_type i = 0; i < size(); ++i) {
      baseFunctionValues[i] *= get(i);
      ret += baseFunctionValues[i];
    }
  }

  void jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    std::vector<JacobianRangeType> baseFunctionJacobianValues(size(), JacobianRangeType(0.0));
    localBaseFunctionSet_.jacobian(x, baseFunctionJacobianValues);
    ret = 0.0;
    for (size_type i = 0; i < size(); ++i) {
      baseFunctionJacobianValues[i] *= get(i);
      ret += baseFunctionJacobianValues[i];
    }
  }

private:
  ThisType& operator=(const ThisType&);

  const DiscreteFunctionType& discreteFunction_;
  const EntityType& entity_;
  const LocalBaseFunctionSetType localBaseFunctionSet_;
  const size_type size_;
  const int order_;
}; // end class LocalConst

template <class DiscreteFunctionImp>
class Local : public LocalConst<DiscreteFunctionImp>
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef Local<DiscreteFunctionType> ThisType;

  typedef LocalConst<DiscreteFunctionType> BaseType;

  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionType::size_type size_type;

  Local(DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : BaseType(const_cast<const DiscreteFunctionType&>(discreteFunction), entity)
    , discreteFunction_(discreteFunction)
  {
  }

  using BaseType::entity;

  using BaseType::get;

  using BaseType::size;

  using BaseType::evaluate;

  using BaseType::jacobian;

  void set(const size_type _localDofNumber, const RangeFieldType& _val)
  {
    assert(_localDofNumber < size());
    const size_type globalDofNumber = discreteFunction_.space().map().toGlobal(BaseType::entity(), _localDofNumber);
    assert(globalDofNumber < discreteFunction_.vector()->size());
    discreteFunction_.vector()->set(globalDofNumber, _val);
  }

private:
  DiscreteFunctionType& discreteFunction_;
}; // end class Local

} // namespace DiscreteFunction
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_LOCAL_HH
