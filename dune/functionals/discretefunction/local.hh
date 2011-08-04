#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH

namespace Dune {

namespace Functionals {

namespace DiscreteFunction {

template <class DiscreteFunctionImp>
class Local
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef Local<DiscreteFunctionType> ThisType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::LocalBaseFunctionSetType LocalBaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  typedef typename DiscreteFunctionType::StorageType StorageType;

  typedef typename DiscreteFunctionType::DomainType DomainType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionType::RangeType RangeType;

  typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;

  Local(DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : discreteFunction_(discreteFunction)
    , space_(discreteFunction.space())
    , entity_(entity)
    , baseFunctionSet_(space_.localBaseFunctionSet(entity_))
  {
  }

  Local(const DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : discreteFunction_(const_cast<DiscreteFunctionType&>(discreteFunction))
    , space_(discreteFunction.space())
    , entity_(entity)
    , baseFunctionSet_(space_.localBaseFunctionSet(entity_))
  {
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  RangeFieldType& operator[](const unsigned int localDofNumber)
  {
    const unsigned int globalDofNumber = space_.mapToGlobal(entity_, localDofNumber);
    return discreteFunction_[globalDofNumber];
  }

  const RangeFieldType& operator[](const unsigned int localDofNumber) const
  {
    const unsigned int globalDofNumber = space_.mapToGlobal(entity_, localDofNumber);
    return discreteFunction_[globalDofNumber];
  }

  int order() const
  {
    return baseFunctionSet_.order();
  }

  unsigned int numDofs() const
  {
    return baseFunctionSet_.numBaseFunctions();
  }

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    std::vector<RangeType> baseFunctionValues(0.0);
    baseFunctionSet_.evaluateAll(x, baseFunctionValues);
    ret = 0.0;
    for (unsigned int i = 0; i < numDofs(); ++i) {
      baseFunctionValues[i] *= operator[](i);
      ret += baseFunctionValues[i];
    }
  }

  void jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    std::vector<JacobianRangeType> baseFunctionJacobianValues(0.0);
    baseFunctionSet_.jacobianAll(x, baseFunctionJacobianValues);
    ret = 0.0;
    for (unsigned int i = 0; i < numDofs(); ++i) {
      baseFunctionJacobianValues[i] *= operator[](i);
      ret += baseFunctionJacobianValues[i];
    }
  }

private:
  DiscreteFunctionType& discreteFunction_;
  const DiscreteFunctionSpaceType& space_;
  const EntityType& entity_;
  const LocalBaseFunctionSetType baseFunctionSet_;

}; // end class Local

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_LOCAL_HH
