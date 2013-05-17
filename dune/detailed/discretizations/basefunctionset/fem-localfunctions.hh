#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim,
          int rangeDimCols = 1>
class FemLocalfunctionsWrapper;


template <class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim,
          int rangeDimCols = 1>
class FemLocalfunctionsWrapperTraits
{
public:
  typedef FemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                   rangeDimCols> derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType BackendType;
  typedef typename BackendType::EntityType EntityType;
  //  typedef typename BackendType::DomainFieldType   DomainFieldType;
  //  static const unsigned int                       dimDomain = BackendType::dimDomain;
  //  typedef typename BackendType::DomainType        DomainType;
  //  typedef typename BackendType::RangeFieldType    RangeFieldType;
  //  static const unsigned int                       dimRange = BackendType::dimRange;
  //  typedef typename BackendType::RangeType         RangeType;
  //  typedef typename BackendType::JacobianRangeType JacobianRangeType;
};


template <class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<FemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp, domainDim,
                                                                     RangeFieldImp, rangeDim, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef FemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  FemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& en)
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
    assert(ret.size() >= size());
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
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
