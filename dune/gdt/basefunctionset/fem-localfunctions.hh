// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <vector>

#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class BaseFunctionSetMapImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
          size_t rangeDimCols = 1>
class FemLocalfunctionsWrapper
{
  static_assert(AlwaysFalse<BaseFunctionSetMapImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class BaseFunctionSetMapImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
          size_t rangeDimCols>
class FemLocalfunctionsWrapperTraits
{
public:
  typedef FemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                   rangeDimCols> derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType BackendType;
  typedef typename BackendType::EntityType EntityType;
};


} // namespace internal


template <class BaseFunctionSetMapImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class FemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::FemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp,
                                                                               domainDim, RangeFieldImp, rangeDim, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef BaseFunctionSetInterface<internal::FemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp,
                                                                            domainDim, RangeFieldImp, rangeDim, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef FemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      ThisType;

public:
  typedef internal::FemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                   rangeDim, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  FemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& ent)
    : BaseType(ent)
    , baseFunctionSetMap_(baseFunctionSetMap)
    , backend_(new BackendType(baseFunctionSetMap_.find(this->entity())))
  {
  }

  FemLocalfunctionsWrapper(ThisType&& source) = default;

  FemLocalfunctionsWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override
  {
    return backend_->size();
  }

  virtual size_t order() const override
  {
    return baseFunctionSetMap_.getOrder(this->entity());
  }

  virtual void evaluate(const DomainType& x, std::vector<RangeType>& ret) const override
  {
    assert(ret.size() >= size());
    backend_->evaluateAll(x, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const override
  {
    assert(ret.size() >= size());
    backend_->jacobianAll(x, this->entity().geometry().jacobianInverseTransposed(x), ret);
  }

  using BaseType::jacobian;

private:
  const BaseFunctionSetMapImp& baseFunctionSetMap_;
  std::unique_ptr<const BackendType> backend_;
}; // class FemLocalfunctionsWrapper


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
