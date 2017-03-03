// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH

#include <vector>

#include <dune/xt/common/type_traits.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class BaseFunctionSetMapImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class DuneFemLocalfunctionsWrapper
{
  static_assert(AlwaysFalse<BaseFunctionSetMapImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class BaseFunctionSetMapImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class DuneFemLocalfunctionsWrapperTraits
{
public:
  typedef DuneFemLocalfunctionsWrapper<BaseFunctionSetMapImp,
                                       DomainFieldImp,
                                       domainDim,
                                       RangeFieldImp,
                                       rangeDim,
                                       rangeDimCols>
      derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType BackendType;
  typedef typename BackendType::EntityType EntityType;
};


} // namespace internal


template <class BaseFunctionSetMapImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class DuneFemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::DuneFemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp,
                                                                                   DomainFieldImp,
                                                                                   domainDim,
                                                                                   RangeFieldImp,
                                                                                   rangeDim,
                                                                                   1>,
                                      DomainFieldImp,
                                      domainDim,
                                      RangeFieldImp,
                                      rangeDim,
                                      1>
{
  typedef BaseFunctionSetInterface<internal::DuneFemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp,
                                                                                DomainFieldImp,
                                                                                domainDim,
                                                                                RangeFieldImp,
                                                                                rangeDim,
                                                                                1>,
                                   DomainFieldImp,
                                   domainDim,
                                   RangeFieldImp,
                                   rangeDim,
                                   1>
      BaseType;
  typedef DuneFemLocalfunctionsWrapper<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      ThisType;

public:
  typedef internal::
      DuneFemLocalfunctionsWrapperTraits<BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
          Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  DuneFemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& ent)
    : BaseType(ent)
    , baseFunctionSetMap_(baseFunctionSetMap)
    , backend_(new BackendType(baseFunctionSetMap_.find(this->entity())))
  {
  }

  DuneFemLocalfunctionsWrapper(ThisType&& source) = default;

  DuneFemLocalfunctionsWrapper(const ThisType& /*other*/) = delete;

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
}; // class DuneFemLocalfunctionsWrapper


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_FEM_LOCALFUNCTIONS_WRAPPER_HH
