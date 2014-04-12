// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH

#include <type_traits>

#include <vector>

#include <dune/common/typetraits.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapper
{
  static_assert(AlwaysFalse< BaseFunctionSetMapImp >::value, "Untested for these dimensions!");
};


template< class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols >
class FemLocalfunctionsWrapperTraits
{
public:
  typedef FemLocalfunctionsWrapper< BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >  derived_type;
  typedef typename BaseFunctionSetMapImp::BaseFunctionSetType BackendType;
  typedef typename BackendType::EntityType                    EntityType;
};


template< class BaseFunctionSetMapImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class FemLocalfunctionsWrapper< BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface
        < FemLocalfunctionsWrapperTraits
            < BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
              DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef BaseFunctionSetInterface
      < FemLocalfunctionsWrapperTraits< BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
        DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
    BaseType;
  typedef FemLocalfunctionsWrapper< BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
    ThisType;
public:
  typedef FemLocalfunctionsWrapperTraits< BaseFunctionSetMapImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > Traits;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  FemLocalfunctionsWrapper(const BaseFunctionSetMapImp& baseFunctionSetMap, const EntityType& ent)
    : BaseType(ent)
    , baseFunctionSetMap_(baseFunctionSetMap)
    , backend_(new BackendType(baseFunctionSetMap_.find(this->entity())))
  {}

  FemLocalfunctionsWrapper(ThisType&& source)
    : BaseType(source.entity())
    , baseFunctionSetMap_(source.baseFunctionSetMap_)
    , backend_(std::move(source.backend_))
  {}

  FemLocalfunctionsWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const DS_OVERRIDE
  {
    return backend_->size();
  }

  virtual size_t order() const DS_OVERRIDE
  {
    return baseFunctionSetMap_.getOrder(this->entity());
  }

  virtual void evaluate(const DomainType& x, std::vector< RangeType >& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= size());
    backend_->evaluateAll(x, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& x, std::vector< JacobianRangeType >& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= size());
    backend_->jacobianAll(x, this->entity().geometry().jacobianInverseTransposed(x), ret);
  }

  using BaseType::jacobian;

private:
  const BaseFunctionSetMapImp& baseFunctionSetMap_;
  std::unique_ptr< const BackendType > backend_;
}; // class FemLocalfunctionsWrapper


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_LOCALFUNCTIONS_HH
