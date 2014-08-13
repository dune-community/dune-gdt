// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH
#define DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH

#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/fmatrix.hh>
# include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/vector.hh>
#include <dune/stuff/common/matrix.hh>

#include "../../basefunctionset/interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FiniteVolume
{
  static_assert(Dune::AlwaysFalse< EntityImp >::value, "Untested for these dimensions!");
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols >
class FiniteVolumeTraits
{
public:
  typedef FiniteVolume< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef RangeFieldImp BackendType;
  typedef EntityImp EntityType;
};


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp >
class FiniteVolume< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
  : public BaseFunctionSetInterface< FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
                                   , DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
{
  typedef FiniteVolume< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > ThisType;
  typedef BaseFunctionSetInterface< FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
                                                      , DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > BaseType;
public:
  typedef FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > Traits;
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

  FiniteVolume(const EntityType& en)
    : BaseType(en)
    , backend_(1)
  {}

  FiniteVolume(ThisType&& source)
    : BaseType(source.entity())
    , backend_(1)
  {}

  FiniteVolume(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return backend_;
  }

  virtual size_t size() const DS_OVERRIDE DS_FINAL
  {
    return 1;
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*xx*/, std::vector< RangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= 0);
    ret[0] = 1.0;
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= 0);
    ret[0] *= 0.0;
  }

  using BaseType::jacobian;

private:
  const BackendType backend_;
}; // class FiniteVolume< ..., 1, 1 >


template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class FiniteVolume< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
                                   , DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef FiniteVolume< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;
  typedef BaseFunctionSetInterface< FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
                                                      , DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > BaseType;
public:
  typedef FiniteVolumeTraits< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > Traits;
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

  FiniteVolume(const EntityType& en)
    : BaseType(en)
    , backend_(1)
  {}

  FiniteVolume(ThisType&& source)
    : BaseType(source.entity())
    , backend_(1)
  {}

  FiniteVolume(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return backend_;
  }

  virtual size_t size() const DS_OVERRIDE DS_FINAL
  {
    return dimRange;
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*xx*/, std::vector< RangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= dimRange);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      ret[ii] *= 0.0;
      ret[ii][ii] = RangeFieldType(1);
    }
  } // ... evaluate(...)

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= dimRange);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] *= 0.0;
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  const BackendType backend_;
}; // class FiniteVolume< ..., rangeDim, 1 >


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH
