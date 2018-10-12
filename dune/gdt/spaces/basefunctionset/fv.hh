// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2017 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_FV_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_FV_HH

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class FiniteVolume
{
  static_assert(Dune::AlwaysFalse<EntityImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class FiniteVolumeTraits
{
public:
  typedef FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef RangeFieldImp BackendType;
  typedef EntityImp EntityType;
};


} // namespace internal


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp>
class FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public BaseFunctionSetInterface<internal::
                                          FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>,
                                      DomainFieldImp,
                                      domainDim,
                                      RangeFieldImp,
                                      1,
                                      1>
{
  typedef FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::
                                       FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>,
                                   DomainFieldImp,
                                   domainDim,
                                   RangeFieldImp,
                                   1,
                                   1>
      BaseType;

public:
  typedef internal::FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  FiniteVolume(const EntityType& en)
    : BaseType(en)
    , backend_(1.)
  {
  }

  FiniteVolume(ThisType&& source) = default;

  FiniteVolume(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return backend_;
  }

  virtual size_t size() const override final
  {
    return 1;
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 0;
  }

  void evaluate(const DomainType& /*xx*/,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    DXT_ASSERT(ret.size() >= 0);
    ret[0] = backend_;
  }

  using BaseType::evaluate;

  void jacobian(const DomainType& /*xx*/,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    DXT_ASSERT(ret.size() >= 0);
    ret[0] *= 0.0;
  }

  using BaseType::jacobian;

private:
  const BackendType backend_;
}; // class FiniteVolume< ..., 1, 1 >


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::FiniteVolumeTraits<EntityImp,
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
  typedef FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::FiniteVolumeTraits<EntityImp,
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

public:
  typedef internal::FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  using typename BaseType::DomainType;
  using BaseType::dimRange;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  FiniteVolume(const EntityType& en)
    : BaseType(en)
    , backend_(1.)
  {
  }

  FiniteVolume(ThisType&& source) = default;

  FiniteVolume(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return backend_;
  }

  virtual size_t size() const override final
  {
    return dimRange;
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 0;
  }

  void evaluate(const DomainType& /*xx*/,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    DXT_ASSERT(ret.size() >= size());
    for (size_t ii = 0; ii < size(); ++ii) {
      ret[ii] *= 0.0;
      ret[ii][ii] = backend_;
    }
  } // ... evaluate(...)

  using BaseType::evaluate;

  void jacobian(const DomainType& /*xx*/,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    DXT_ASSERT(ret.size() >= size());
    for (size_t ii = 0; ii < size(); ++ii)
      ret[ii] *= 0.0;
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  const BackendType backend_;
}; // class FiniteVolume< ..., rangeDim, 1 >


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_FV_HH
