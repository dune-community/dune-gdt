// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SAPCES_BASEFUNCTIONSET_DUNE_FEM_WRAPPER_HH
#define DUNE_GDT_SAPCES_BASEFUNCTIONSET_DUNE_FEM_WRAPPER_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#endif

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class BasisFunctionSetImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class DuneFemWrapper
{
  static_assert(Dune::AlwaysFalse<BasisFunctionSetImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class BasisFunctionSetImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class DuneFemWrapperTraits
{
public:
  typedef DuneFemWrapper<BasisFunctionSetImp,
                         EntityImp,
                         DomainFieldImp,
                         domainDim,
                         RangeFieldImp,
                         rangeDim,
                         rangeDimCols>
      derived_type;
  typedef typename Dune::Fem::DefaultBasisFunctionSet<EntityImp, typename BasisFunctionSetImp::ShapeFunctionSetType>
      BackendType;
  typedef EntityImp EntityType;
};


} // namespace internal


template <class BasisFunctionSetImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim>
class DuneFemWrapper<BasisFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::DuneFemWrapperTraits<BasisFunctionSetImp,
                                                                     EntityImp,
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
  typedef DuneFemWrapper<BasisFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      ThisType;
  typedef BaseFunctionSetInterface<internal::DuneFemWrapperTraits<BasisFunctionSetImp,
                                                                  EntityImp,
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
  typedef internal::
      DuneFemWrapperTraits<BasisFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
          Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  template <class S>
  DuneFemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface<S>& femSpace, const EntityType& ent)
    : BaseType(ent)
    , backend_(new BackendType(ent, femSpace.basisFunctionSet(ent).shapeFunctionSet()))
  {
  }

  DuneFemWrapper(ThisType&& source) = default;

  DuneFemWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size();
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    assert(backend_->order() >= 0);
    return backend_->order();
  }

  void evaluate(const DomainType& xx,
                std::vector<RangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = XT::Common::Parameter()) const override final
  {
    assert(ret.size() >= size());
    backend_->evaluateAll(xx, ret);
  }

  using BaseType::evaluate;

  void jacobian(const DomainType& xx,
                std::vector<JacobianRangeType>& ret,
                const XT::Common::Parameter& /*mu*/ = XT::Common::Parameter()) const override final
  {
    assert(ret.size() >= size());
    backend_->jacobianAll(xx, ret);
  }

  using BaseType::jacobian;

private:
  std::unique_ptr<const BackendType> backend_;
}; // class DuneFemWrapper


#else // HAVE_DUNE_FEM


template <class BasisFunctionSetImp,
          class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class DuneFemWrapper
{
  static_assert(Dune::AlwaysFalse<BasisFunctionSetImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SAPCES_BASEFUNCTIONSET_DUNE_FEM_WRAPPER_HH
