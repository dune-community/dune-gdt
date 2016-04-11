// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#endif

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class ShapeFunctionSetImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class FemWrapper
{
  static_assert(Dune::AlwaysFalse<ShapeFunctionSetImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class ShapeFunctionSetImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols>
class FemWrapperTraits
{
public:
  typedef FemWrapper<ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      derived_type;
  typedef typename Dune::Fem::DefaultBasisFunctionSet<EntityImp, ShapeFunctionSetImp> BackendType;
  typedef EntityImp EntityType;
};


} // namespace internal


template <class ShapeFunctionSetImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim>
class FemWrapper<ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::FemWrapperTraits<ShapeFunctionSetImp, EntityImp, DomainFieldImp,
                                                                 domainDim, RangeFieldImp, rangeDim, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef FemWrapper<ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::FemWrapperTraits<ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim,
                                                              RangeFieldImp, rangeDim, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

public:
  typedef internal::FemWrapperTraits<ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                     1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  template <class S>
  FemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface<S>& femSpace, const EntityType& ent)
    : BaseType(ent)
    , backend_(new BackendType(femSpace.basisFunctionSet(this->entity())))
  {
  }

  FemWrapper(ThisType&& source) = default;

  FemWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size();
  }

  virtual size_t order() const override final
  {
    assert(backend_->order() >= 0);
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const override final
  {
    assert(ret.size() >= size());
    backend_->evaluateAll(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const override final
  {
    assert(ret.size() >= size());
    backend_->jacobianAll(xx, ret);
  }

  using BaseType::jacobian;

private:
  std::unique_ptr<const BackendType> backend_;
}; // class FemWrapper


#else // HAVE_DUNE_FEM


template <class ShapeFunctionSetImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class FemWrapper
{
  static_assert(Dune::AlwaysFalse<ShapeFunctionSetImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_HH
