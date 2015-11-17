// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_PDELAB_VIRTUAL_HH
#define DUNE_GDT_BASEFUNCTIONSET_PDELAB_VIRTUAL_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#if HAVE_DUNE_PDELAB
# include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#endif

#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#if HAVE_DUNE_PDELAB
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class VirtualRefinedPdelabWrapper
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "Untested for arbitrary dimension!");
};

namespace internal {
// forward, to allow for specialization
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim, size_t rangeDimCols >
class VirtualRefinedPdelabWrapperTraits;

template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp >
class VirtualRefinedPdelabWrapperTraits< PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
{
public:
  typedef VirtualRefinedPdelabWrapper < PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > derived_type;
private:
  typedef PDELab::LocalFunctionSpace< PdelabSpaceImp, PDELab::TrialSpaceTag > PdelabLFSType;
  typedef FiniteElementInterfaceSwitch< typename PdelabSpaceImp::Traits::FiniteElementType > FESwitchType;
public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;
private:
  friend class VirtualRefinedPdelabWrapper< PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >;
};

} //namespace internal {


template< class PdelabSpaceType, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp >
class VirtualRefinedPdelabWrapper< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
  : public BaseFunctionSetInterface< internal::PdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                       DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >,
                                     DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
{
  typedef VirtualRefinedPdelabWrapper < PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > ThisType;
  typedef BaseFunctionSetInterface
      < internal::PdelabWrapperTraits< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >,
        DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
      BaseType;
public:
  typedef internal::VirtualRefinedPdelabWrapperTraits< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > Traits;
  typedef typename Traits::BackendType   BackendType;
  typedef typename Traits::EntityType    EntityType;
private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType  FESwitchType;

public:
  typedef typename BaseType::DomainType        DomainType;
  typedef typename BaseType::RangeType         RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  VirtualRefinedPdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(0)
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_ = std::unique_ptr< PdelabLFSType >(lfs_ptr);
    backend_ = std::unique_ptr< BackendType >(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
  } // PdelabWrapper(...)

  VirtualRefinedPdelabWrapper(ThisType&& source) = default;
  VirtualRefinedPdelabWrapper(const ThisType& /*other*/) = delete;

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
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const override final
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& /*ret*/) const override final
  {
    DUNE_THROW(NotImplemented, "");
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  std::unique_ptr< const PdelabLFSType > lfs_;
  std::unique_ptr< const BackendType > backend_;
}; // class PdelabWrapper
#endif // HAVE_DUNE_PDELAB

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_PDELAB_VIRTUAL_HH
