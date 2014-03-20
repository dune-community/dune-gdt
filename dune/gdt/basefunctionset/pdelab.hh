// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH
#define DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <dune/stuff/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabWrapper;


template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabWrapperTraits
{
public:
  typedef PdelabWrapper
      < PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
    derived_type;
private:
  typedef PDELab::LocalFunctionSpace< PdelabSpaceImp, PDELab::TrialSpaceTag > PdelabLFSType;
  typedef FiniteElementInterfaceSwitch< typename PdelabSpaceImp::Traits::FiniteElementType > FESwitchType;
public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;
private:
  friend class PdelabWrapper
  < PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >;
};


template< class PdelabSpaceType, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim >
class PdelabWrapper< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< PdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                       DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
                                     DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef PdelabWrapper
      < PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;
  typedef BaseFunctionSetInterface
      < PdelabWrapperTraits< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
        DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
    BaseType;
public:
  typedef PdelabWrapperTraits
      < PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > Traits;
  typedef typename Traits::BackendType    BackendType;
  typedef typename Traits::EntityType     EntityType;
private:
  typedef typename Traits::PdelabLFSType  PdelabLFSType;
  typedef typename Traits::FESwitchType   FESwitchType;

public:
  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  PdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , order_(0)
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_ = std::unique_ptr< PdelabLFSType >(lfs_ptr);
    backend_ = std::unique_ptr< BackendType >(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
    order_ = backend_->order();
  }

  PdelabWrapper(ThisType&& source)
    : BaseType(source.entity())
    , order_(std::move(source.order_))
    , lfs_(std::move(source.lfs_))
    , backend_(std::move(source.backend_))
  {}

  PdelabWrapper(const ThisType& /*other*/) = delete;

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
    return order_;
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, ret);
    DomainType tmp(0.0);
    const auto jocobian_inverse_transposed = this->entity().geometry().jacobianInverseTransposed(xx);
    for (size_t ii = 0; ii < ret.size(); ++ii) {
      jocobian_inverse_transposed.mv(ret[ii][0], tmp);
      ret[ii][0] = tmp;
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  size_t order_;
  std::unique_ptr< const PdelabLFSType > lfs_;
  std::unique_ptr< const BackendType > backend_;
}; // class PdelabWrapper


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH
