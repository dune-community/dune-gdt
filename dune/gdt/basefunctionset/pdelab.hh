// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH
#define DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/common/typetraits.hh>

#ifdef HAVE_DUNE_PDELAB
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#endif

#include <dune/stuff/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#ifdef HAVE_DUNE_PDELAB

// forward, to be used in the traits and to allow for specialization
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim,
          int rangeDimCols = 1>
class PdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "Untested for arbitrary dimension!");
};


// forward, to allow for specialization
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim,
          int rangeDimCols = 1>
class PdelabWrapperTraits;


//! Specialization for dimRange = 1, dimRangeRows = 1
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp>
class PdelabWrapperTraits<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
public:
  typedef PdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> derived_type;

private:
  typedef PDELab::LocalFunctionSpace<PdelabSpaceImp, PDELab::TrialSpaceTag> PdelabLFSType;
  typedef FiniteElementInterfaceSwitch<typename PdelabSpaceImp::Traits::FiniteElementType> FESwitchType;

public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;

private:
  friend class PdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>;
};


//! Specialization for dimRange = 1, dimRangeRows = 1
template <class PdelabSpaceType, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp>
class PdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public BaseFunctionSetInterface<PdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim,
                                                          RangeFieldImp, 1, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
  typedef PdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;
  typedef BaseFunctionSetInterface<PdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim,
                                                       RangeFieldImp, 1, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, 1, 1> BaseType;

public:
  typedef PdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType FESwitchType;

public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = 1;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  PdelabWrapper(const PdelabSpaceType& space, const EntityType& ent, const size_t oo)
    : BaseType(ent)
    , order_(oo)
    , tmp_domain_(RangeFieldType(0))
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_     = std::unique_ptr<PdelabLFSType>(lfs_ptr);
    backend_ = std::unique_ptr<BackendType>(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
    order_   = backend_->order();
  } // PdelabWrapper(...)

  PdelabWrapper(ThisType&& source)
    : BaseType(source.entity())
    , order_(std::move(source.order_))
    , tmp_domain_(RangeFieldType(0))
    , lfs_(std::move(source.lfs_))
    , backend_(std::move(source.backend_))
  {
  }

  PdelabWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const DS_FINAL
  {
    return backend_->size();
  }

  virtual size_t order() const DS_FINAL
  {
    return order_;
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const DS_FINAL
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, ret);
  } // ... evaluate(...)

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const DS_FINAL
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, ret);
    const auto jacobian_inverse_transposed = this->entity().geometry().jacobianInverseTransposed(xx);
    for (size_t ii = 0; ii < ret.size(); ++ii) {
      jacobian_inverse_transposed.mv(ret[ii][0], tmp_domain_);
      ret[ii][0] = tmp_domain_;
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  size_t order_;
  mutable DomainType tmp_domain_;
  std::unique_ptr<const PdelabLFSType> lfs_;
  std::unique_ptr<const BackendType> backend_;
}; // class PdelabWrapper


#else // HAVE_DUNE_PDELAB


template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim,
          int rangeDimCols = 1>
class PdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_PDELAB_HH
