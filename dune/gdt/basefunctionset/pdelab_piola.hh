// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_PDELAB_PIOLA_HH
#define DUNE_GDT_BASEFUNCTIONSET_PDELAB_PIOLA_HH

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
class PiolaTransformedPdelabWrapper
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "Untested for these dimensions!");
};

namespace internal {
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim >
class PiolaTransformedPdelabWrapperTraits
{
  static_assert(domainDim == rangeDim, "Untested!");
public:
  typedef PiolaTransformedPdelabWrapper < PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      derived_type;
private:
  typedef PDELab::LocalFunctionSpace< PdelabSpaceImp, PDELab::TrialSpaceTag > PdelabLFSType;
  typedef FiniteElementInterfaceSwitch< typename PdelabSpaceImp::Traits::FiniteElementType > FESwitchType;
public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;
private:
  friend class PiolaTransformedPdelabWrapper < PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >;
};
} // namespace internal {

template< class PdelabSpaceType, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim >
class PiolaTransformedPdelabWrapper< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< internal::PiolaTransformedPdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                                       DomainFieldImp, domainDim,
                                                                       RangeFieldImp, rangeDim >,
                                     DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef PiolaTransformedPdelabWrapper
    < PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;
  typedef BaseFunctionSetInterface< internal::PiolaTransformedPdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                                      DomainFieldImp, domainDim,
                                                                      RangeFieldImp, rangeDim >,
                                    DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
      BaseType;
public:
  typedef internal::PiolaTransformedPdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                            DomainFieldImp, domainDim,
                                            RangeFieldImp, rangeDim >
                                         Traits;
  typedef typename Traits::BackendType   BackendType;
  typedef typename Traits::EntityType    EntityType;
private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType  FESwitchType;

public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  PiolaTransformedPdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(DomainFieldType(0))
    , tmp_jacobian_transposed_(DomainFieldType(0))
    , tmp_jacobian_inverse_transposed_(DomainFieldType(0))
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_ = std::unique_ptr< PdelabLFSType >(lfs_ptr);
    backend_ = std::unique_ptr< BackendType >(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
    tmp_ranges_ = std::vector< RangeType >(backend_->size(), RangeType(0));
    tmp_jacobian_ranges_ = std::vector< JacobianRangeType >(backend_->size(), JacobianRangeType(0));
  } // PdelabWrapper(...)

  PiolaTransformedPdelabWrapper(ThisType&& source) = default;

  PiolaTransformedPdelabWrapper(const ThisType& /*other*/) = delete;

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
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(tmp_ranges_.size() >= backend_->size());
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, tmp_ranges_);
    const auto geometry = this->entity().geometry();
    tmp_jacobian_transposed_ = geometry.jacobianTransposed(xx);
    const DomainFieldType integration_element = geometry.integrationElement(xx);
    for (size_t ii = 0; ii < backend_->size(); ++ii) {
      tmp_jacobian_transposed_.mtv(tmp_ranges_[ii], ret[ii]);
      ret[ii] /= integration_element;
    }
  } // ... evaluate(...)

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, tmp_jacobian_ranges_);
    const auto geometry = this->entity().geometry();
    tmp_jacobian_transposed_ = geometry.jacobianTransposed(xx);
    tmp_jacobian_inverse_transposed_ = geometry.jacobianInverseTransposed(xx);
    const DomainFieldType integration_element = geometry.integrationElement(xx);
    for (size_t ii = 0; ii < backend_->size(); ++ii) {
      for (size_t jj = 0; jj < dimDomain; ++jj) {
        tmp_jacobian_inverse_transposed_.mv(tmp_jacobian_ranges_[ii][jj], ret[ii][jj]);
        tmp_jacobian_transposed_.mv(ret[ii][jj], tmp_jacobian_ranges_[ii][jj]);
        tmp_jacobian_ranges_[ii][jj] /= integration_element;
        ret[ii][jj] = tmp_jacobian_ranges_[ii][jj];
      }
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  mutable typename EntityType::Geometry::JacobianTransposed tmp_jacobian_transposed_;
  mutable typename EntityType::Geometry::JacobianInverseTransposed tmp_jacobian_inverse_transposed_;
  std::unique_ptr< const PdelabLFSType > lfs_;
  std::unique_ptr< const BackendType > backend_;
  mutable std::vector< RangeType > tmp_ranges_;
  mutable std::vector< JacobianRangeType > tmp_jacobian_ranges_;
}; // class PiolaTransformedPdelabWrapper
#else // HAVE_DUNE_PDELAB
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class PiolaTransformedPdelabWrapper
{
  static_assert(AlwaysFalse< PdelabSpaceImp >::value, "You are missing dune-pdelab!");
};
#endif // HAVE_DUNE_PDELAB

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_PDELAB_PIOLA_HH
