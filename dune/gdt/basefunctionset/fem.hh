// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/stuff/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp,
          int rangeDim, int rangeDimCols = 1>
class FemWrapper;


template <class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp,
          int rangeDim, int rangeDimCols = 1>
class FemWrapperTraits
{
public:
  typedef FemWrapper<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                     rangeDimCols> derived_type;
  typedef typename Dune::Fem::BaseFunctionSetInterface<FemBaseFunctionSetTraits>::BaseFunctionSetType BackendType;
  typedef EntityImp EntityType;
};


template <class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp,
          int rangeDim>
class FemWrapper<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<FemWrapperTraits<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim,
                                                       RangeFieldImp, rangeDim, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef FemWrapper<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      ThisType;
  typedef BaseFunctionSetInterface<FemWrapperTraits<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim,
                                                    RangeFieldImp, rangeDim, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

public:
  typedef FemWrapperTraits<FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  template <class S>
  FemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface<S>& femSpace, const EntityType& ent)
    : BaseType(ent)
    , order_(femSpace.order())
    , backend_(new BackendType(femSpace.baseFunctionSet(this->entity())))
  {
  }

  FemWrapper(ThisType&& source)
    : BaseType(source.entity())
    , order_(std::move(source.order_))
    , backend_(std::move(source.backend_))
  {
  }

  FemWrapper(const ThisType& /*other*/) = delete;

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

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateAll(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const DS_OVERRIDE
  {
    assert(ret.size() >= backend_->size());
    backend_->jacobianAll(xx, this->entity().geometry().jacobianInverseTransposed(xx), ret);
  }

  using BaseType::jacobian;

private:
  const size_t order_;
  std::unique_ptr<const BackendType> backend_;
}; // class FemWrapper


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_HH
