// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_HH

#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/space/basisfunctionset/default.hh>
# include <dune/fem/space/common/discretefunctionspace.hh>
#endif // HAVE_DUNE_FEM

#include <dune/stuff/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template< class ShapeFunctionSetImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemWrapper
{
  static_assert(Dune::AlwaysFalse< ShapeFunctionSetImp >::value, "Untested for these dimensions!");
};


template< class ShapeFunctionSetImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim, int rangeDimCols >
class FemWrapperTraits
{
public:
  typedef FemWrapper
      < ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
    derived_type;
  typedef typename Dune::Fem::DefaultBasisFunctionSet< EntityImp, ShapeFunctionSetImp > BackendType;
  typedef EntityImp EntityType;
};


template< class ShapeFunctionSetImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim >
class FemWrapper< ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< FemWrapperTraits< ShapeFunctionSetImp, EntityImp,
                                                       DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
                                     DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef FemWrapper
      < ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;
  typedef BaseFunctionSetInterface
      < FemWrapperTraits< ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >,
        DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
    BaseType;
public:
  typedef FemWrapperTraits
      < ShapeFunctionSetImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > Traits;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  template< class S >
  FemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface< S >& femSpace, const EntityType& ent)
    : BaseType(ent)
    , backend_(new BackendType(femSpace.basisFunctionSet(this->entity())))
  {}

  FemWrapper(ThisType&& source)
    : BaseType(source.entity())
    , backend_(std::move(source.backend_))
  {}

  FemWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const DS_OVERRIDE DS_FINAL
  {
    return backend_->size();
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    assert(backend_->order() >= 0);
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= size());
    backend_->evaluateAll(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= size());
    backend_->jacobianAll(xx, ret);
  }

  using BaseType::jacobian;

private:
  std::unique_ptr< const BackendType > backend_;
}; // class FemWrapper


#else // HAVE_DUNE_FEM


template< class ShapeFunctionSetImp, class EntityImp,
          class DomainFieldImp, int domainDim,
          class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemWrapper
{
  static_assert(Dune::AlwaysFalse< FemBaseFunctionSetTraits >::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_HH
