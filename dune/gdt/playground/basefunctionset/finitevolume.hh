// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH
#define DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH

#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/stuff/common/memory.hh>

#include "../../basefunctionset/interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class FiniteVolume
{
  static_assert(Dune::AlwaysFalse<EntityImp>::value, "Untested for these dimensions!");
};


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols>
class FiniteVolumeTraits
{
public:
  typedef FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef RangeFieldImp BackendType;
  typedef EntityImp EntityType;
};


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp>
class FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public BaseFunctionSetInterface<FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
  typedef FiniteVolume<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;
  typedef BaseFunctionSetInterface<FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, 1, 1> BaseType;

public:
  typedef FiniteVolumeTraits<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = 1;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  FiniteVolume(const EntityType& en)
    : BaseType(en)
  {
  }

  FiniteVolume(ThisType&& source)
    : BaseType(source.entity())
  {
  }

  FiniteVolume(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return 1.0;
  }

  virtual size_t size() const DS_OVERRIDE DS_FINAL
  {
    return 1;
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*xx*/, std::vector<RangeType>& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= 0);
    ret[0] = 1.0;
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& /*xx*/, std::vector<JacobianRangeType>& ret) const DS_OVERRIDE DS_FINAL
  {
    assert(ret.size() >= 0);
    ret[0] *= 0.0;
  }

  using BaseType::jacobian;
}; // class FiniteVolume


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FINITEVOLUME_HH
