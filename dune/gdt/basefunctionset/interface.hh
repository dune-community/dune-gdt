// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH
#define DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/crtp.hh>

namespace Dune {
namespace GDT {


template< class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class BaseFunctionSetInterface
  : public Stuff::LocalfunctionSetInterface< typename Traits::EntityType, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
  , public Stuff::CRTPInterface< BaseFunctionSetInterface< Traits, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >,
                                 Traits >
{
  typedef Stuff::LocalfunctionSetInterface
      < typename Traits::EntityType, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  BaseFunctionSetInterface(const EntityType& ent)
    : BaseType(ent)
  {}

  const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp(*this).backend());
    return this->as_imp(*this).backend();
  }
}; // class BaseFunctionSetInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH
