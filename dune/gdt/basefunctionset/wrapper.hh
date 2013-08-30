// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_BASEFUNCTIONSET_WRAPPER_HH
#define DUNE_GDT_BASEFUNCTIONSET_WRAPPER_HH

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


template< class LocalFunctionImp, class D, int d, class R, int rR, int rC = 1 >
class LocalFunctionWrapper;


template< class LocalFunctionImp, class D, int d, class R, int rR, int rC = 1 >
class LocalFunctionWrapperTraits
{
  static_assert(std::is_base_of< Stuff::LocalFunctionInterface< typename LocalFunctionImp::Traits, D, d, R, rR, rC >,
                                                                LocalFunctionImp >::value,
                "LocalFunctionImp has to be derived from Stuff::LocalFunctionInterface< ..., D, d, R, rR, rC >");
public:
  typedef LocalFunctionWrapper< LocalFunctionImp, D, d, R, rR, rC > derived_type;
  typedef LocalFunctionImp BackendType;
  typedef typename LocalFunctionImp::EntityType EntityType;
};


template< class LocalFunctionImp, class D, int d, class R, int r >
class LocalFunctionWrapper< LocalFunctionImp, D, d, R, r, 1 >
  : public BaseFunctionSetInterface< LocalFunctionWrapperTraits< LocalFunctionImp, D, d, R, r, 1 >, D, d, R, r, 1 >
{
public:
  typedef LocalFunctionWrapperTraits< LocalFunctionImp, D, d, R, r, 1 > Traits;

  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  typedef D                                               DomainFieldType;
  static const unsigned int                               dimDomain = d;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef R                                             RangeFieldType;
  static const unsigned int                             dimRange = r;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  LocalFunctionWrapper(const BackendType& backend)
    : backend_(backend)
  {}

  const EntityType& entity() const
  {
    return backend_.entity();
  }

  size_t size() const
  {
    return 1;
  }

  size_t order() const
  {
    return backend_.order();
  }

  void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const
  {
    backend_.evaluate(xx, ret[0]);
  }

  void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const
  {
    backend_.jacobian(xx, ret[0]);
  }

private:
  const BackendType& backend_;
};


} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_WRAPPER_HH
