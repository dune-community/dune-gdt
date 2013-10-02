// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH
#define DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 *  \brief Interface for matrix valued basis functions.
 *
 *  \note   see specialization for rangeDimCols = 1 for vector and scalar valued basis functions.
 */
template< class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
class BaseFunctionSetInterface;


/**
 *  \brief Interface for scalar and vector valued basis functions.
 */
template< class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class BaseFunctionSetInterface< Traits, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public Stuff::LocalfunctionSetInterface< typename Traits::EntityType, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
{
  typedef Stuff::LocalfunctionSetInterface
      < typename Traits::EntityType, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
    BaseType;
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  static const unsigned int                             dimRangeRows = rangeDim;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class BaseFunctionSetInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_INTERFACE_HH
