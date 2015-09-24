// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#warning This header is deprecated, include <dune/gdt/localfunctional/interfaces.hh> instead (24.09.2015)!

#ifndef DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH
#define DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH

#include <dune/common/deprecated.hh>

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/basefunctionset/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace LocalFunctional {


template< class Traits >
class
  DUNE_DEPRECATED_MSG("Use LocalVolumeFunctionalInterface instead (24.09.2015)!")
      Codim0Interface
  : public Stuff::CRTPInterface< Codim0Interface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_CRTP(this->as_imp().numTmpObjectsRequired());
    return this->as_imp().numTmpObjectsRequired();
  }

  /**
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the test BaseFunctionSetInterface
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange
   *  \tparam rC  dimRangeCols
   */
  template< class T, class D, size_t d, class R, size_t r, size_t rC >
  void apply(const BaseFunctionSetInterface< T, D, d, R, r, rC >& testBase,
             Dune::DynamicVector< R >& ret,
             std::vector< Dune::DynamicVector< R > >& tmpLocalVectors) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBase, ret, tmpLocalVectors));
  }
}; // class Codim0Interface


template< class Traits >
class
  DUNE_DEPRECATED_MSG("Use LocalFaceFunctionalInterface instead (24.09.2015)!")
      Codim1Interface
  : public Stuff::CRTPInterface< Codim1Interface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_CRTP(this->as_imp().numTmpObjectsRequired());
    return this->as_imp().numTmpObjectsRequired();
  }

  /**
   *  \brief Applies the local functional.
   *  \tparam T                 Traits of the test BaseFunctionSetInterface implementation
   *  \tparam IntersectionType  A model of Dune::Intersection< ... >
   *  \tparam D                 DomainFieldType
   *  \tparam d                 dimDomain
   *  \tparam R                 RangeFieldType
   *  \tparam r                 dimRange of the of the testBase
   *  \tparam rC                dimRangeCols of the testBase
   */
  template< class T, class IntersectionType, class D, size_t d, class R, size_t r, size_t rC >
  void apply(const BaseFunctionSetInterface< T, D, d, R, r, rC >& testBase,
             const IntersectionType& intersection,
             Dune::DynamicVector< R >& ret,
             std::vector< Dune::DynamicVector< R > >& tmpLocalVectors) const
  {
    CHECK_AND_CALL_CR(this->as_imp().apply(testBase, intersection, ret, tmpLocalVectors));
  }
}; // class Codim1Interface


} // namespace LocalFunctional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH
