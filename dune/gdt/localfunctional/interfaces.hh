// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALFUNCTIONAL_INTERFACES_HH
#define DUNE_GDT_LOCALFUNCTIONAL_INTERFACES_HH

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class LocalFunctionalInterface
  : public Stuff::CRTPInterface< LocalFunctionalInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  template< class E, class D, size_t d, class R, size_t r, size_t rC >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
             Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBase, ret));
  }
}; // class LocalFunctionalInterface


template< class Traits >
class LocalFaceFunctionalInterface
  : public Stuff::CRTPInterface< LocalFaceFunctionalInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  template< class E, class IntersectionType, class D, size_t d, class R, size_t r, size_t rC >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
             const IntersectionType& intersection,
             Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_CR(this->as_imp().apply(testBase, intersection, ret));
  }
}; // class LocalFaceFunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_INTERFACES_HH
