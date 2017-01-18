// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/spaces/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class LocalVolumeFunctionalInterface : public XT::CRTPInterface<LocalVolumeFunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  void apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_basis,
             Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(test_basis, ret));
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  Dune::DynamicVector<R> apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_basis) const
  {
    Dune::DynamicVector<R> ret(test_basis.size(), 0.);
    apply(test_basis, ret);
    return ret;
  }
}; // class LocalFunctionalInterface


template <class Traits>
class LocalFaceFunctionalInterface : public XT::CRTPInterface<LocalFaceFunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  void apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_basis,
             const IntersectionType& intersection,
             Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(test_basis, intersection, ret));
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  Dune::DynamicVector<R> apply(const XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>& test_basis,
                               const IntersectionType& intersection) const
  {
    Dune::DynamicVector<R> ret(test_basis.size(), 0.);
    apply(test_basis, intersection, ret);
    return ret;
  }
}; // class LocalFaceFunctionalInterface


namespace internal {


template <class Tt>
struct is_local_volume_functional_helper
{
  DXTC_has_typedef_initialize_once(Traits);

  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


} // namespace internal


template <class T, bool candidate = internal::is_local_volume_functional_helper<T>::is_candidate>
struct is_local_volume_functional : public std::is_base_of<LocalVolumeFunctionalInterface<typename T::Traits>, T>
{
};

template <class T>
struct is_local_volume_functional<T, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
