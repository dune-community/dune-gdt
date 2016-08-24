// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/spaces/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class LocalVolumeFunctionalInterface : public XT::CRTPInterface<LocalVolumeFunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  void apply(const Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase, Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBase, ret));
  }
}; // class LocalFunctionalInterface


template <class Traits>
class LocalFaceFunctionalInterface : public XT::CRTPInterface<LocalFaceFunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;

  template <class E, class IntersectionType, class D, size_t d, class R, size_t r, size_t rC>
  void apply(const Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase, const IntersectionType& intersection,
             Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_CR(this->as_imp().apply(testBase, intersection, ret));
  }
}; // class LocalFaceFunctionalInterface


namespace internal {


template <class Tt>
struct is_local_volume_functional_helper
{
  DXTC_has_typedef_initialize_once(Traits)

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
