// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TYPE_TRAITS_HH
#define DUNE_GDT_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/xt/common/type_traits.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

namespace Dune {
namespace GDT {


// ///////////////////////////////////////////
// // all of local/integrands/interfaces.hh //
// ///////////////////////////////////////////

namespace internal {


template <class Tt>
struct is_unary_volume_integrand_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


template <class Tt>
struct is_binary_volume_integrand_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


template <class Tt>
struct is_unary_face_integrand_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


template <class Tt>
struct is_binary_face_integrand_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


template <class Tt>
struct is_quaternary_face_integrand_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


} // namespace internal


template <class T, bool candidate = internal::is_unary_volume_integrand_helper<T>::is_candidate>
struct is_unary_volume_integrand : public std::is_base_of<LocalVolumeIntegrandInterface<typename T::Traits, 1>, T>
{
};

template <class T>
struct is_unary_volume_integrand<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_binary_volume_integrand_helper<T>::is_candidate>
struct is_binary_volume_integrand : public std::is_base_of<LocalVolumeIntegrandInterface<typename T::Traits, 2>, T>
{
};

template <class T>
struct is_binary_volume_integrand<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_unary_face_integrand_helper<T>::is_candidate>
struct is_unary_face_integrand : public std::is_base_of<LocalFaceIntegrandInterface<typename T::Traits, 1>, T>
{
};

template <class T>
struct is_unary_face_integrand<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_binary_face_integrand_helper<T>::is_candidate>
struct is_binary_face_integrand : public std::is_base_of<LocalFaceIntegrandInterface<typename T::Traits, 2>, T>
{
};

template <class T>
struct is_binary_face_integrand<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_quaternary_face_integrand_helper<T>::is_candidate>
struct is_quaternary_face_integrand : public std::is_base_of<LocalFaceIntegrandInterface<typename T::Traits, 4>, T>
{
};

template <class T>
struct is_quaternary_face_integrand<T, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TYPE_TRAITS_HH
