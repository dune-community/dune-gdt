// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TYPE_TRAITS_HH
#define DUNE_GDT_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace GDT {

// forwards
// from #include <dune/gdt/spaces/interface.hh>
enum class ChooseSpaceBackend;

enum class ChoosePattern;

template <class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols>
class SpaceInterface;


// from #include <dune/gdt/local/integrands/interfaces.hh>
template <class Traits, size_t numArguments>
class LocalVolumeIntegrandInterface;

template <class Traits, size_t numArguments>
class LocalFaceIntegrandInterface;

// from #include <dune/gdt/operators/interfaces.hh>
template <class Traits>
class OperatorInterface;

// from #include <dune/gdt/operators/base.hh>
template <class M, class RS, class GV, class SS, class F, ChoosePattern pt>
class MatrixOperatorBase;

template <class GridViewImp, class RangeImp, class SourceImp, class FieldImp>
class LocalizableProductBase;

template <class GridViewImp, class SourceImp, class RangeImp>
class LocalizableOperatorBase;


namespace internal {


// helper structs
// from #include <dune/gdt/local/integrands/interfaces.hh>
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


// from #include <dune/gdt/operators/interfaces.hh>
template <class Tt>
struct is_operator_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  static const bool is_candidate = DXTC_has_typedef(Traits)<Tt>::value;
};


// from #include <dune/gdt/operators/base.hh>
template <class Tt>
struct is_localizable_product_helper
{
  DXTC_has_typedef_initialize_once(GridViewType);
  DXTC_has_typedef_initialize_once(RangeType);
  DXTC_has_typedef_initialize_once(SourceType);
  DXTC_has_typedef_initialize_once(FieldType);
  static const bool is_candidate = DXTC_has_typedef(GridViewType)<Tt>::value && DXTC_has_typedef(RangeType)<Tt>::value
                                   && DXTC_has_typedef(SourceType)<Tt>::value && DXTC_has_typedef(FieldType)<Tt>::value;
};


template <class Tt>
struct is_matrix_operator_helper
{
  DXTC_has_typedef_initialize_once(MatrixType);
  DXTC_has_typedef_initialize_once(RangeSpaceType);
  DXTC_has_typedef_initialize_once(GridViewType);
  DXTC_has_typedef_initialize_once(SourceSpaceType);
  DXTC_has_typedef_initialize_once(FieldType);
  DXTC_has_static_member_initialize_once(pattern_type);
  static const bool is_candidate =
      DXTC_has_typedef(MatrixType)<Tt>::value && DXTC_has_typedef(RangeSpaceType)<Tt>::value
      && DXTC_has_typedef(GridViewType)<Tt>::value && DXTC_has_typedef(SourceSpaceType)<Tt>::value
      && DXTC_has_typedef(FieldType)<Tt>::value && DXTC_has_static_member(pattern_type)<Tt>::value;
};


} // namespace internal


// actual structs
// from #include <dune/gdt/local/integrands/interfaces.hh>
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


// from #include <dune/gdt/operators/interfaces.hh>
template <class T, bool candidate = internal::is_operator_helper<T>::is_candidate>
struct is_operator : public std::is_base_of<OperatorInterface<typename T::Traits>, T>
{
};

template <class T>
struct is_operator<T, false> : public std::false_type
{
};


// from #include <dune/gdt/operators/base.hh>
template <class T, bool candidate = internal::is_localizable_product_helper<T>::is_candidate>
struct is_localizable_product : public std::is_base_of<LocalizableProductBase<typename T::GridViewType,
                                                                              typename T::RangeType,
                                                                              typename T::SourceType,
                                                                              typename T::FieldType>,
                                                       T>
{
};

template <class T>
struct is_localizable_product<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_matrix_operator_helper<T>::is_candidate>
struct is_matrix_operator : public std::is_base_of<MatrixOperatorBase<typename T::MatrixType,
                                                                      typename T::RangeSpaceType,
                                                                      typename T::GridViewType,
                                                                      typename T::SourceSpaceType,
                                                                      typename T::FieldType,
                                                                      T::pattern_type>,
                                                   T>
{
};

template <class T>
struct is_matrix_operator<T, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TYPE_TRAITS_HH
