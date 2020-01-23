// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   René Fritze     (2016, 2018)
//   Tobias Leibner  (2016, 2018)

#ifndef DUNE_GDT_TYPE_TRAITS_HH
#define DUNE_GDT_TYPE_TRAITS_HH

#include <ostream>
#include <type_traits>

#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace GDT {


enum class SpaceType
{
  continuous_flattop,
  continuous_lagrange,
  discontinuous_lagrange,
  finite_volume,
  raviart_thomas
};


inline std::ostream& operator<<(std::ostream& out, const SpaceType& space_type)
{
  if (space_type == SpaceType::continuous_lagrange)
    out << "continuous_lagrange";
  else if (space_type == SpaceType::continuous_flattop)
    out << "continuous_flattop";
  else if (space_type == SpaceType::discontinuous_lagrange)
    out << "discontinuous_lagrange";
  else if (space_type == SpaceType::finite_volume)
    out << "finite_volume";
  else if (space_type == SpaceType::raviart_thomas)
    out << "finite_volume";
  else
    out << "SpaceType(" << int(space_type) << ")";
  return out;
}


enum class Stencil
{
  element,
  intersection,
  element_and_intersection
};


inline std::ostream& operator<<(std::ostream& out, const Stencil& stencil)
{
  if (stencil == Stencil::element)
    out << "element";
  else if (stencil == Stencil::intersection)
    out << "intersection";
  else if (stencil == Stencil::element_and_intersection)
    out << "element_and_intersection";
  else
    out << "Stencil(" << int(stencil) << ")";
  return out;
}


// forwards
// from #include <dune/gdt/local/finite-elements/interfaces.hh>
template <class D, size_t d, class R, size_t r, size_t rC>
class LocalFiniteElementInterface;

template <class D, size_t d, class R, size_t r, size_t rC>
class LocalFiniteElementFamilyInterface;

// from #include <dune/gdt/spaces/interface.hh>
template <class GV, size_t r, size_t rC, class R>
class SpaceInterface;

// from #include <dune/gdt/test/momentmodels/basisfunctions.hh>
enum class EntropyType;

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          size_t dimFlux,
          EntropyType entropy>
class HatFunctionMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols, EntropyType entropy>
class LegendreMomentBasis;

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          size_t dimFlux,
          size_t order,
          EntropyType entropy>
class PartialMomentBasis;

template <class DomainFieldType,
          class RangeFieldType,
          size_t order,
          size_t fluxDim,
          bool only_positive,
          EntropyType entropy>
class SphericalHarmonicsMomentBasis;

template <class DomainFieldType,
          class RangeFieldType,
          size_t order,
          size_t fluxDim,
          bool only_even,
          EntropyType entropy>
class RealSphericalHarmonicsMomentBasis;


namespace internal {


// helper structs
// from #include <dune/gdt/spaces/interface.hh>
template <class S>
struct is_space_helper
{
  DXTC_has_typedef_initialize_once(GV);
  DXTC_has_static_member_initialize_once(r);
  DXTC_has_static_member_initialize_once(rC);
  DXTC_has_typedef_initialize_once(R);

  static const constexpr bool is_candidate = DXTC_has_typedef(GV)<S>::value && DXTC_has_static_member(r)<S>::value
                                             && DXTC_has_static_member(rC)<S>::value && DXTC_has_typedef(R)<S>::value;
};


} // namespace internal


// actual structs
// from #include <dune/gdt/local/finite-elements/interfaces.hh>
template <class T>
struct is_local_finite_element : public std::false_type
{};

template <class D, size_t d, class R, size_t r, size_t rC>
struct is_local_finite_element<LocalFiniteElementInterface<D, d, R, r, rC>> : public std::true_type
{};


template <class T>
struct is_local_finite_element_family : public std::false_type
{};

template <class D, size_t d, class R, size_t r, size_t rC>
struct is_local_finite_element_family<LocalFiniteElementFamilyInterface<D, d, R, r, rC>> : public std::true_type
{};


// from #include <dune/gdt/spaces/interface.hh>
template <class S, bool is_candidate = internal::is_space_helper<S>::is_candidate>
struct is_space : public std::false_type
{};

template <class S>
struct is_space<S, true> : public std::is_base_of<SpaceInterface<typename S::GV, S::r, S::rC, typename S::R>, S>
{};


// from #include <dune/gdt/test/momentmodels/basisfunctions.hh>
template <class T>
struct is_hatfunction_basis : public std::false_type
{};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t dimFlux,
          EntropyType entropy>
struct is_hatfunction_basis<HatFunctionMomentBasis<DomainFieldType,
                                                   dimDomain,
                                                   RangeFieldType,
                                                   dimRange_or_refinements,
                                                   dimRangeCols,
                                                   dimFlux,
                                                   entropy>> : public std::true_type
{};

template <class T>
struct is_legendre_basis : public std::false_type
{};

template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols, EntropyType entropy>
struct is_legendre_basis<LegendreMomentBasis<DomainFieldType, RangeFieldType, order, dimRangeCols, entropy>>
  : public std::true_type
{};

template <class T>
struct is_partial_moment_basis : public std::false_type
{};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t dimFlux,
          size_t order,
          EntropyType entropy>
struct is_partial_moment_basis<PartialMomentBasis<DomainFieldType,
                                                  dimDomain,
                                                  RangeFieldType,
                                                  dimRange_or_refinements,
                                                  dimRangeCols,
                                                  dimFlux,
                                                  order,
                                                  entropy>> : public std::true_type
{};

template <class T>
struct is_1d_partial_moment_basis : public std::false_type
{};

template <class DomainFieldType,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t order,
          EntropyType entropy>
struct is_1d_partial_moment_basis<
    PartialMomentBasis<DomainFieldType, 1, RangeFieldType, dimRange_or_refinements, dimRangeCols, 1, order, entropy>>
  : public std::true_type
{};

template <class T>
struct is_3d_partial_moment_basis : public std::false_type
{};

template <class DomainFieldType,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t order,
          EntropyType entropy>
struct is_3d_partial_moment_basis<
    PartialMomentBasis<DomainFieldType, 3, RangeFieldType, dimRange_or_refinements, dimRangeCols, 3, order, entropy>>
  : public std::true_type
{};

template <class T>
struct is_spherical_harmonics_basis : public std::false_type
{};

template <class DomainFieldType,
          class RangeFieldType,
          size_t order,
          size_t fluxDim,
          bool only_positive,
          EntropyType entropy>
struct is_spherical_harmonics_basis<
    SphericalHarmonicsMomentBasis<DomainFieldType, RangeFieldType, order, fluxDim, only_positive, entropy>>
  : public std::true_type
{};

template <class T>
struct is_real_spherical_harmonics_basis : public std::false_type
{};

template <class DomainFieldType,
          class RangeFieldType,
          size_t order,
          size_t fluxDim,
          bool only_even,
          EntropyType entropy>
struct is_real_spherical_harmonics_basis<
    RealSphericalHarmonicsMomentBasis<DomainFieldType, RangeFieldType, order, fluxDim, only_even, entropy>>
  : public std::true_type
{};

template <class T>
struct is_full_moment_basis
{
  static constexpr bool value = is_legendre_basis<T>::value || is_spherical_harmonics_basis<T>::value
                                || is_real_spherical_harmonics_basis<T>::value;
};


#if 0
// from #include <dune/gdt/playground/spaces/restricted.hh>
template <class S, bool candidate = internal::is_restricted_space_helper<S>::is_candidate>
struct is_restricted_space
    : public std::is_base_of<RestrictedSpace<typename S::UnrestrictedSpaceType, typename S::RestrictionGridLayerType>,
                             S>
{
};

template <class S>
struct is_restricted_space<S, false> : public std::false_type
{
};


// from #include <dune/gdt/spaces/cg/interface.hh>
template <class S,
          bool space_candidate = internal::is_space_helper<S>::is_candidate,
          bool restricted = is_restricted_space<S>::value>
struct is_cg_space : public std::false_type
{
};

template <class S>
struct is_cg_space<S, true, false>
    : public std::is_base_of<CgSpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};

template <class S>
struct is_cg_space<S, true, true> : public is_cg_space<typename S::UnrestrictedSpaceType>
{
};


// from #include <dune/gdt/spaces/rt/interface.hh>
template <class S,
          bool space_candidate = internal::is_space_helper<S>::is_candidate,
          bool restricted = is_restricted_space<S>::value>
struct is_rt_space : public std::false_type
{
};

template <class S>
struct is_rt_space<S, true, false>
    : public std::is_base_of<RtSpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};

template <class S>
struct is_rt_space<S, true, true> : public is_rt_space<typename S::UnrestrictedSpaceType>
{
};


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
struct is_localizable_product : public std::is_base_of<LocalizableProductBase<typename T::GridLayerType,
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
                                                                      typename T::GridLayerType,
                                                                      typename T::SourceSpaceType,
                                                                      typename T::FieldType,
                                                                      T::pattern_type,
                                                                      typename T::OuterRangeSpaceType,
                                                                      typename T::OuterSourceSpaceType>,
                                                   T>
{
};

template <class T>
struct is_matrix_operator<T, false> : public std::false_type
{
};
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TYPE_TRAITS_HH
