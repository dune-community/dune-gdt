// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TYPE_TRAITS_HH
#define DUNE_GDT_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace GDT {

// forwards
// from #include <dune/gdt/spaces/interface.hh>
enum class Backends;

enum class ChoosePattern;

template <class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols>
class SpaceInterface;

template <class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols>
class ProductSpaceInterface;

// from #include <dune/gdt/playground/spaces/restricted.hh>
template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedSpace;

// from #include <dune/gdt/spaces/cg/interface.hh>
template <class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols>
class CgSpaceInterface;

// from #include <dune/gdt/spaces/rt/interface.hh>
template <class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols>
class RtSpaceInterface;

// from #include <dune/gdt/local/fluxes/interfaces.hh>
template <class Traits>
class LocalNumericalCouplingFluxInterface;

template <class Traits>
class LocalNumericalBoundaryFluxInterface;

// from #include <dune/gdt/local/integrands/interfaces.hh>
template <class Traits, size_t numArguments>
class LocalVolumeIntegrandInterface;

template <class Traits, size_t numArguments>
class LocalFaceIntegrandInterface;

// from #include <dune/gdt/operators/interfaces.hh>
template <class Traits>
class OperatorInterface;

// from #include <dune/gdt/operators/base.hh>
template <class M, class RS, class GL, class SS, class F, ChoosePattern pt, class ORS, class OSS>
class MatrixOperatorBase;

// from #include <dune/gdt/discretefunction/default.hh>
template <class SpaceImp, class VectorImp>
class ConstDiscreteFunction;

template <class SpaceImp, class VectorImp>
class DiscreteFunction;

// from #include <dune/gdt/operators/fv/boundary.hh>
template <class EntityImp, class IntersectionImp, class RangeImp>
class LocalizableBoundaryValueInterface;

// from #include <dune/gdt/operators/fv/reconstructed_function.hh>
template <class GridLayerImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols>
class ReconstructedLocalizableFunction;

template <class GridLayerImp, class RangeImp, class SourceImp, class FieldImp>
class LocalizableProductBase;

template <class GridLayerImp, class SourceImp, class RangeImp>
class LocalizableOperatorBase;

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          size_t dimFlux>
class HatFunctionMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols>
class LegendreMomentBasis;

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols,
          size_t dimFlux,
          size_t order>
class PartialMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive>
class SphericalHarmonicsMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_even>
class RealSphericalHarmonicsMomentBasis;


namespace internal {


// helper structs
// from #include <dune/gdt/spaces/interface.hh>
template <class S>
struct is_space_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);

  static const bool is_candidate = DXTC_has_typedef(Traits)<S>::value && DXTC_has_static_member(dimDomain)<S>::value
                                   && DXTC_has_static_member(dimRange)<S>::value
                                   && DXTC_has_static_member(dimRangeCols)<S>::value;
}; // class is_space_helper


// from #include <dune/gdt/playground/spaces/restricted.hh>
template <class S>
struct is_restricted_space_helper
{
  DXTC_has_typedef_initialize_once(UnrestrictedSpaceType);
  DXTC_has_typedef_initialize_once(RestrictionGridLayerType);

  static const bool is_candidate =
      DXTC_has_typedef(UnrestrictedSpaceType)<S>::value && DXTC_has_typedef(RestrictionGridLayerType)<S>::value;
}; // class is_restricted_space_helper


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
  DXTC_has_typedef_initialize_once(GridLayerType);
  DXTC_has_typedef_initialize_once(RangeType);
  DXTC_has_typedef_initialize_once(SourceType);
  DXTC_has_typedef_initialize_once(FieldType);
  static const bool is_candidate = DXTC_has_typedef(GridLayerType)<Tt>::value && DXTC_has_typedef(RangeType)<Tt>::value
                                   && DXTC_has_typedef(SourceType)<Tt>::value && DXTC_has_typedef(FieldType)<Tt>::value;
};


template <class Tt>
struct is_matrix_operator_helper
{
  DXTC_has_typedef_initialize_once(MatrixType);
  DXTC_has_typedef_initialize_once(RangeSpaceType);
  DXTC_has_typedef_initialize_once(GridLayerType);
  DXTC_has_typedef_initialize_once(SourceSpaceType);
  DXTC_has_typedef_initialize_once(FieldType);
  DXTC_has_static_member_initialize_once(pattern_type);
  DXTC_has_typedef_initialize_once(OuterRangeSpaceType);
  DXTC_has_typedef_initialize_once(OuterSourceSpaceType);
  static const bool is_candidate =
      DXTC_has_typedef(MatrixType)<Tt>::value && DXTC_has_typedef(RangeSpaceType)<Tt>::value
      && DXTC_has_typedef(GridLayerType)<Tt>::value && DXTC_has_typedef(SourceSpaceType)<Tt>::value
      && DXTC_has_typedef(FieldType)<Tt>::value && DXTC_has_static_member(pattern_type)<Tt>::value
      && DXTC_has_typedef(OuterRangeSpaceType)<Tt>::value && DXTC_has_typedef(OuterSourceSpaceType)<Tt>::value;
};


template <class Tt>
struct is_const_discrete_function_helper
{
  DXTC_has_typedef_initialize_once(SpaceType);
  DXTC_has_typedef_initialize_once(VectorType);

  static const bool is_candidate = DXTC_has_typedef(SpaceType)<Tt>::value && DXTC_has_typedef(VectorType)<Tt>::value;
};


template <class Tt>
struct is_localizable_boundary_value_helper
{
  DXTC_has_typedef_initialize_once(EntityType);
  DXTC_has_typedef_initialize_once(IntersectionType);
  DXTC_has_typedef_initialize_once(RangeType);

  static const bool is_candidate = DXTC_has_typedef(EntityType)<Tt>::value
                                   && DXTC_has_typedef(IntersectionType)<Tt>::value
                                   && DXTC_has_typedef(RangeType)<Tt>::value;
};


template <class Tt>
struct is_reconstructed_localizable_function_helper
{
  DXTC_has_typedef_initialize_once(GridLayerType);
  DXTC_has_typedef_initialize_once(DomainFieldType);
  DXTC_has_typedef_initialize_once(RangeFieldType);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);
  static const bool is_candidate =
      DXTC_has_typedef(GridLayerType)<Tt>::value && DXTC_has_typedef(DomainFieldType)<Tt>::value
      && DXTC_has_typedef(RangeFieldType)<Tt>::value && DXTC_has_static_member(dimDomain)<Tt>::value
      && DXTC_has_static_member(dimRange)<Tt>::value && DXTC_has_static_member(dimRangeCols)<Tt>::value;
};


} // namespace internal


// actual structs
// from #include <dune/gdt/spaces/interface.hh>
template <class S, bool candidate = internal::is_space_helper<S>::is_candidate>
struct is_space
    : public std::is_base_of<SpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};

template <class S>
struct is_space<S, false> : public std::false_type
{
};

template <class S, bool candidate = internal::is_space_helper<S>::is_candidate>
struct is_product_space
    : public std::is_base_of<ProductSpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};

template <class S>
struct is_product_space<S, false> : public std::false_type
{
};


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

// from #include <dune/gdt/local/fluxes/interfaces.hh>
template <class T, bool candidate = internal::is_operator_helper<T>::is_candidate>
struct is_local_numerical_coupling_flux : std::is_base_of<LocalNumericalCouplingFluxInterface<typename T::Traits>, T>
{
};

template <class T, bool candidate = internal::is_operator_helper<T>::is_candidate>
struct is_local_numerical_boundary_flux : std::is_base_of<LocalNumericalBoundaryFluxInterface<typename T::Traits>, T>
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


// from #include <dune/gdt/discretefunction/default.hh>
template <class T, bool candidate = internal::is_const_discrete_function_helper<T>::is_candidate>
struct is_const_discrete_function
    : public std::is_base_of<ConstDiscreteFunction<typename T::SpaceType, typename T::VectorType>, T>
{
};

template <class T>
struct is_const_discrete_function<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_const_discrete_function_helper<T>::is_candidate>
struct is_discrete_function : public std::is_base_of<DiscreteFunction<typename T::SpaceType, typename T::VectorType>, T>
{
};

template <class T>
struct is_discrete_function<T, false> : public std::false_type
{
};


// from #include <dune/gdt/operators/fv/boundary.hh>
template <class T, bool candidate = internal::is_localizable_boundary_value_helper<T>::is_candidate>
struct is_localizable_boundary_value
    : public std::is_base_of<LocalizableBoundaryValueInterface<typename T::EntityType,
                                                               typename T::IntersectionType,
                                                               typename T::RangeType>,
                             T>
{
};

template <class T>
struct is_localizable_boundary_value<T, false> : public std::false_type
{
};

// from #include <dune/gdt/operators/fv/reconstructed_function.hh>
template <class T, bool candidate = internal::is_reconstructed_localizable_function_helper<T>::is_candidate>
struct is_reconstructed_localizable_function
    : public std::is_base_of<ReconstructedLocalizableFunction<typename T::GridLayerType,
                                                              typename T::DomainFieldType,
                                                              T::dimDomain,
                                                              typename T::RangeFieldType,
                                                              T::dimRange,
                                                              T::dimRangeCols>,
                             T>
{
};

template <class T>
struct is_reconstructed_localizable_function<T, false> : public std::false_type
{
};

// from #include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>
template <class T>
struct is_hatfunction_basis : public std::false_type
{
};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t dimFlux>
struct is_hatfunction_basis<HatFunctionMomentBasis<DomainFieldType,
                                                   dimDomain,
                                                   RangeFieldType,
                                                   dimRange_or_refinements,
                                                   dimRangeCols,
                                                   dimFlux>> : public std::true_type
{
};

template <class T>
struct is_legendre_basis : public std::false_type
{
};

template <class DomainFieldType, class RangeFieldType, size_t order, size_t dimRangeCols>
struct is_legendre_basis<LegendreMomentBasis<DomainFieldType, RangeFieldType, order, dimRangeCols>>
    : public std::true_type
{
};

template <class T>
struct is_partial_moment_basis : public std::false_type
{
};

template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          size_t dimRangeCols,
          size_t dimFlux,
          size_t order>
struct is_partial_moment_basis<PartialMomentBasis<DomainFieldType,
                                                  dimDomain,
                                                  RangeFieldType,
                                                  dimRange_or_refinements,
                                                  dimRangeCols,
                                                  dimFlux,
                                                  order>> : public std::true_type
{
};

template <class T>
struct is_spherical_harmonics_basis : public std::false_type
{
};

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive>
class SphericalHarmonicsMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive>
struct is_spherical_harmonics_basis<SphericalHarmonicsMomentBasis<DomainFieldType,
                                                                  RangeFieldType,
                                                                  order,
                                                                  fluxDim,
                                                                  only_positive>> : public std::true_type
{
};

template <class T>
struct is_real_spherical_harmonics_basis : public std::false_type
{
};

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_positive>
class SphericalHarmonicsMomentBasis;

template <class DomainFieldType, class RangeFieldType, size_t order, size_t fluxDim, bool only_even>
struct is_real_spherical_harmonics_basis<RealSphericalHarmonicsMomentBasis<DomainFieldType,
                                                                           RangeFieldType,
                                                                           order,
                                                                           fluxDim,
                                                                           only_even>> : public std::true_type
{
};

template <class T>
struct is_full_moment_basis
{
  static constexpr bool value = is_legendre_basis<T>::value || is_spherical_harmonics_basis<T>::value
                                || is_real_spherical_harmonics_basis<T>::value;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TYPE_TRAITS_HH
