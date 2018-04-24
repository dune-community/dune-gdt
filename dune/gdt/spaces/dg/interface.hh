// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_DG_INTERFACE_HH
#define DUNE_GDT_SPACES_DG_INTERFACE_HH

#include <dune/xt/common/type_traits.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {


static constexpr Backends default_dg_backend = default_space_backend;


template <class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1>
class DgSpaceInterface : public SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols>
{
  typedef SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols> BaseType;

public:
  typedef ImpTraits Traits;
  using typename BaseType::EntityType;
  using typename BaseType::PatternType;

  using BaseType::compute_pattern;

  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return BaseType::compute_face_and_volume_pattern(grd_layr, ansatz_space);
  }

  using BaseType::local_constraints;

  template <class S, size_t d, size_t r, size_t rC, class C, class R>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*other*/,
                         const EntityType& /*entity*/,
                         ConstraintsInterface<C>& /*ret*/) const
  {
    static_assert(AlwaysFalse<S>::value, "DG spaces do not implement constraints!");
  }

  static constexpr bool associates_data_with(int codim)
  {
    return codim == 0;
  }
}; // class DgSpaceInterface


namespace internal {


template <class S>
struct is_dg_space_helper
{
  DXTC_has_typedef_initialize_once(Traits);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);

  static const bool is_candidate = DXTC_has_typedef(Traits)<S>::value && DXTC_has_static_member(dimDomain)<S>::value
                                   && DXTC_has_static_member(dimRange)<S>::value
                                   && DXTC_has_static_member(dimRangeCols)<S>::value;
}; // class is_dg_space_helper


} // namespace internal


template <class S, bool candidate = internal::is_dg_space_helper<S>::is_candidate>
struct is_dg_space
    : public std::is_base_of<DgSpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};


template <class S>
struct is_dg_space<S, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_DG_INTERFACE_HH
