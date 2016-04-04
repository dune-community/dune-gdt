// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_FV_INTERFACE_HH
#define DUNE_GDT_SPACES_FV_INTERFACE_HH

#include <dune/stuff/common/type_utils.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {


template< class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1 >
class FVInterface
  : public SpaceInterface< Traits, domainDim, rangeDim, rangeDimCols >
{
  typedef SpaceInterface< Traits, domainDim, rangeDim, rangeDimCols > BaseType;
public:
  using typename BaseType::EntityType;
  using typename BaseType::PatternType;

  using BaseType::compute_pattern;

  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return BaseType::compute_face_and_volume_pattern(local_grid_view, ansatz_space);
  }

  using BaseType::local_constraints;

  template< class S, size_t d, size_t r, size_t rC, class C, class R >
  void local_constraints(const SpaceInterface< S, d, r, rC >& /*other*/,
                         const EntityType& /*entity*/,
                         Spaces::ConstraintsInterface< C >& /*ret*/) const
  {
    static_assert(AlwaysFalse< S >::value, "FV spaces do not implement constraints!");
  }
}; // class FVInterface


} // namespace Spaces
namespace internal {


template< class S >
struct is_fv_space_helper
{
  DSC_has_typedef_initialize_once(Traits)
  DSC_has_static_member_initialize_once(dimDomain)
  DSC_has_static_member_initialize_once(dimRange)
  DSC_has_static_member_initialize_once(dimRangeCols)

  static const bool is_candidate = DSC_has_typedef(Traits)< S >::value
                                   && DSC_has_static_member(dimDomain)< S >::value
                                   && DSC_has_static_member(dimRange)< S >::value
                                   && DSC_has_static_member(dimRangeCols)< S >::value;
}; // class is_fv_space_helper


} // namespace internal


template< class S, bool candidate = internal::is_fv_space_helper< S >::is_candidate >
struct is_fv_space
  : public std::is_base_of< Spaces::FVInterface< typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols >
                          , S >
{};

template< class S >
struct is_fv_space< S, false >
  : public std::false_type
{};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_INTERFACE_HH
