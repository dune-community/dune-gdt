// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#ifndef DUNE_GDT_SPACES_RESTRICTED_HH
#define DUNE_GDT_SPACES_RESTRICTED_HH

#include <sstream>

#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include "mapper/restricted.hh"

namespace Dune {
namespace GDT {


// forward
template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedSpace;


namespace internal {


template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedSpaceTraits
{
  static_assert(XT::Grid::is_layer<RestrictionGridLayer>::value, "");
  static_assert(is_space<UnrestrictedSpace>::value, "");

public:
  typedef RestrictedSpace<UnrestrictedSpace, RestrictionGridLayer> derived_type;
  static const int polOrder = UnrestrictedSpace::polOrder;
  static const bool continuous = UnrestrictedSpace::continuous;
  typedef UnrestrictedSpace BackendType;
  typedef RestrictedMapper<UnrestrictedSpace, RestrictionGridLayer> MapperType;
  typedef typename UnrestrictedSpace::BaseFunctionSetType BaseFunctionSetType;
  typedef typename UnrestrictedSpace::DofCommunicatorType DofCommunicatorType;
  typedef RestrictionGridLayer GridLayerType;
  typedef typename UnrestrictedSpace::RangeFieldType RangeFieldType;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static constexpr const GDT::Backends backend_type{UnrestrictedSpace::backend_type};
}; // class RestrictedSpaceTraits


} // namespace internal


template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedSpace : public SpaceInterface<internal::RestrictedSpaceTraits<UnrestrictedSpace, RestrictionGridLayer>,
                                              UnrestrictedSpace::dimDomain,
                                              UnrestrictedSpace::dimRange,
                                              UnrestrictedSpace::dimRangeCols>
{
  typedef SpaceInterface<internal::RestrictedSpaceTraits<UnrestrictedSpace, RestrictionGridLayer>,
                         UnrestrictedSpace::dimDomain,
                         UnrestrictedSpace::dimRange,
                         UnrestrictedSpace::dimRangeCols>
      BaseType;
  typedef RestrictedSpace<UnrestrictedSpace, RestrictionGridLayer> ThisType;

public:
  typedef internal::RestrictedSpaceTraits<UnrestrictedSpace, RestrictionGridLayer> Traits;
  typedef UnrestrictedSpace UnrestrictedSpaceType; //       These are mainly here to detect this space type from the
  typedef RestrictionGridLayer RestrictionGridLayerType; // outside (see also type_traits.hh).
  using typename BaseType::MapperType;
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::DofCommunicatorType;
  using typename BaseType::PatternType;
  using typename BaseType::DomainType;

  RestrictedSpace(const UnrestrictedSpace& unrestricted_space, RestrictionGridLayer restriction_grid_layer)
    : unrestricted_space_(unrestricted_space)
    , grid_layer_(restriction_grid_layer)
    , mapper_(unrestricted_space_, grid_layer_)
  {
  }

  RestrictedSpace(const ThisType& other) = default;
  RestrictedSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    return unrestricted_space_;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    check_entity(entity);
    return unrestricted_space_.base_function_set(entity);
  }

  DofCommunicatorType& dof_communicator() const
  {
    return unrestricted_space_.dof_communicator();
  }

  using BaseType::local_constraints;

  template <class S, size_t d, size_t r, size_t rC, class C>
  void local_constraints(const SpaceInterface<S, d, r, rC>& ansatz_space,
                         const EntityType& entity,
                         ConstraintsInterface<C>& ret) const
  {
    check_entity(entity);
    return unrestricted_space_.local_constraints(ansatz_space, entity, ret);
  }

  using BaseType::compute_pattern;

  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_pattern(const GL& /*grd_layr*/, const SpaceInterface<S, d, r, rC>& /*ansatz_space*/) const
  {
    DUNE_THROW(NotImplemented, "Yet");
  }

  // if we are CG
  template <class E>
  std::vector<DomainType> lagrange_points(const E& entity) const
  {
    check_entity(entity);
    return unrestricted_space_.lagrange_points(entity);
  }

  template <class E, class I>
  std::set<size_t> local_dirichlet_DoFs(const E& entity, const XT::Grid::BoundaryInfo<I>& boundaryInfo) const
  {
    check_entity(entity);
    return unrestricted_space_.local_dirichlet_DoFs(entity, boundaryInfo);
  }

  // if we are RT
  template <class E>
  std::vector<size_t> local_DoF_indices(const E& entity) const
  {
    check_entity(entity);
    return unrestricted_space_.local_DoF_indices(entity);
  }

private:
  void check_entity(const EntityType& entity) const
  {
    if (grid_layer_.indexSet().contains(entity))
      return;
    if (unrestricted_space_.grid_layer().indexSet().contains(entity))
      DUNE_THROW(restricted_space_error,
                 "Entity not contained in restriction grid layer, but contained in the unrestricted grid layer "
                     << "with index "
                     << unrestricted_space_.grid_layer().indexSet().index(entity)
                     << "!");
    else
      DUNE_THROW(restricted_space_error,
                 "Entity neither contained in restriction grid layer nor in the unrestricted grid layer!");

  } // ... check_entity(...)

  const UnrestrictedSpace unrestricted_space_;
  const GridLayerType grid_layer_;
  const MapperType mapper_;
}; // class RestrictedSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RESTRICTED_HH
