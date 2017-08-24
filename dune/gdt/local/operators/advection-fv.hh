// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH

#include <functional>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class SpaceType>
class LocalAdvectionFvInnerOperator;


namespace internal {


template <class SpaceType>
class LocalAdvectionFvInnerOperatorTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgInnerOperator instead!");

public:
  using derived_type = LocalAdvectionFvInnerOperator<SpaceType>;
};


} // namespace internal


/**
 * \note Presumes that the basis evaluates to 1.
 */
template <class SpaceType>
class LocalAdvectionFvInnerOperator
    : public LocalCouplingOperatorInterface<internal::LocalAdvectionFvInnerOperatorTraits<SpaceType>>
{
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = SpaceType::dimRangeCols;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>;

public:
  using FluxType = XT::Functions::LocalizableFluxFunctionInterface<E, D, d, StateType, 0, R, d>;
  using NumericalFluxType = std::function<R(const typename FluxType::LocalfunctionType&,
                                            const typename StateType::RangeType&,
                                            const typename StateType::RangeType&,
                                            const typename StateType::DomainType&)>;

  LocalAdvectionFvInnerOperator(const FluxType& flux, NumericalFluxType numerical_flux)
    : flux_(flux)
    , numerical_flux_(numerical_flux)
  {
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range_entity.entity();
    const auto& neighbor = local_range_neighbor.entity();
    const auto normal = intersection.centerUnitOuterNormal();
    const auto u_inside = source.local_discrete_function(entity);
    const auto u_outside = source.local_discrete_function(neighbor);
    assert(u_inside->vector().size() == 1);
    assert(u_outside->vector().size() == 1);
    assert(local_range_entity.vector().size() == 1);
    assert(local_range_neighbor.vector().size() == 1);
    // If we do not assume that the local basis evaluates to 1 (but is still constant), we need to
    // * use u_inside->evaluate(x_entity) instead of u_inside->vector().get(0)
    // * use u_outside->evaluate(x_neighbor) instead of u_outside->vector().get(0)
    // * \int_entity basis^2 \dx instead of h
    // where x_entity and x_neighbor are the corresponding coordinates of the intersections midpoint.
    const auto g =
        numerical_flux_(*flux_.local_function(entity), u_inside->vector().get(0), u_outside->vector().get(0), normal);
    const auto h = local_range_entity.entity().geometry().volume();
    local_range_entity.vector().add(0, (g * intersection.geometry().volume()) / h);
    local_range_neighbor.vector().add(0, (-1.0 * g * intersection.geometry().volume()) / h);
  }

private:
  const FluxType& flux_;
  const NumericalFluxType numerical_flux_;
}; // class LocalAdvectionFvInnerOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
