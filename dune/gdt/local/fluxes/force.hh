// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_FORCE_HH
#define DUNE_GDT_LOCAL_FLUXES_FORCE_HH

#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"
#include "laxfriedrichs.hh"
#include "laxwendroff.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class LocalizableFunctionImp, class Traits>
class ForceLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp, class Traits>
class ForceLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class LocalizableFunctionImp>
class ForceLocalNumericalCouplingFluxTraits
  : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>
{
  typedef ForceLocalNumericalCouplingFluxTraits ThisType;

public:
  typedef ForceLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp, ThisType> derived_type;
}; // class ForceLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp>
class ForceLocalDirichletNumericalBoundaryFluxTraits
  : public LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                  BoundaryValueImp,
                                                                  LocalizableFunctionImp>
{
  typedef ForceLocalDirichletNumericalBoundaryFluxTraits ThisType;

public:
  typedef ForceLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                   BoundaryValueImp,
                                                   LocalizableFunctionImp,
                                                   ThisType>
      derived_type;
}; // class ForceLocalDirichletNumericalBoundaryFluxTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class LocalizableFunctionImp,
          class Traits = internal::ForceLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>>
class ForceLocalNumericalCouplingFlux : public LocalNumericalCouplingFluxInterface<Traits>
{
  typedef LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp> LaxFriedrichsLocalFluxType;
  typedef LaxWendroffLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp> LaxWendroffLocalFluxType;

public:
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename LocalizableFunctionType::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit ForceLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                           const XT::Common::Parameter& param,
                                           const LocalizableFunctionType& dx,
                                           const bool is_linear = false,
                                           const RangeFieldType alpha = dimDomain)
    : lax_friedrichs_flux_(analytical_flux, param, dx, false, is_linear, alpha)
    , lax_wendroff_flux_(analytical_flux, param, dx, is_linear, alpha)
  {}

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return lax_friedrichs_flux_.local_functions(entity);
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple_entity,
      const LocalfunctionTupleType& local_functions_tuple_neighbor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_neighbor,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    RangeType ret = lax_friedrichs_flux_.evaluate(local_functions_tuple_entity,
                                                  local_functions_tuple_neighbor,
                                                  local_source_entity,
                                                  local_source_neighbor,
                                                  intersection,
                                                  x_in_intersection_coords);
    ret += lax_wendroff_flux_.evaluate(local_functions_tuple_entity,
                                       local_functions_tuple_neighbor,
                                       local_source_entity,
                                       local_source_neighbor,
                                       intersection,
                                       x_in_intersection_coords);
    ret *= 0.5;
    return ret;
  } // RangeType evaluate(...) const

private:
  const LaxFriedrichsLocalFluxType lax_friedrichs_flux_;
  const LaxWendroffLocalFluxType lax_wendroff_flux_;
}; // class ForceLocalNumericalCouplingFlux

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          class Traits = internal::ForceLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                  BoundaryValueImp,
                                                                                  LocalizableFunctionImp>>
class ForceLocalDirichletNumericalBoundaryFlux : public LocalNumericalBoundaryFluxInterface<Traits>
{
  typedef LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, LocalizableFunctionImp>
      LaxFriedrichsLocalFluxType;
  typedef LaxWendroffLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, LocalizableFunctionImp>
      LaxWendroffLocalFluxType;

public:
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit ForceLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                    const BoundaryValueType& boundary_values,
                                                    const XT::Common::Parameter& param,
                                                    const LocalizableFunctionType& dx,
                                                    const bool is_linear = false,
                                                    const RangeFieldType alpha = dimDomain)
    : lax_friedrichs_flux_(analytical_flux, boundary_values, param, dx, false, is_linear, alpha)
    , lax_wendroff_flux_(analytical_flux, boundary_values, param, dx, is_linear, alpha)
  {}

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return lax_friedrichs_flux_.local_functions(entity);
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    RangeType ret = lax_friedrichs_flux_.evaluate(
        local_functions_tuple, local_source_entity, intersection, x_in_intersection_coords);
    ret +=
        lax_wendroff_flux_.evaluate(local_functions_tuple, local_source_entity, intersection, x_in_intersection_coords);
    ret *= 0.5;
    return ret;
  } // RangeType evaluate(...) const

private:
  const LaxFriedrichsLocalFluxType lax_friedrichs_flux_;
  const LaxWendroffLocalFluxType lax_wendroff_flux_;
}; // class ForceLocalDirichletNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_FORCE_HH
