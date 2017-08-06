// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FLUXES_FORCE_HH
#define DUNE_GDT_LOCAL_FLUXES_FORCE_HH

#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"
#include "laxfriedrichs.hh"
#include "laxwendroff.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp,
          class LocalizableFunctionImp>
class ForceLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp>
class ForceLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class LocalizableFunctionImp>
class ForceLocalNumericalCouplingFluxTraits
    : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>
{
public:
  typedef ForceLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp> derived_type;
}; // class ForceLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp>
class ForceLocalDirichletNumericalBoundaryFluxTraits
    : public LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                    BoundaryValueImp,
                                                                    LocalizableFunctionImp>
{
public:
  typedef ForceLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                   BoundaryValueImp,
                                                   LocalizableFunctionImp>
      derived_type;
}; // class ForceLocalDirichletNumericalBoundaryFluxTraits


} // namespace internal


/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 *
 *  The Lax-Friedrichs flux is an approximation to the integral
 *  \int_{S_{ij}} \mathbf{F}(\mathbf{u}) \cdot \mathbf{n}_{ij},
 *  where S_{ij} is the intersection between the entities i and j, \mathbf{F}(\mathbf{u}) is the analytical flux
 *  (evaluated at \mathbf{u}) and \mathbf{n}_{ij} is the unit outer normal of S_{ij}.
 *  The Lax-Friedrichs flux takes the form
 *  \mathbf{g}_{ij}^{LF}(\mathbf{u}_i, \mathbf{u}_j)
 *  = \int_{S_{ij}} \frac{1}{2}(\mathbf{F}(\mathbf{u}_i) + \mathbf{F}(\mathbf{u}_j) \cdot \mathbf{n}_{ij}
 *  - \frac{1}{\alpha_i \lambda_{ij}} (\mathbf{u}_j - \mathbf{u}_i),
 *  where \alpha_i is the number of neighbors (i.e. intersections) of the entity i and lambda_{ij} is a local
 *  constant fulfilling
 *  \lambda_{ij} \sup_{\mathbf{u}} (\mathbf{F}(\mathbf{u} \cdot \mathbf{n}_{ij})^\prime \leq 1.
 *  The integration is done numerically and implemented in the LocalCouplingFvOperator. This class implements
 *  the evaluation of the integrand. As we are restricting ourselves to axis-parallel cubic grids, only one component of
 *  \mathbf{n}_{ij} is non-zero, denote this component by k. Then the Lax-Friedrichs flux evaluation reduces to
 *  \frac{1}{2}(\mathbf{f}^k(\mathbf{u}_i) + \mathbf{f}^k(\mathbf{u}_j) n_{ij,k}
 *  - \frac{1}{\alpha_i \lambda_{ij}} (\mathbf{u}_j - \mathbf{u}_i),
 *  where \mathbf{f}^k is the k-th column of the analytical flux.
 *  For the classical Lax-Friedrichs flux, \lambda_{ij} is chosen as dt/dx_i, where dt is the current time
 *  step length and dx_i is the width of entity i. This fulfills the equation above as long as the CFL condition
 *  is fulfilled.
 *  The local Lax-Friedrichs flux can be chosen by setting \param use_local to true, here \lambda_{ij} is chosen
 *  as the inverse of the maximal eigenvalue of \mathbf{f}^k(\mathbf{u}_i) and \mathbf{f}^k(\mathbf{u}_j). In this
 *  case, you should also specify whether your analytical flux is linear by setting \param is_linear, which avoids
 *  recalculating the eigenvalues on every intersection in the linear case.
 *  You can also provide a user-defined \param lambda that is used as \lambda_{ij} on all intersections. You need to set
 *  use_local to false, otherwise lambda will not be used.
 */
template <class AnalyticalFluxImp, class LocalizableFunctionImp>
class ForceLocalNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::ForceLocalNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                 LocalizableFunctionImp>>
{
  typedef LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp>
      LaxFriedrichsLocalFluxType;
  typedef LaxWendroffLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp>
      LaxWendroffLocalFluxType;

public:
  typedef internal::ForceLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>
      Traits;
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
  {
  }

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

/**
*  \brief  Lax-Friedrichs flux evaluation for Dirichlet boundary intersections.
*  \see    LaxFriedrichsLocalNumericalCouplingFlux
*/
template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp>
class ForceLocalDirichletNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     ForceLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                    BoundaryValueImp,
                                                                                                    LocalizableFunctionImp>>
{
  typedef LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                           BoundaryValueImp,
                                                           LocalizableFunctionImp>
      LaxFriedrichsLocalFluxType;
  typedef LaxWendroffLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                         BoundaryValueImp,
                                                         LocalizableFunctionImp>
      LaxWendroffLocalFluxType;

public:
  typedef internal::ForceLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                   BoundaryValueImp,
                                                                   LocalizableFunctionImp>
      Traits;
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
  {
  }

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
