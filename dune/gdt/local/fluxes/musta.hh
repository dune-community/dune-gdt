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

#ifndef DUNE_GDT_LOCAL_FLUXES_MUSTA_HH
#define DUNE_GDT_LOCAL_FLUXES_MUSTA_HH

#include <tuple>

#include "interfaces.hh"
#include "laxfriedrichs.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp,
          class LocalizableFunctionImp,
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::RangeFieldType,
                                                    AnalyticalFluxImp::dimRange,
                                                    AnalyticalFluxImp::dimRangeCols>>
class MustaLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          class EigenSolverImp = DefaultEigenSolver<typename AnalyticalFluxImp::RangeFieldType,
                                                    AnalyticalFluxImp::dimRange,
                                                    AnalyticalFluxImp::dimRangeCols>>
class MustaLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class LocalizableFunctionImp, class EigenSolverImp>
class MustaLocalNumericalCouplingFluxTraits
    : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, EigenSolverImp>
{
public:
  typedef MustaLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp, EigenSolverImp> derived_type;
}; // class MustaLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp, class EigenSolverImp>
class MustaLocalDirichletNumericalBoundaryFluxTraits
    : public LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                    BoundaryValueImp,
                                                                    LocalizableFunctionImp,
                                                                    EigenSolverImp>
{
public:
  typedef MustaLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                   BoundaryValueImp,
                                                   LocalizableFunctionImp,
                                                   EigenSolverImp>
      derived_type;
}; // class MustaLocalDirichletNumericalBoundaryFluxTraits

template <class Traits>
class MustaFluxImplementation
{
public:
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::AnalyticalFluxLocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename Traits::EigenSolverType EigenSolverType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> JacobianRangeType;

  explicit MustaFluxImplementation(const AnalyticalFluxType& analytical_flux,
                                   XT::Common::Parameter param,
                                   const size_t num_stages,
                                   const bool boundary)
    : analytical_flux_(analytical_flux)
    , param_inside_(param)
    , param_outside_(param)
    , dt_(param.get("dt")[0])
    , num_stages_(num_stages)
  {
    param_inside_.set("boundary", {0.}, true);
    param_outside_.set("boundary", {double(boundary)}, true);
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple_entity,
                     const LocalfunctionTupleType& local_functions_tuple_neighbor,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords,
                     const DomainType& x_in_inside_coords,
                     const DomainType& x_in_outside_coords,
                     const RangeType& u_i,
                     const RangeType& u_j) const
  {
    // find direction of unit outer normal
    auto n_ij = intersection.unitOuterNormal(x_in_intersection_coords);
    size_t direction = intersection.indexInInside() / 2;
    bool positive_dir = intersection.indexInInside() % 2;

    // get data on left and right entity
    const RangeType& u_R = positive_dir ? u_j : u_i;
    const RangeType& u_L = positive_dir ? u_i : u_j;
    const DomainType& x_in_R_coords = positive_dir ? x_in_outside_coords : x_in_inside_coords;
    const DomainType& x_in_L_coords = positive_dir ? x_in_inside_coords : x_in_outside_coords;
    const auto& local_flux_R =
        std::get<0>(positive_dir ? local_functions_tuple_neighbor : local_functions_tuple_entity);
    const auto& local_flux_L =
        std::get<0>(positive_dir ? local_functions_tuple_entity : local_functions_tuple_neighbor);
    const auto& local_flux_inside = positive_dir ? local_flux_L : local_flux_R;
    const auto& param_L = positive_dir ? param_outside_ : param_inside_;
    const auto& param_R = positive_dir ? param_inside_ : param_outside_;

    const RangeFieldType dx = std::get<1>(local_functions_tuple_entity)->evaluate(x_in_inside_coords)[0];

    RangeType u_L_l = u_L;
    RangeType u_R_l = u_R;
    RangeType ret;
    for (size_t ll = 0; ll <= num_stages_; ++ll)
      calculate_stage(ret,
                      u_L_l,
                      u_R_l,
                      *local_flux_L,
                      *local_flux_R,
                      *local_flux_inside,
                      x_in_L_coords,
                      x_in_R_coords,
                      x_in_inside_coords,
                      param_L,
                      param_R,
                      dx,
                      direction,
                      ll != num_stages_);
    ret *= n_ij[direction];
    return ret;
  }

  const AnalyticalFluxType& analytical_flux() const
  {
    return analytical_flux_;
  }

private:
  void calculate_stage(RangeType& ret,
                       RangeType& u_L,
                       RangeType& u_R,
                       const AnalyticalFluxLocalfunctionType& local_flux_L,
                       const AnalyticalFluxLocalfunctionType& local_flux_R,
                       const AnalyticalFluxLocalfunctionType& local_flux_inside,
                       const DomainType& x_in_L_coords,
                       const DomainType& x_in_R_coords,
                       const DomainType& x_in_inside_coords,
                       const XT::Common::Parameter& param_L,
                       const XT::Common::Parameter& param_R,
                       const RangeFieldType& dx,
                       const size_t direction,
                       const bool update) const
  {
    // get data of current_stage
    const RangeType f_L = local_flux_L.evaluate_col(direction, x_in_L_coords, u_L, param_L);
    const RangeType f_R = local_flux_R.evaluate_col(direction, x_in_R_coords, u_R, param_R);

    // calculate u^{1/2}
    RangeType u_half = u_L;
    u_half += u_R;
    RangeType f_diff = f_R;
    f_diff -= f_L;
    f_diff *= dx / dt_;
    u_half -= f_diff;
    u_half *= 0.5;

    // calculate f_M
    const RangeType f_M = local_flux_inside.evaluate_col(direction, x_in_inside_coords, u_half, param_inside_);

    // calculate force flux
    auto& f_force = f_diff;
    f_force = f_M;
    f_force *= 2.;
    f_force += f_L;
    f_force += f_R;
    auto& u_diff = u_half;
    u_diff = u_R;
    u_diff -= u_L;
    u_diff *= dx / dt_;
    f_force -= u_diff;
    f_force *= 0.25;

    // update data
    if (update) {
      auto update_L = f_force;
      update_L -= f_L;
      update_L *= dt_ / dx;
      u_L -= update_L;
      auto& update_R = update_L;
      update_R = f_R;
      update_R -= f_force;
      update_R *= dt_ / dx;
      u_R -= update_R;
    } else {
      ret = f_force;
    }
  }

  const AnalyticalFluxType& analytical_flux_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;
  const double dt_;
  const size_t num_stages_;
}; // class MustaFluxImplementation<...>


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
template <class AnalyticalFluxImp, class LocalizableFunctionImp, class EigenSolverImp>
class MustaLocalNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::MustaLocalNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                 LocalizableFunctionImp,
                                                                                                 EigenSolverImp>>
{
public:
  typedef internal::MustaLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, EigenSolverImp>
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

  explicit MustaLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                           const XT::Common::Parameter& param,
                                           const LocalizableFunctionType& dx,
                                           const size_t num_stages = 2)
    : dx_(dx)
    , implementation_(analytical_flux, param, num_stages, false)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(implementation_.analytical_flux().local_function(entity), dx_.local_function(entity));
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
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const auto x_in_outside_coords = intersection.geometryInOutside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = local_source_neighbor.evaluate(x_in_outside_coords);
    return implementation_.evaluate(local_functions_tuple_entity,
                                    local_functions_tuple_neighbor,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_outside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const LocalizableFunctionType& dx_;
  const internal::MustaFluxImplementation<Traits> implementation_;
}; // class MustaLocalNumericalCouplingFlux

/**
*  \brief  MUSTA flux evaluation for Dirichlet boundary intersections.
*  \see    MustaLocalNumericalCouplingFlux
*/
template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp, class EigenSolverImp>
class MustaLocalDirichletNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     MustaLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                    BoundaryValueImp,
                                                                                                    LocalizableFunctionImp,
                                                                                                    EigenSolverImp>>
{
public:
  typedef internal::MustaLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                   BoundaryValueImp,
                                                                   LocalizableFunctionImp,
                                                                   EigenSolverImp>
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

  explicit MustaLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                    const BoundaryValueType& boundary_values,
                                                    const XT::Common::Parameter& param,
                                                    const LocalizableFunctionType& dx,
                                                    const size_t num_stages = 2)
    : boundary_values_(boundary_values)
    , dx_(dx)
    , implementation_(analytical_flux, param, num_stages, true)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(implementation_.analytical_flux().local_function(entity),
                           dx_.local_function(entity),
                           boundary_values_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const

  {
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = std::get<2>(local_functions_tuple)->evaluate(x_in_inside_coords);
    return implementation_.evaluate(local_functions_tuple,
                                    local_functions_tuple,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_inside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const BoundaryValueType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const internal::MustaFluxImplementation<Traits> implementation_;
}; // class MustaLocalDirichletNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_MUSTA_HH
