// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_LAXWENDROFF_HH
#define DUNE_GDT_LOCAL_FLUXES_LAXWENDROFF_HH

#include <tuple>

#include "interfaces.hh"
#include "laxfriedrichs.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class LocalizableFunctionImp, class Traits>
class LaxWendroffLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          class LocalizableFunctionImp,
          class Traits>
class LaxWendroffLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class LocalizableFunctionImp>
class LaxWendroffLocalNumericalCouplingFluxTraits
    : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>
{
  typedef LaxWendroffLocalNumericalCouplingFluxTraits ThisType;

public:
  typedef LaxWendroffLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionImp, ThisType> derived_type;
}; // class LaxWendroffLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class BoundaryInfoImp, class LocalizableFunctionImp>
class LaxWendroffLocalDirichletNumericalBoundaryFluxTraits
    : public LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                    BoundaryValueImp,
                                                                    BoundaryInfoImp,
                                                                    LocalizableFunctionImp>
{
  typedef LaxWendroffLocalDirichletNumericalBoundaryFluxTraits ThisType;

public:
  typedef LaxWendroffLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                         BoundaryValueImp,
                                                         BoundaryInfoImp,
                                                         LocalizableFunctionImp,
                                                         ThisType>
      derived_type;
}; // class LaxWendroffLocalDirichletNumericalBoundaryFluxTraits

template <class Traits>
class LaxWendroffFluxImplementation
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
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> JacobianRangeType;

  explicit LaxWendroffFluxImplementation(const AnalyticalFluxType& analytical_flux,
                                         XT::Common::Parameter param,
                                         const bool is_linear = false,
                                         const RangeFieldType alpha = dimDomain,
                                         const bool boundary = false)
    : analytical_flux_(analytical_flux)
    , param_inside_(param)
    , param_outside_(param)
    , dt_(param.get("dt")[0])
    , is_linear_(is_linear)
    , alpha_(alpha)
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

    const auto& local_flux_inside = std::get<0>(local_functions_tuple_entity);
    const auto& local_flux_outside = std::get<0>(local_functions_tuple_neighbor);

    // calculate flux evaluation as
    // ret = f((u_i + u_j)*0.5 - (f_u_j - f_u_i)*0.5*dt/dx*dimDomain)*n_ij[direction]
    const RangeFieldType dx = std::get<1>(local_functions_tuple_entity)->evaluate(x_in_inside_coords)[0];
    RangeType ret = u_i;
    ret += u_j;
    ret *= 0.5;
    RangeType second_part = local_flux_outside->evaluate_col(direction, x_in_outside_coords, u_j, param_outside_);
    second_part -= local_flux_inside->evaluate_col(direction, x_in_inside_coords, u_i, param_inside_);
    second_part *= n_ij[direction] * 0.5 * dt_ / dx * alpha_;
    ret -= second_part;
    ret = local_flux_inside->evaluate_col(direction, x_in_inside_coords, ret, param_inside_);
    ret *= n_ij[direction];
    return ret;
  } // ... evaluate(...)

  const AnalyticalFluxType& analytical_flux() const
  {
    return analytical_flux_;
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;
  const double dt_;
  const bool is_linear_;
  const RangeFieldType alpha_;
}; // class LaxWendroffFluxImplementation<...>


} // namespace internal


template <class AnalyticalFluxImp,
          class LocalizableFunctionImp,
          class Traits =
              internal::LaxWendroffLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp>>
class LaxWendroffLocalNumericalCouplingFlux : public LocalNumericalCouplingFluxInterface<Traits>
{
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

  explicit LaxWendroffLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                                 const XT::Common::Parameter& param,
                                                 const LocalizableFunctionType& dx,
                                                 const bool is_linear = false,
                                                 const RangeFieldType alpha = dimDomain)
    : dx_(dx)
    , implementation_(analytical_flux, param, is_linear, alpha, false)
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
  const internal::LaxWendroffFluxImplementation<Traits> implementation_;
}; // class LaxWendroffLocalNumericalCouplingFlux

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          class LocalizableFunctionImp,
          class Traits = internal::LaxWendroffLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                        BoundaryValueImp,
                                                                                        BoundaryInfoImp,
                                                                                        LocalizableFunctionImp>>
class LaxWendroffLocalDirichletNumericalBoundaryFlux : public LocalNumericalBoundaryFluxInterface<Traits>
{
  typedef LocalNumericalBoundaryFluxInterface<Traits> InterfaceType;

public:
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;
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

  explicit LaxWendroffLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                          const BoundaryValueType& boundary_values,
                                                          const BoundaryInfoType& boundary_info,
                                                          const XT::Common::Parameter& param,
                                                          const LocalizableFunctionType& dx,
                                                          const bool is_linear = false,
                                                          const RangeFieldType alpha = dimDomain)
    : InterfaceType(boundary_values, boundary_info)
    , dx_(dx)
    , implementation_(analytical_flux, param, is_linear, alpha, true)
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
    const auto values =
        InterfaceType::get_values(local_functions_tuple, local_source_entity, intersection, x_in_intersection_coords);
    return implementation_.evaluate(local_functions_tuple,
                                    local_functions_tuple,
                                    intersection,
                                    x_in_intersection_coords,
                                    std::get<2>(values),
                                    std::get<2>(values),
                                    std::get<0>(values),
                                    std::get<1>(values));
  } // RangeType evaluate(...) const

private:
  const LocalizableFunctionType& dx_;
  const internal::LaxWendroffFluxImplementation<Traits> implementation_;
  using InterfaceType::boundary_values_;
}; // class LaxWendroffLocalDirichletNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_LAXWENDROFF_HH
