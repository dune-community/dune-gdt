// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH
#define DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/crtp.hh>

#include <dune/xt/grid/boundaryinfo/types.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace GDT {
namespace internal {


class IsNumericalCouplingFlux
{
};

class IsNumericalBoundaryFlux
{
};

class IsAnalyticalFlux
{
};

class IsRHSEvaluation
{
};


} // namespace internal


template <class Traits>
class LocalNumericalCouplingFluxInterface
    : public XT::CRTPInterface<LocalNumericalCouplingFluxInterface<Traits>, Traits>,
      internal::IsNumericalCouplingFlux
{
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;

public:
  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity))
    return this->as_imp().local_functions(entity);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  auto evaluate(const LocalfunctionTupleType& local_functions_tuple_entity,
                const LocalfunctionTupleType& local_functions_tuple_neighbor,
                const XT::Functions::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_entity,
                const XT::Functions::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_neighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& x_intersection) const ->
      typename XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple_entity,
                                       local_functions_tuple_neighbor,
                                       local_source_entity,
                                       local_source_neighbor,
                                       intersection,
                                       x_intersection));
    this->as_imp().evaluate(local_functions_tuple_entity,
                            local_functions_tuple_neighbor,
                            local_source_entity,
                            local_source_neighbor,
                            intersection,
                            x_intersection);
  }
}; // class LocalNumericalCouplingFluxInterface


template <class Traits>
class LocalNumericalBoundaryFluxInterface
    : public XT::CRTPInterface<LocalNumericalBoundaryFluxInterface<Traits>, Traits>,
      internal::IsNumericalBoundaryFlux
{
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;


public:
  LocalNumericalBoundaryFluxInterface(const BoundaryValueType& boundary_values, const BoundaryInfoType& boundary_info)
    : boundary_values_(boundary_values)
    , boundary_info_(boundary_info)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity));
    return this->as_imp().local_functions(entity);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  auto evaluate(const LocalfunctionTupleType& local_functions_tuple,
                const XT::Functions::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_entity,
                const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& x_intersection) const ->
      typename XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple, local_source_entity, intersection, x_intersection))
    return this->as_imp().evaluate(local_functions_tuple, local_source_entity, intersection, x_intersection);
  }

protected:
  // returns u_i, u_j and x_in_inside_coords
  template <class IntersectionType, size_t boundary_value_pos>
  std::tuple<RangeType, RangeType, DomainType> get_values(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    std::tuple<RangeType, RangeType, DomainType> ret;
    auto& u_i = std::get<0>(ret);
    auto& u_j = std::get<1>(ret);
    auto& x_in_inside_coords = std::get<2>(ret);
    x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    u_i = local_source_entity.evaluate(x_in_inside_coords);
    const auto& boundary_type = boundary_info_.type(intersection);
    if (boundary_type == dirichlet_boundary_)
      u_j = std::get<boundary_value_pos>(local_functions_tuple)->evaluate(x_in_inside_coords);
    else if (boundary_type == reflecting_boundary_)
      u_j = std::get<1>(boundary_values_).evaluate(intersection, u_i);
    else if (boundary_type == absorbing_boundary_)
      u_j = u_i;
    else
      DUNE_THROW(NotImplemented, "Can't handle boundarys other than dirichlet, reflecting and absorbing!");
    return ret;
  }

  const BoundaryValueType& boundary_values_;
  const BoundaryInfoType& boundary_info_;
  static constexpr XT::Grid::DirichletBoundary dirichlet_boundary_{};
  static constexpr XT::Grid::ReflectingBoundary reflecting_boundary_{};
  static constexpr XT::Grid::AbsorbingBoundary absorbing_boundary_{};
}; // class LocalNumericalBoundaryFluxInterface

template <class Traits>
constexpr XT::Grid::DirichletBoundary LocalNumericalBoundaryFluxInterface<Traits>::dirichlet_boundary_;

template <class Traits>
constexpr XT::Grid::ReflectingBoundary LocalNumericalBoundaryFluxInterface<Traits>::reflecting_boundary_;

template <class Traits>
constexpr XT::Grid::AbsorbingBoundary LocalNumericalBoundaryFluxInterface<Traits>::absorbing_boundary_;

template <class T>
struct is_local_numerical_coupling_flux : std::is_base_of<internal::IsNumericalCouplingFlux, T>
{
};

template <class T>
struct is_local_numerical_boundary_flux : std::is_base_of<internal::IsNumericalBoundaryFlux, T>
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH
