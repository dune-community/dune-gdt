// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_KINETIC_HH
#define DUNE_GDT_LOCAL_FLUXES_KINETIC_HH

#include <tuple>
#include <memory>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/functions/interfaces.hh>

#include "base.hh"
#include "interfaces.hh"
#include "entropybased.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class BasisfunctionImp, class GridLayerImp, class Traits>
class KineticLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp, class GridLayerImp, class Traits>
class KineticLocalNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class BasisfunctionImp, class GridLayerImp>
class KineticLocalNumericalCouplingFluxTraits : public NumericalCouplingFluxTraitsBase<AnalyticalFluxImp>
{
public:
  typedef BasisfunctionImp BasisfunctionType;
  typedef GridLayerImp GridLayerType;
  typedef std::tuple<> LocalfunctionTupleType;
  typedef KineticLocalNumericalCouplingFlux<AnalyticalFluxImp,
                                            BasisfunctionImp,
                                            GridLayerImp,
                                            KineticLocalNumericalCouplingFluxTraits>
      derived_type;
}; // class KineticLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp, class GridLayerImp>
class KineticLocalNumericalBoundaryFluxTraits
    : public NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  typedef NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  using typename BaseType::LocalBoundaryValueType;
  typedef BasisfunctionImp BasisfunctionType;
  typedef GridLayerImp GridLayerType;
  typedef std::tuple<std::unique_ptr<LocalBoundaryValueType>> LocalfunctionTupleType;
  typedef KineticLocalNumericalBoundaryFlux<AnalyticalFluxImp,
                                            BoundaryValueImp,
                                            BasisfunctionImp,
                                            GridLayerImp,
                                            KineticLocalNumericalBoundaryFluxTraits>
      derived_type;
}; // class KineticLocalNumericalBoundaryFluxTraits


template <class Traits>
class KineticFluxImplementation
{
public:
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::AnalyticalFluxLocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename Traits::BasisfunctionType BasisfunctionType;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename AnalyticalFluxType::StateType StateType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename XT::LA::CommonSparseOrDenseMatrixCsr<RangeFieldType> SparseMatrixType;
  typedef EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, StateType> EntropyFluxType;

  explicit KineticFluxImplementation(const AnalyticalFluxType* analytical_flux,
                                     const XT::Common::Parameter& param,
                                     const BasisfunctionType& basis_functions,
                                     const bool boundary)
    : analytical_flux_(analytical_flux)
    , param_inside_(param)
    , param_outside_(param)
    , basis_functions_(basis_functions)
  {
    param_inside_.set("boundary", {0.}, true);
    param_outside_.set("boundary", {double(boundary)}, true);
  }

  template <class IntersectionType>
  RangeType evaluate(const IntersectionType& intersection,
                     const EntityType& entity,
                     const EntityType& neighbor,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords,
                     const DomainType& x_in_inside_coords,
                     const DomainType& x_in_outside_coords,
                     const RangeType& u_i,
                     const RangeType& u_j) const
  {
    // find direction of unit outer normal
    assert(intersection.indexInInside() >= 0);
    size_t direction = static_cast<size_t>(intersection.indexInInside()) / 2;
    auto n_ij = intersection.unitOuterNormal(x_in_intersection_coords);
    if (dynamic_cast<const EntropyFluxType*>(analytical_flux_) != nullptr) {
      return dynamic_cast<const EntropyFluxType*>(analytical_flux_)
          ->evaluate_kinetic_flux(entity,
                                  x_in_inside_coords,
                                  u_i,
                                  neighbor,
                                  x_in_outside_coords,
                                  u_j,
                                  n_ij,
                                  direction,
                                  param_inside_,
                                  param_outside_);
    } else {
      static const auto flux_matrices = initialize_flux_matrices(basis_functions_);
      RangeType ret(0);
      auto tmp_vec = ret;
      const auto& inner_flux_matrix = flux_matrices[direction][n_ij[direction] > 0 ? 1 : 0];
      const auto& outer_flux_matrix = flux_matrices[direction][n_ij[direction] > 0 ? 0 : 1];
      inner_flux_matrix.mv(u_i, tmp_vec);
      outer_flux_matrix.mv(u_j, ret);
      ret += tmp_vec;
      ret *= n_ij[direction];
      return ret;
    }
  } // ... evaluate(...)

  const AnalyticalFluxType& analytical_flux() const
  {
    return analytical_flux_;
  }

private:
  static FieldVector<FieldVector<SparseMatrixType, 2>, dimDomain>
  initialize_flux_matrices(const BasisfunctionType& basis_functions)
  {
    // calculate < v_i b b^T >_- M^{-1} and < v_i b b^T >_+ M^{-1}
    auto kinetic_flux_matrices = basis_functions.kinetic_flux_matrices();
    auto mass_matrix = basis_functions.mass_matrix();
    auto flux_matrices_dense = kinetic_flux_matrices;
    FieldVector<FieldVector<SparseMatrixType, 2>, dimDomain> flux_matrices(
        FieldVector<SparseMatrixType, 2>(SparseMatrixType(dimRange, dimRange, size_t(0))));
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t kk = 0; kk < 2; ++kk) {
        for (size_t rr = 0; rr < dimRange; ++rr)
          mass_matrix.solve(flux_matrices_dense[dd][kk][rr], kinetic_flux_matrices[dd][kk][rr]);
        flux_matrices[dd][kk] = flux_matrices_dense[dd][kk];
      } // kk
    } // dd
    return flux_matrices;
  }

  const AnalyticalFluxType* analytical_flux_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;
  const BasisfunctionType& basis_functions_;
}; // class KineticFluxImplementation<...>


} // namespace internal


/**
 *  \brief  Kinetic flux evaluation for inner intersections and periodic boundary intersections.
 */
template <class AnalyticalFluxImp,
          class BasisfunctionImp,
          class GridLayerImp,
          class Traits =
              internal::KineticLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, BasisfunctionImp, GridLayerImp>>
class KineticLocalNumericalCouplingFlux : public LocalNumericalCouplingFluxInterface<Traits>
{
public:
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::BasisfunctionType BasisfunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit KineticLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                             const XT::Common::Parameter& param,
                                             const BasisfunctionType& basis_functions)
    : implementation_(&analytical_flux, param, basis_functions, false)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return std::make_tuple();
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& /*local_functions_tuple_entity*/,
      const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_neighbor,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    // get function values
    const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
    const auto x_intersection_neighbor_coords = intersection.geometryInOutside().global(x_intersection);
    const RangeType u_i = local_source_entity.evaluate(x_intersection_entity_coords);
    RangeType u_j = local_source_neighbor.evaluate(x_intersection_neighbor_coords);
    return implementation_.evaluate(intersection,
                                    intersection.inside(),
                                    intersection.outside(),
                                    x_intersection,
                                    x_intersection_entity_coords,
                                    x_intersection_neighbor_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const internal::KineticFluxImplementation<Traits> implementation_;
}; // class KineticLocalNumericalCouplingFlux

/**
*  \brief  Kinetic flux evaluation for boundary intersections.
*/
template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BasisfunctionImp,
          class GridLayerImp,
          class Traits = internal::KineticLocalNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                           BoundaryValueImp,
                                                                           BasisfunctionImp,
                                                                           GridLayerImp>>
class KineticLocalNumericalBoundaryFlux : public LocalNumericalBoundaryFluxInterface<Traits>
{
  typedef LocalNumericalBoundaryFluxInterface<Traits> InterfaceType;

public:
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::BasisfunctionType BasisfunctionType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit KineticLocalNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const XT::Common::Parameter& param,
                                             const BasisfunctionType& basis_functions)
    : InterfaceType(boundary_values)
    , implementation_(&analytical_flux, param, basis_functions, true)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(boundary_values_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    const auto values = InterfaceType::template get_values<IntersectionType, 0>(
        local_functions_tuple, local_source_entity, intersection, x_in_intersection_coords);
    return implementation_.evaluate(intersection,
                                    intersection.inside(),
                                    intersection.inside(),
                                    x_in_intersection_coords,
                                    std::get<2>(values),
                                    std::get<2>(values),
                                    std::get<0>(values),
                                    std::get<1>(values));
  } // RangeType evaluate(...) const

private:
  const internal::KineticFluxImplementation<Traits> implementation_;
  using InterfaceType::boundary_values_;
}; // class KineticLocalNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_KINETIC_HH
