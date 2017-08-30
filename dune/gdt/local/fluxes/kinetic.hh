// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FLUXES_KINETIC_HH
#define DUNE_GDT_LOCAL_FLUXES_KINETIC_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/la/container/eigen.hh>

#include "interfaces.hh"
#include "entropybased.hh"
#include "godunov.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class BasisfunctionImp>
class KineticLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueType, class BasisfunctionImp>
class KineticLocalNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp, class BasisfunctionImp>
class KineticLocalNumericalCouplingFluxTraits : public GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp>
{
public:
  typedef BasisfunctionImp BasisfunctionType;
  typedef std::tuple<> LocalfunctionTupleType;
  typedef KineticLocalNumericalCouplingFlux<AnalyticalFluxImp, BasisfunctionImp> derived_type;
}; // class KineticLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp>
class KineticLocalNumericalBoundaryFluxTraits
    : public KineticLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, BasisfunctionImp>
{
  typedef KineticLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, BasisfunctionImp> BaseType;

public:
  typedef BoundaryValueImp BoundaryValueType;
  typedef typename BoundaryValueType::LocalfunctionType BoundaryValueLocalfunctionType;
  typedef KineticLocalNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, BasisfunctionImp> derived_type;
  typedef std::tuple<std::shared_ptr<BoundaryValueLocalfunctionType>> LocalfunctionTupleType;
}; // class KineticLocalNumericalBoundaryFluxTraits


template <class AnalyticalFluxType>
struct is_entropy_flux_candidate
{
  DXTC_has_typedef_initialize_once(BasisfunctionType);
  DXTC_has_typedef_initialize_once(GridLayerType);
  DXTC_has_typedef_initialize_once(StateType);
  DXTC_has_static_member_initialize_once(dimQuadrature);

  static const bool value = DXTC_has_typedef(BasisfunctionType)<AnalyticalFluxType>::value
                            && DXTC_has_typedef(GridLayerType)<AnalyticalFluxType>::value
                            && DXTC_has_typedef(StateType)<AnalyticalFluxType>::value
                            && DXTC_has_static_member(dimQuadrature)<AnalyticalFluxType>::value;
}; // class is_entropy_flux_candidate

template <class AnalyticalFluxType, bool candidate = is_entropy_flux_candidate<AnalyticalFluxType>::value>
struct is_entropy_flux : public std::is_base_of<EntropyBasedLocalFlux<typename AnalyticalFluxType::BasisfunctionType,
                                                                      typename AnalyticalFluxType::GridLayerType,
                                                                      typename AnalyticalFluxType::StateType,
                                                                      AnalyticalFluxType::dimQuadrature>,
                                                AnalyticalFluxType>
{
};

template <class AnalyticalFluxType>
struct is_entropy_flux<AnalyticalFluxType, false> : public std::false_type
{
};


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
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename XT::LA::CommonSparseMatrix<RangeFieldType> SparseMatrixType;

  explicit KineticFluxImplementation(const AnalyticalFluxType& analytical_flux,
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
    size_t direction = intersection.indexInInside() / 2;
    auto n_ij = intersection.unitOuterNormal(x_in_intersection_coords);
    return helper<>::evaluate_kinetic_flux(analytical_flux_,
                                           basis_functions_,
                                           entity,
                                           x_in_inside_coords,
                                           u_i,
                                           neighbor,
                                           x_in_outside_coords,
                                           u_j,
                                           n_ij,
                                           direction,
                                           param_inside_,
                                           param_outside_);
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
    // transpose
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = ii + 1; jj < dimRange; ++jj) {
        RangeFieldType tmp_entry = mass_matrix[ii][jj];
        mass_matrix[ii][jj] = mass_matrix[jj][ii];
        mass_matrix[jj][ii] = tmp_entry;
      }
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

  template <bool entropy_flux = is_entropy_flux<AnalyticalFluxType>::value, class anything = void>
  struct helper
  {
    static RangeType evaluate_kinetic_flux(const AnalyticalFluxType& analytical_flux,
                                           const BasisfunctionType& /*basis_functions*/,
                                           const EntityType& entity,
                                           const DomainType& x_local_entity,
                                           const RangeType& u_i,
                                           const EntityType& neighbor,
                                           const DomainType& x_local_neighbor,
                                           const RangeType u_j,
                                           const DomainType& n_ij,
                                           size_t direction,
                                           const XT::Common::Parameter& param,
                                           const XT::Common::Parameter& param_neighbor)
    {
      static_assert(std::is_same<BasisfunctionType, typename AnalyticalFluxType::BasisfunctionType>::value, "");
      return analytical_flux.evaluate_kinetic_flux(
          entity, x_local_entity, u_i, neighbor, x_local_neighbor, u_j, n_ij, direction, param, param_neighbor);
    }
  }; // struct helper<true, ...>

  template <class anything>
  struct helper<false, anything>
  {
    static RangeType evaluate_kinetic_flux(const AnalyticalFluxType& /*analytical_flux*/,
                                           const BasisfunctionType& basis_functions,
                                           const EntityType& /*entity*/,
                                           const DomainType& /*x_local_entity*/,
                                           const RangeType& u_i,
                                           const EntityType& /*neighbor*/,
                                           const DomainType& /*x_local_neighbor*/,
                                           const RangeType u_j,
                                           const DomainType& n_ij,
                                           const size_t direction,
                                           const XT::Common::Parameter& /*param*/,
                                           const XT::Common::Parameter& /*param_neighbor*/)
    {
      static auto flux_matrices = initialize_flux_matrices(basis_functions);
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
  }; // struct helper<false, ...>

  const AnalyticalFluxType& analytical_flux_;
  const BasisfunctionType& basis_functions_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;
}; // class KineticFluxImplementation<...>


} // namespace internal


/**
 *  \brief  Kinetic flux evaluation for inner intersections and periodic boundary intersections.
 */
template <class AnalyticalFluxImp, class BasisfunctionImp>
class KineticLocalNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::KineticLocalNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                   BasisfunctionImp>>
{
public:
  typedef internal::KineticLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, BasisfunctionImp> Traits;
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
                                             const XT::Common::Parameter param,
                                             const BasisfunctionType& basis_functions)
    : implementation_(analytical_flux, param, basis_functions, false)
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
*  \brief  Kinetic flux evaluation for Dirichlet boundary intersections.
*/
template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp>
class KineticLocalNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::KineticLocalNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                   BoundaryValueImp,
                                                                                                   BasisfunctionImp>>
{
public:
  typedef internal::KineticLocalNumericalBoundaryFluxTraits<AnalyticalFluxImp, BoundaryValueImp, BasisfunctionImp>
      Traits;
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
                                             const XT::Common::Parameter param,
                                             const BasisfunctionType& basis_functions)
    : implementation_(analytical_flux, param, basis_functions, true)
    , boundary_values_(boundary_values)
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
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    // get function values
    const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
    const RangeType u_i = local_source_entity.evaluate(x_intersection_entity_coords);
    auto u_j = std::get<0>(local_functions_tuple)->evaluate(x_intersection_entity_coords);
    return implementation_.evaluate(intersection,
                                    intersection.inside(),
                                    intersection.inside(),
                                    x_intersection,
                                    x_intersection_entity_coords,
                                    x_intersection_entity_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const internal::KineticFluxImplementation<Traits> implementation_;
  const BoundaryValueType& boundary_values_;
}; // class KineticLocalNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_KINETIC_HH
