// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_OPERATORS_DARCY_HH
#define DUNE_GDT_OPERATORS_DARCY_HH

#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/cg/interface.hh>
#include <dune/gdt/spaces/rt/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits
template <class GridViewImp, class FunctionImp>
class DarcyOperator;


namespace internal {


template <class GridViewImp, class FunctionImp>
class DarcyOperatorTraits
{
  static_assert(Stuff::is_localizable_function<FunctionImp>::value,
                "FunctionImp has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename GridViewImp::ctype, typename FunctionImp::DomainFieldType>::value,
                "Types do not match!");
  static_assert(GridViewImp::dimension == FunctionImp::dimDomain, "Dimensions do not match!");
  static_assert(FunctionImp::dimRange == FunctionImp::dimRangeCols, "Dimensions do not match!");

public:
  typedef DarcyOperator<GridViewImp, FunctionImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef typename FunctionImp::RangeFieldType FieldType;
}; // class DarcyOperatorTraits


} // namespace internal


/**
  * \note Only works for scalar valued function atm.
  * \todo add make_darcy_operator
  **/
template <class GridViewImp, class FunctionImp>
class DarcyOperator : public OperatorInterface<internal::DarcyOperatorTraits<GridViewImp, FunctionImp>>
{
public:
  typedef internal::DarcyOperatorTraits<GridViewImp, FunctionImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

  DarcyOperator(const GridViewType& grd_vw, const FunctionImp& function)
    : grid_view_(grd_vw)
    , function_(function)
  {
  }

  /**
   * \brief Applies the operator.
   * \note  See redirect_apply for the implementation (depending on the type of the range space).
   * \sa    redirect_apply
   */
  template <class S, class V, size_t r, size_t rC>
  void
  apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, r, rC>& source,
        DiscreteFunction<S, V>& range) const
  {
    redirect_apply(range.space(), source, range);
  }

private:
  template <class R, size_t r, size_t rC>
  struct Helper
  {
    static_assert(AlwaysFalse<R>::value, "This should not happen (see check in Traits)!");
  };

  template <class R, size_t r>
  struct Helper<R, r, r>
  {
    typedef Stuff::Common::FieldMatrix<R, r, r> type;
  };

  template <class R>
  struct Helper<R, 1, 1>
  {
    typedef Stuff::Common::FieldVector<R, 1> type;
  };

  typedef typename Helper<typename FunctionImp::RangeFieldType, FunctionImp::dimRange, FunctionImp::dimRangeCols>::type
      ValueType;

  /**
   * \brief Does an L2 projection of '- function * \gradient source' onto range.
   */
  template <class T, class S, class V>
  void redirect_apply(
      const CgSpaceInterface<T, dimDomain, dimDomain, 1>& /*space*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1, 1>& source,
      DiscreteFunction<S, V>& range) const
  {
    typedef typename Stuff::LA::Container<FieldType, V::sparse_matrix_type>::MatrixType MatrixType;
    MatrixType lhs(
        range.space().mapper().size(), range.space().mapper().size(), range.space().compute_volume_pattern());
    V rhs(range.space().mapper().size());

    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const auto local_function = function_.local_function(entity);
      const auto local_source   = source.local_function(entity);
      const auto basis          = range.space().base_function_set(entity);
      // do a volume quadrature
      const size_t integrand_order =
          std::max(local_function->order() + ssize_t(local_source->order()) - 1, basis.order()) + basis.order();
      const auto& quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
      const auto quadrature_it_end = quadrature.end();
      for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
        const auto xx                  = quadrature_it->position();
        const auto quadrature_weight   = quadrature_it->weight();
        const auto integration_element = entity.geometry().integrationElement(xx);
        const ValueType function_value = local_function->evaluate(xx);
        const auto source_gradient     = local_source->jacobian(xx);
        const auto basis_value         = basis.evaluate(xx);
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          const size_t global_ii = range.space().mapper().mapToGlobal(entity, ii);
          rhs.add_to_entry(global_ii,
                           integration_element * quadrature_weight * -1.0
                               * ((function_value * source_gradient[0]) * basis_value[ii]));
          for (size_t jj = 0; jj < basis.size(); ++jj) {
            const size_t global_jj = range.space().mapper().mapToGlobal(entity, jj);
            lhs.add_to_entry(
                global_ii, global_jj, integration_element * quadrature_weight * (basis_value[ii] * basis_value[jj]));
          }
        }
      } // do a volume quadrature
    } // walk the grid

    // solve
    try {
      Stuff::LA::Solver<MatrixType>(lhs).apply(rhs, range.vector());
    } catch (Stuff::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(operator_error,
                 "Application of the Darcy operator failed because a matrix could not be inverted!\n\n"
                     << "This was the original error: "
                     << ee.what());
    }
  } // ... redirect_apply(...)

  template <class T, class S, class V>
  void redirect_apply(
      const RtSpaceInterface<T, dimDomain, dimDomain, 1>& /*space*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1>& source,
      DiscreteFunction<S, V>& range) const
  {
    static_assert(RtSpaceInterface<T, dimDomain, 1>::polOrder == 0, "Untested!");
    const auto& rtn0_space = range.space();
    auto& range_vector     = range.vector();
    const auto infinity    = std::numeric_limits<FieldType>::infinity();
    for (size_t ii     = 0; ii < range_vector.size(); ++ii)
      range_vector[ii] = infinity;
    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity            = *entity_it;
      const auto local_DoF_indices  = rtn0_space.local_DoF_indices(entity);
      const auto global_DoF_indices = rtn0_space.mapper().globalIndices(entity);
      assert(global_DoF_indices.size() == local_DoF_indices.size());
      const auto local_function = function_.local_function(entity);
      const auto local_source   = source.local_function(entity);
      const auto local_basis    = rtn0_space.base_function_set(entity);
      // walk the intersections
      const auto intersection_it_end = grid_view_.iend(entity);
      for (auto intersection_it = grid_view_.ibegin(entity); intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbor_ptr = intersection.outside();
          const auto& neighbor    = *neighbor_ptr;
          if (grid_view_.indexSet().index(entity) < grid_view_.indexSet().index(neighbor)) {
            const auto local_function_neighbor    = function_.local_function(neighbor);
            const auto local_source_neighbor      = source.local_function(neighbor);
            const size_t local_intersection_index = intersection.indexInInside();
            const size_t local_DoF_index          = local_DoF_indices[local_intersection_index];
            // do a face quadrature
            FieldType lhs                = 0;
            FieldType rhs                = 0;
            const size_t integrand_order = local_function->order();
            const auto& quadrature       = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
                intersection.type(), boost::numeric_cast<int>(integrand_order));
            const auto quadrature_it_end = quadrature.end();
            for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
              const auto xx_intersection         = quadrature_it->position();
              const auto normal                  = intersection.unitOuterNormal(xx_intersection);
              const FieldType integration_factor = intersection.geometry().integrationElement(xx_intersection);
              const FieldType weight             = quadrature_it->weight();
              const auto xx_entity               = intersection.geometryInInside().global(xx_intersection);
              const auto xx_neighbor             = intersection.geometryInOutside().global(xx_intersection);
              // evaluate
              ValueType function_value = local_function->evaluate(xx_entity);
              function_value *= 0.5;
              ValueType function_value_neighbor = local_function_neighbor->evaluate(xx_neighbor);
              function_value_neighbor *= 0.5;
              function_value += function_value_neighbor;
              auto source_gradient = local_source->jacobian(xx_entity)[0];
              source_gradient *= 0.5;
              auto source_gradient_neighbor = local_source_neighbor->jacobian(xx_neighbor)[0];
              source_gradient_neighbor *= 0.5;
              source_gradient += source_gradient_neighbor;
              const auto basis_values = local_basis.evaluate(xx_entity);
              const auto basis_value  = basis_values[local_DoF_index];
              // compute integrals
              lhs += integration_factor * weight * (basis_value * normal);
              rhs += integration_factor * weight * -1.0 * compute_value(function_value, source_gradient, normal);
            } // do a face quadrature
            // set DoF
            const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
            assert(!(range_vector[global_DoF_index] < infinity));
            range_vector[global_DoF_index] = rhs / lhs;
          }
        } else if (intersection.boundary() && !intersection.neighbor()) {
          const size_t local_intersection_index = intersection.indexInInside();
          const size_t local_DoF_index          = local_DoF_indices[local_intersection_index];
          // do a face quadrature
          FieldType lhs                = 0;
          FieldType rhs                = 0;
          const size_t integrand_order = local_function->order();
          const auto& quadrature       = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
              intersection.type(), boost::numeric_cast<int>(integrand_order));
          const auto quadrature_it_end = quadrature.end();
          for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
            const auto xx_intersection    = quadrature_it->position();
            const auto normal             = intersection.unitOuterNormal(xx_intersection);
            const auto integration_factor = intersection.geometry().integrationElement(xx_intersection);
            const auto weight             = quadrature_it->weight();
            const auto xx_entity          = intersection.geometryInInside().global(xx_intersection);
            // evalaute
            const ValueType function_value = local_function->evaluate(xx_entity);
            const auto source_gradient     = local_source->jacobian(xx_entity)[0];
            const auto basis_values        = local_basis.evaluate(xx_entity);
            const auto basis_value         = basis_values[local_DoF_index];
            // compute integrals
            lhs += integration_factor * weight * (basis_value * normal);
            rhs += integration_factor * weight * -1.0 * compute_value(function_value, source_gradient, normal);
          } // do a face quadrature
          // set DoF
          const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
          assert(!(range_vector[global_DoF_index] < infinity));
          range_vector[global_DoF_index] = rhs / lhs;
        } else
          DUNE_THROW(Stuff::Exceptions::internal_error, "Unknown intersection type!");
      } // walk the intersections
    } // walk the grid
  } // ... redirect_apply(...)

  template <class R, int d>
  R compute_value(const FieldVector<R, 1>& function_value, const FieldVector<R, d>& source_gradient,
                  const FieldVector<R, d>& normal) const
  {
    return function_value * (source_gradient * normal);
  }

  template <class R, int d>
  typename std::enable_if<(d > 1), R>::type compute_value(const Stuff::Common::FieldMatrix<R, d, d>& function_value,
                                                          const FieldVector<R, d>& source_gradient,
                                                          const FieldVector<R, d>& normal) const
  {
    return (function_value * source_gradient) * normal;
  }

  const GridViewType& grid_view_;
  const FunctionImp& function_;
}; // class DarcyOperator


template <class G, class F>
std::unique_ptr<DarcyOperator<G, F>>
make_darcy(const G& grid_view, const F& function)
{
  return std::unique_ptr<DarcyOperator<G, F>>(new DarcyOperator<G, F>(grid_view, function));
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_DARCY_HH
