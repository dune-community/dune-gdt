// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_DARCY_HH
#define DUNE_GDT_OPERATORS_DARCY_HH

#include <limits>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/spaces/cg/fem.hh>
#include <dune/gdt/playground/spaces/rt/pdelab.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forward, to be used in the traits
template <class GridViewImp, class FunctionImp>
class Darcy;


template <class GridViewImp, class FunctionImp>
class DarcyTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, FunctionImp>::value,
                "FunctionImp has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename GridViewImp::ctype, typename FunctionImp::DomainFieldType>::value,
                "Types do not match!");
  static_assert(GridViewImp::dimension == FunctionImp::dimDomain, "Dimensions do not match!");

public:
  typedef Darcy<GridViewImp, FunctionImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef typename FunctionImp::RangeFieldType FieldType;
}; // class DarcyTraits


/**
  * \note Only works for scalar valued function atm.
  **/
template <class GridViewImp, class FunctionImp>
class Darcy : public OperatorInterface<DarcyTraits<GridViewImp, FunctionImp>>
{
public:
  typedef DarcyTraits<GridViewImp, FunctionImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  Darcy(const GridViewType& grd_vw, const FunctionImp& function)
    : grid_view_(grd_vw)
    , function_(function)
  {
  }

  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void apply(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
             DiscreteFunction<SpaceInterface<T, d, r, rC>, V>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<E>::value, "Not implemented for this combination of source and range!");
  }

  /**
   * \brief Does an L2 projection of '- function * \gradient source' onto range.
   */
  template <class GP, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1, 1>& source,
             DiscreteFunction<Spaces::CG::FemBased<GP, p, FieldType, dimDomain, 1>, V>& range) const
  {
    typedef typename Stuff::LA::Container<FieldType>::MatrixType MatrixType;
    typedef typename Stuff::LA::Container<FieldType>::VectorType VectorType;
    MatrixType lhs(
        range.space().mapper().size(), range.space().mapper().size(), range.space().compute_volume_pattern());
    VectorType rhs(range.space().mapper().size());

    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const auto local_function = function_.local_function(entity);
      const auto local_source   = source.local_function(entity);
      const auto basis          = range.space().base_function_set(entity);
      // do a volume quadrature
      const size_t integrand_order =
          std::max(local_function->order() + size_t(local_source->order() - 1), basis.order()) + basis.order();
      assert(integrand_order < std::numeric_limits<int>::max());
      const auto& quadrature       = QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), int(integrand_order));
      const auto quadrature_it_end = quadrature.end();
      for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
        const auto xx                             = quadrature_it->position();
        const FieldType quadrature_weight         = quadrature_it->weight();
        const DomainFieldType integration_element = entity.geometry().integrationElement(xx);
        const auto function_value                 = local_function->evaluate(xx);
        const auto source_gradient                = local_source->jacobian(xx);
        const auto basis_value = basis.evaluate(xx);
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          const size_t global_ii = range.space().mapper().mapToGlobal(entity, ii);
          rhs.add_to_entry(global_ii,
                           integration_element * quadrature_weight * -1.0 * function_value
                               * (source_gradient[0] * basis_value[ii]));
          for (size_t jj = 0; jj < basis.size(); ++jj) {
            const size_t global_jj = range.space().mapper().mapToGlobal(entity, jj);
            lhs.add_to_entry(
                global_ii, global_jj, integration_element * quadrature_weight * (basis_value[ii] * basis_value[jj]));
          }
        }
      } // do a volume quadrature
    } // walk the grid

    // solve
    Stuff::LA::Solver<MatrixType>(lhs).apply(rhs, range.vector());
  } // ... apply(...)

  template <class GP, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1>& source,
             DiscreteFunction<Spaces::RT::PdelabBased<GP, 0, FieldType, dimDomain>, V>& range) const
  {
    const auto& rtn0_space   = range.space();
    auto& range_vector       = range.vector();
    const FieldType infinity = std::numeric_limits<FieldType>::infinity();
    for (size_t ii = 0; ii < range_vector.size(); ++ii)
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
          const auto& neighbor = *neighbor_ptr;
          if (grid_view_.indexSet().index(entity) < grid_view_.indexSet().index(neighbor)) {
            const auto local_function_neighbor    = function_.local_function(neighbor);
            const auto local_source_neighbor      = source.local_function(neighbor);
            const size_t local_intersection_index = intersection.indexInInside();
            const size_t local_DoF_index          = local_DoF_indices[local_intersection_index];
            // do a face quadrature
            FieldType lhs                = 0;
            FieldType rhs                = 0;
            const size_t integrand_order = local_function->order();
            assert(integrand_order < std::numeric_limits<int>::max());
            const auto& quadrature =
                QuadratureRules<DomainFieldType, dimDomain - 1>::rule(intersection.type(), int(integrand_order));
            const auto quadrature_it_end = quadrature.end();
            for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
              const auto xx_intersection         = quadrature_it->position();
              const auto normal                  = intersection.unitOuterNormal(xx_intersection);
              const FieldType integration_factor = intersection.geometry().integrationElement(xx_intersection);
              const FieldType weigth             = quadrature_it->weight();
              const auto xx_entity               = intersection.geometryInInside().global(xx_intersection);
              const auto xx_neighbor             = intersection.geometryInOutside().global(xx_intersection);
              // evalaute
              auto function_value = local_function->evaluate(xx_entity);
              function_value *= 0.5;
              auto function_value_neighbor = local_function_neighbor->evaluate(xx_neighbor);
              function_value_neighbor *= 0.5;
              function_value += function_value_neighbor;
              auto source_gradient = local_source->jacobian(xx_entity);
              source_gradient *= 0.5;
              auto source_gradient_neighbor = local_source_neighbor->jacobian(xx_neighbor);
              source_gradient_neighbor *= 0.5;
              source_gradient += source_gradient_neighbor;
              const auto basis_values = local_basis.evaluate(xx_entity);
              const auto basis_value  = basis_values[local_DoF_index];
              // compute integrals
              lhs += integration_factor * weigth * (basis_value * normal);
              rhs += integration_factor * weigth * -1.0 * function_value * (source_gradient[0] * normal);
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
          assert(integrand_order < std::numeric_limits<int>::max());
          const auto& quadrature =
              QuadratureRules<DomainFieldType, dimDomain - 1>::rule(intersection.type(), int(integrand_order));
          const auto quadrature_it_end = quadrature.end();
          for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
            const auto xx_intersection         = quadrature_it->position();
            const auto normal                  = intersection.unitOuterNormal(xx_intersection);
            const FieldType integration_factor = intersection.geometry().integrationElement(xx_intersection);
            const FieldType weigth             = quadrature_it->weight();
            const auto xx_entity               = intersection.geometryInInside().global(xx_intersection);
            // evalaute
            const auto function_value  = local_function->evaluate(xx_entity);
            const auto source_gradient = local_source->jacobian(xx_entity);
            const auto basis_values    = local_basis.evaluate(xx_entity);
            const auto basis_value     = basis_values[local_DoF_index];
            // compute integrals
            lhs += integration_factor * weigth * (basis_value * normal);
            rhs += integration_factor * weigth * -1.0 * function_value * (source_gradient[0] * normal);
          } // do a face quadrature
          // set DoF
          const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
          assert(!(range_vector[global_DoF_index] < infinity));
          range_vector[global_DoF_index] = rhs / lhs;
        } else
          DUNE_THROW(Stuff::Exceptions::internal_error, "Unknown intersection type!");
      } // walk the intersections
    } // walk the grid
  } // ... apply(...)

private:
  const GridViewType& grid_view_;
  const FunctionImp& function_;
}; // class Darcy


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_DARCY_HH
