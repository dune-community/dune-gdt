// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_RECONSTRUCTIONS_HH
#define DUNE_GDT_OPERATOR_RECONSTRUCTIONS_HH

#include <type_traits>
#include <limits>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/raviartthomas/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localevaluation/swipdg.hh>

namespace Dune {
namespace GDT {
namespace ReconstructionOperator {


template <class GridPartType, class LocalizableFunctionType>
class DiffusiveFlux
{
  static_assert(GridPartType::dimension == 2, "Only implemented for dimDomain 2 at the moment!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableFunctionType>::value,
                "LocalizableFunctionType has to be tagged as Stuff::IsLocalizableFunction!");

public:
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  DiffusiveFlux(const GridPartType& grid_part, const LocalizableFunctionType& diffusion)
    : grid_part_(grid_part)
    , diffusion_(diffusion)
  {
  }

  template <class SGP, class R, class SV, class RGP, class RV>
  void apply(const ConstDiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<SGP, 1, R, 1>, SV>& source,
             DiscreteFunction<RaviartThomasSpace::FemLocalfunctionsWrapper<RGP, 0, R, dimDomain>, RV>& range) const
  {
    typedef ConstDiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<SGP, 1, R, 1>, SV> SourceType;
    typedef DiscreteFunction<RaviartThomasSpace::FemLocalfunctionsWrapper<RGP, 0, R, dimDomain>, RV> RangeType;
    static_assert(dimDomain == 2, "Only implemented for dimDomain 2!");
    // prepare tmp storage
    std::vector<size_t> intersection_id_to_local_DoF_id_map(
        3, 0); // <- there should be 3 intersections on a simplex in 2d!
    const auto& rtn_space = range.space();
    std::vector<typename RangeType::RangeType> rtn_basis_values(rtn_space.mapper().maxNumDofs(),
                                                                typename RangeType::RangeType(0));
    typename SourceType::DomainType unit_outer_normal(0);
    typename SourceType::DomainType point_entity(0);
    Dune::DynamicMatrix<R> tmp_result(1, 1, 0);
    Dune::DynamicMatrix<R> tmp_result_en_en(1, 1, 0);
    Dune::DynamicMatrix<R> tmp_result_en_ne(1, 1, 0);
    // create local evaluations
    const LocalEvaluation::SWIPDG::BoundaryLHS<LocalizableFunctionType> boundary_evaluation(diffusion_);
    const LocalEvaluation::SWIPDG::Inner<LocalizableFunctionType> inner_evaluation(diffusion_);
    const Stuff::Function::Constant<EntityType, DomainFieldType, dimDomain, R, 1> constant_one(1);
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_source       = source.local_discrete_function(entity);
      auto local_range              = range.local_discrete_function(entity);
      auto local_range_vector       = local_range.vector();
      const auto local_diffusion    = diffusion_.local_function(entity);
      const auto local_constant_one = constant_one.local_function(entity);
      const auto rtn_basis          = rtn_space.baseFunctionSet(entity);
      // get the local finite element
      const auto rtn_finite_element      = rtn_space.backend().finiteElement(entity);
      const auto& rtn_local_coefficients = rtn_finite_element.localCoefficients();
      assert(intersection_id_to_local_DoF_id_map.size() >= rtn_local_coefficients.size());
      assert(rtn_local_coefficients.size() < std::numeric_limits<int>::max());
      for (size_t ii = 0; ii < rtn_local_coefficients.size(); ++ii)
        intersection_id_to_local_DoF_id_map[rtn_local_coefficients.localKey(int(ii)).subEntity()] = ii;
      // loop over all intersections
      const auto intersection_end_it = grid_part_.iend(entity);
      for (auto intersection_it = grid_part_.ibegin(entity); intersection_it != intersection_end_it;
           ++intersection_it) {
        // get intersection information
        const auto& intersection        = *intersection_it;
        const size_t intersection_index = intersection.indexInInside();
        assert(intersection_index < intersection_id_to_local_DoF_id_map.size() && "This should not happen!");
        const size_t intersection_DoF_index = intersection_id_to_local_DoF_id_map[intersection_index];
        // prepare quadrature
        const size_t quadrature_order =
            std::max(rtn_basis.order(),
                     std::max(inner_evaluation.order(*local_diffusion,
                                                     *local_diffusion,
                                                     local_source,
                                                     *local_constant_one,
                                                     local_source,
                                                     *local_constant_one),
                              boundary_evaluation.order(*local_diffusion, local_source, *local_constant_one)));
        assert((2 * quadrature_order + 1) < std::numeric_limits<int>::max());
        const auto& face_quadrature =
            QuadratureRules<DomainFieldType, dimDomain - 1>::rule(intersection.type(), int(2 * quadrature_order + 1));
        R lhs_integral = 0;
        R rhs_integral = 0;
        if (intersection.boundary() && !intersection.neighbor()) {
          // we are on the domain boundary
          // loop over all quadrature points
          for (auto quadrature_point : face_quadrature) {
            const FieldVector<DomainFieldType, dimDomain - 1>& point_intersection = quadrature_point.position();
            point_entity                                                          = intersection.geometryInInside().global(point_intersection);
            const R integration_factor                                            = intersection.geometry().integrationElement(point_intersection);
            const R quadrature_weight                                             = quadrature_point.weight();
            // evaluate
            unit_outer_normal = intersection.unitOuterNormal(point_intersection);
            rtn_basis.evaluate(point_entity, rtn_basis_values);
            assert(rtn_basis_values.size() > intersection_DoF_index);
            const auto& rtn_basis_value = rtn_basis_values[intersection_DoF_index];
            // compute integrals
            lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
            tmp_result *= R(0);
            boundary_evaluation.evaluate(
                *local_diffusion, *local_constant_one, local_source, intersection, point_intersection, tmp_result);
            assert(tmp_result.rows() == 1);
            assert(tmp_result.cols() == 1);
            rhs_integral += integration_factor * quadrature_weight * tmp_result[0][0];
          } // loop over all quadrature points
          // set local DoF
          local_range_vector.set(intersection_DoF_index, rhs_integral / lhs_integral);
        } else if (intersection.neighbor() && !intersection.boundary()) {
          // we are on the inside
          // get the neighbour
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour    = *neighbour_ptr;
          // work only once on each intersection
          if (grid_part_.indexSet().index(entity) < grid_part_.indexSet().index(neighbour)) {
            // get the neighbours local functions
            const auto local_source_neighbour       = source.local_discrete_function(neighbour);
            const auto local_diffusion_neighbour    = diffusion_.local_function(neighbour);
            const auto local_constant_one_neighbour = constant_one.local_function(neighbour);
            // loop over all quadrature points
            for (auto quadrature_point : face_quadrature) {
              const FieldVector<DomainFieldType, dimDomain - 1>& point_intersection = quadrature_point.position();
              point_entity                                                          = intersection.geometryInInside().global(point_intersection);
              const R integration_factor                                            = intersection.geometry().integrationElement(point_intersection);
              const R quadrature_weight                                             = quadrature_point.weight();
              // evaluate
              unit_outer_normal = intersection.unitOuterNormal(point_intersection);
              rtn_basis.evaluate(point_entity, rtn_basis_values);
              const auto& rtn_basis_value = rtn_basis_values[intersection_DoF_index];
              // compute integrals
              lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
              tmp_result *= R(0);
              tmp_result_en_en *= R(0);
              tmp_result_en_ne *= R(0);
              inner_evaluation.evaluate(*local_diffusion,
                                        *local_diffusion_neighbour,
                                        *local_constant_one,
                                        local_source,
                                        *local_constant_one_neighbour,
                                        local_source_neighbour,
                                        intersection,
                                        point_intersection,
                                        tmp_result_en_en, // <- we are interested in this one
                                        tmp_result,
                                        tmp_result_en_ne, // <- and this one
                                        tmp_result);
              assert(tmp_result_en_en.rows() == 1);
              assert(tmp_result_en_en.cols() == 1);
              assert(tmp_result_en_ne.rows() == 1);
              assert(tmp_result_en_ne.cols() == 1);
              rhs_integral +=
                  integration_factor * quadrature_weight * (tmp_result_en_en[0][0] + tmp_result_en_ne[0][0]);
            } // loop over all quadrature points
            // set local DoF
            local_range_vector.set(intersection_DoF_index, rhs_integral / lhs_integral);
          } // work only once on each intersection
        } else
          DUNE_THROW(InvalidStateException, "Unknwon intersection type encountered!");
      } // loop over all intersections
    } // walk the grid
  } // ... apply(...)

private:
  const GridPartType& grid_part_;
  const LocalizableFunctionType& diffusion_;
}; // class DiffusiveFlux


} // namespace ReconstructionOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_RECONSTRUCTIONS_HH
