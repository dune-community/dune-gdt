// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PROLONGATIONS_HH
#define DUNE_GDT_OPERATOR_PROLONGATIONS_HH

#include <type_traits>
#include <vector>
#include <limits>

#include <dune/common/dynmatrix.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/grid/search.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>

namespace Dune {
namespace GDT {
namespace ProlongationOperator {


class Generic
{
public:
  template< class SGP, int sp, class R, int dimRange, class SV, class RGP, int rp, class RV >
  static void apply(const DiscreteFunctionDefaultConst< ContinuousLagrangeSpace::FemWrapper< SGP, sp, R, dimRange >, SV >& source,
                    DiscreteFunctionDefault< ContinuousLagrangeSpace::FemWrapper< RGP, rp, R, dimRange >, RV >& range)
  {
    //checks
    if (source.space().mapper().size() > range.space().mapper().size())
      DUNE_THROW(RangeError, "A prolongation from a larger to a smaller space does not make sense!");

    // create search in the source grid part
    typedef DiscreteFunctionDefaultConst< ContinuousLagrangeSpace::FemWrapper< SGP, sp, R, dimRange >, SV > SourceType;
    typedef DiscreteFunctionDefault< ContinuousLagrangeSpace::FemWrapper< RGP, rp, R, dimRange >, RV >      TargetType;

    typedef typename SourceType::SpaceType::GridPartType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch< SourceGridViewType > EntitySearch;
    EntitySearch entity_search(source.space().gridPart().gridView());

    // set all dofs to infinity
    const auto infinity = std::numeric_limits< R >::infinity();
    for (size_t ii = 0; ii < range.space().mapper().size(); ++ii)
      range.vector()->backend()[ii] = infinity;

    // walk the grid
    const auto entity_it_end = range.space().gridPart().template end< 0 >();
    for (auto entity_it = range.space().gridPart().template begin< 0 >();
         entity_it != entity_it_end;
         ++entity_it) {
      const auto& entity = *entity_it;

      // get global lagrange point coordinates
      const auto lagrange_point_set = range.space().backend().lagrangePointSet(entity);
      typedef FieldVector< typename SourceGridViewType::ctype, SourceGridViewType::dimension > DomainType;
      std::vector< DomainType > lagrange_points(lagrange_point_set.nop());
      for (size_t ii = 0; ii < lagrange_point_set.nop(); ++ii)
        lagrange_points[ii] = entity.geometry().global(lagrange_point_set.point(ii));

      // get source entities
      const auto source_entities = entity_search(lagrange_points);
      assert(source_entities.size() == lagrange_points.size());

      auto local_range = range.localFunction(entity);
      auto local_range_vector = local_range.vector();

      size_t kk = 0;
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        if (std::isinf(local_range_vector.get(kk))) {
          const auto& global_point = lagrange_points[ii];
          // evaluate source function
          const auto source_entity_ptr = source_entities[ii];
          const auto& source_entity = *source_entity_ptr;
          const auto local_source_point = source_entity.geometry().local(global_point);
          const auto local_source = source.localFunction(source_entity);
          const auto source_value = local_source.evaluate(local_source_point);
          for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
            local_range_vector.set(kk, source_value[jj]);
        }
        else
          kk += dimRange;
      }
    } // walk the grid
  } // ... apply(...)

  template< class SGP, int sp, class R, int dimRange, class SV, class RGP, int rp, class RV >
  static void apply(const DiscreteFunctionDefaultConst< DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< SGP, sp, R, dimRange >, SV >& source,
                    DiscreteFunctionDefault< DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< RGP, rp, R, dimRange >, RV >& range)
  {
    //checks
    static_assert(rp >= sp, "A prolongation from a higher to a lower polynomial order does not make sense!");
    if (source.space().mapper().size() > range.space().mapper().size())
      DUNE_THROW(RangeError, "A prolongation from a larger to a smaller space does not make sense!");

    // create search in the source grid part
    typedef DiscreteFunctionDefaultConst< DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< SGP, sp, R, dimRange >, SV > SourceType;
    typedef DiscreteFunctionDefault< DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< RGP, rp, R, dimRange >, RV >      RangeType;

    typedef typename SourceType::SpaceType::GridPartType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch< SourceGridViewType > EntitySearch;
    EntitySearch entity_search(source.space().gridPart().gridView());

    // walk the grid
    const auto entity_it_end = range.space().gridPart().template end< 0 >();
    for (auto entity_it = range.space().gridPart().template begin< 0 >();
         entity_it != entity_it_end;
         ++entity_it) {
      // prepare
      const auto& entity = *entity_it;
      const auto local_basis = range.space().baseFunctionSet(entity);
      DynamicMatrix< R > local_matrix(local_basis.size(), local_basis.size(), R(0));
      DynamicVector< R > local_vector(local_basis.size(), R(0));
      std::vector< typename RangeType::SpaceType::BaseFunctionSetType::RangeType > basis_values(local_basis.size());

      // create quadrature
      const size_t quadrature_order = 2*rp;
      typedef typename SourceGridViewType::ctype DomainFieldType;
      const unsigned int dimDomain = SourceGridViewType::dimension;
      const auto& quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                   2*quadrature_order + 1);

      // get global quadrature coordinates
      typedef FieldVector< DomainFieldType, dimDomain > DomainType;
      std::vector< DomainType > global_quadrature_points(quadrature.size());
      size_t pp = 0;
      for (const auto& quadrature_point : quadrature) {
        global_quadrature_points[pp] = entity.geometry().global(quadrature_point.position());
        ++pp;
      }

      // find source entities
      const auto source_entities = entity_search(global_quadrature_points);
      assert(source_entities.size() == global_quadrature_points.size());

      // loop over all quadrature points
      pp = 0;
      for (const auto& quadrature_point : quadrature) {
        const auto local_point = quadrature_point.position();
        const auto global_point = global_quadrature_points[pp];
        const auto quadrature_weight = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        const auto source_entity_ptr = source_entities[pp];
        const auto& source_entity = *source_entity_ptr;
        const auto local_source = source.localFunction(source_entity);
        const auto local_point_source = source_entity.geometry().local(global_point);
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        const auto local_source_value = local_source.evaluate(local_point_source);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (local_source_value * basis_values[ii]);
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix[ii][jj] += integration_element * quadrature_weight
                                    * (basis_values[ii] * basis_values[jj]);
          }
        }
        ++pp;
      } // loop over all quadrature points

      // compute local DoFs
      DynamicVector< R > local_DoFs(local_basis.size(), R(0));
      local_matrix.solve(local_DoFs, local_vector);

      // set local DoFs
      auto local_range = range.localFunction(entity);
      auto local_range_vector = local_range.vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs[ii]);

    } // walk the grid
  } // ... apply(...)
}; // class Generic


} // namespace ProlongationOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROLONGATIONS_HH
