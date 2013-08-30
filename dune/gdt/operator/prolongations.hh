// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PROLONGATIONS_HH
#define DUNE_GDT_OPERATOR_PROLONGATIONS_HH

#include <type_traits>
#include <vector>
#include <limits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>

#include <dune/stuff/grid/search.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>

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
  }
}; // class Generic


} // namespace ProlongationOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROLONGATIONS_HH
