// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PROJECTIONS_HH
#define DUNE_GDT_OPERATOR_PROJECTIONS_HH

#include <type_traits>
#include <vector>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>

namespace Dune {
namespace GDT {
namespace Operator {


template <class SourceImp, class RangeImp>
class DirichletProjection
{
  static_assert(std::is_base_of<Stuff::LocalizableFunction, SourceImp>::value,
                "SourceImp has to be derived from Stuff::LocalizableFunction");
  static_assert(std::is_same<RangeImp, DiscreteFunctionDefault<ContinuousLagrangeSpace::FemWrapper<
                                                                   typename RangeImp::SpaceType::GridPartType, 1,
                                                                   typename RangeImp::SpaceType::RangeFieldType, 1>,
                                                               typename RangeImp::VectorType>>::value,
                "RangeImp has wrong type!");

public:
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;

  typedef Stuff::GridboundaryInterface<typename RangeType::SpaceType::GridPartType::GridViewType> BoundaryInfoType;

  DirichletProjection(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  void apply(const SourceType& source, RangeType& range) const
  {
    // clear range
    range.vector()->backend() *= 0.0;
    // walk the grid
    const auto entity_it_end = range.space().gridPart().template end<0>();
    for (auto entity_it = range.space().gridPart().template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // only consider entities with boundary intersections
      if (entity.hasBoundaryIntersections()) {
        const auto source_local_function = source.localFunction(entity);
        auto range_local_function        = range.localFunction(entity);
        auto range_local_DoF_vector      = range_local_function.vector();
        const auto lagrangePointSet      = range.space().backend().lagrangePointSet(entity);
        // get the lagrange points' coordinates
        typedef typename RangeType::SpaceType::BackendType::LagrangePointSetType::CoordinateType
            LagrangePointCoordinateType;
        std::vector<LagrangePointCoordinateType> lagrangePoints(lagrangePointSet.nop(), LagrangePointCoordinateType(0));
        for (size_t ii = 0; ii < lagrangePointSet.nop(); ++ii)
          lagrangePoints[ii] = entity.geometry().global(lagrangePointSet.point(ii));
        // walk the intersections
        const auto intersection_it_end = range.space().gridPart().iend(entity);
        for (auto intersection_it = range.space().gridPart().ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          // only consider dirichlet boundary intersection
          if (boundary_info_.dirichlet(intersection)) {
            // loop over all lagrange points
            for (size_t ii = 0; ii < lagrangePointSet.nop(); ++ii) {
              // if dof lies on the boundary intersection
              if (Dune::Stuff::Grid::intersectionContains(intersection, lagrangePoints[ii])) {
                // set the corresponding target dof
                range_local_DoF_vector.set(ii, source_local_function.evaluate(lagrangePointSet.point(ii)));
              } // if dof lies on the boundary intersection
            } // loop over all lagrange points
          } // only consider dirichlet boundary intersection
        } // walk the intersections
      } // only consider entities with boundary intersection
    } // walk the grid
  } // ... apply(...) const

private:
  const BoundaryInfoType& boundary_info_;
}; // class DirichletProjection


} // namespace Operator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROJECTIONS_HH
