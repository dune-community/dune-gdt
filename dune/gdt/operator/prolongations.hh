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

namespace Dune {
namespace GDT {
namespace ProlongationOperator {


template< class GridPartType >
class L2
{
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  L2(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {}

  template< class SourceSpaceType, class VS, class RangeSpaceType, class VR >
  void apply(const ConstDiscreteFunction< SourceSpaceType, VS >& source,
             DiscreteFunction< RangeSpaceType, VR >& range) const
  {
    // check
    static_assert(std::is_same< typename SourceSpaceType::RangeFieldType,
                                typename RangeSpaceType::RangeFieldType >::value, "Types do not match!");
    static_assert(SourceSpaceType::dimRange == RangeSpaceType::dimRange, "Dimensions do not match!");
    static_assert(SourceSpaceType::dimRangeCols == RangeSpaceType::dimRangeCols, "Dimensions do not match!");
    static_assert(SourceSpaceType::dimRangeCols == 1, "Not implemented yet!");
    typedef typename RangeSpaceType::BaseFunctionSetType::DomainType DomainType;
    typedef typename RangeSpaceType::BaseFunctionSetType::RangeType RangeType;
    typedef typename RangeSpaceType::BaseFunctionSetType::RangeFieldType RangeFieldType;
    // clear
    range.vector().backend() *= RangeFieldType(0);
    // create search in the source grid part
    typedef typename SourceSpaceType::GridPartType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch< SourceGridViewType > EntitySearch;
    EntitySearch entity_search(source.space().gridPart()->gridView());
    // walk the grid
    RangeType source_value(0);
    std::vector< RangeType > basis_values(range.space().mapper().maxNumDofs());
    const auto entity_it_end = grid_part_.template end< 0 >();
    for (auto entity_it = grid_part_.template begin< 0 >();
         entity_it != entity_it_end;
         ++entity_it) {
      // prepare
      const auto& entity = *entity_it;
      const auto local_basis = range.space().baseFunctionSet(entity);
      auto local_range = range.local_discrete_function(entity);
      DynamicMatrix< RangeFieldType > local_matrix(local_basis.size(), local_basis.size(), RangeFieldType(0));
      DynamicVector< RangeFieldType > local_vector(local_basis.size(), RangeFieldType(0));
      // create quadrature
      const size_t quadrature_order = local_range.order();
      const auto& quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                   2*quadrature_order + 1);
      // get global quadrature points
      std::vector< DomainType > quadrature_points;
      for (const auto& quadrature_point : quadrature)
        quadrature_points.emplace_back(entity.geometry().global(quadrature_point.position()));
      // get source entities
      const auto source_entities = entity_search(quadrature_points);
      assert(source_entities.size() == quadrature_points.size());
      // loop over all quadrature points
      size_t pp = 0;
      for (const auto& quadrature_point : quadrature) {
        const auto local_point = quadrature_point.position();
        const auto quadrature_weight = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        // evaluate source
        const auto source_entity_ptr = source_entities[pp];
        const auto& source_entity = *source_entity_ptr;
        const auto local_source = source.local_function(source_entity);
        local_source->evaluate(source_entity.geometry().local(entity.geometry().global(local_point)), source_value);
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (source_value * basis_values[ii]);
          auto& local_matrix_row = local_matrix[ii];
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix_row[jj] += integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]);
          }
        }
        ++pp;
      } // loop over all quadrature points
      // compute local DoFs
      DynamicVector< RangeFieldType > local_DoFs(local_basis.size(), 0);
      local_matrix.solve(local_DoFs, local_vector);
      // set local DoFs
      auto local_range_vector = local_range.vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs[ii]);
    } // walk the grid
  } // ... apply(...)

private:
  const GridPartType& grid_part_;
}; // class L2


template< class GridPartType >
class Generic
{
public:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  Generic(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {}

  template< class SGP, int sp, class R, int r, class SV, class RGP, int rp, class RV >
  void apply(const ConstDiscreteFunction< ContinuousLagrangeSpace::FemWrapper< SGP, sp, R, r >, SV >& source,
             DiscreteFunction< ContinuousLagrangeSpace::FemWrapper< RGP, rp, R, r >, RV >& range) const
  {
    // create search in the source grid part
    typedef ConstDiscreteFunction< ContinuousLagrangeSpace::FemWrapper< SGP, sp, R, r >, SV > SourceType;
    typedef DiscreteFunction< ContinuousLagrangeSpace::FemWrapper< RGP, rp, R, r >, RV >      TargetType;
    typedef typename SourceType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange = SourceType::dimRangeRows;

    typedef typename SourceType::SpaceType::GridPartType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch< SourceGridViewType > EntitySearch;
    EntitySearch entity_search(source.space().gridPart()->gridView());

    // set all dofs to infinity
    const auto infinity = std::numeric_limits< RangeFieldType >::infinity();
    for (size_t ii = 0; ii < range.vector().size(); ++ii)
      range.vector().set(ii, infinity);

    // walk the grid
    const auto entity_it_end = grid_part_.template end< 0 >();
    for (auto entity_it = grid_part_.template begin< 0 >();
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

      auto local_range = range.local_discrete_function(entity);
      auto local_range_vector = local_range.vector();

      size_t kk = 0;
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        if (std::isinf(local_range_vector.get(kk))) {
          const auto& global_point = lagrange_points[ii];
          // evaluate source function
          const auto source_entity_ptr = source_entities[ii];
          const auto& source_entity = *source_entity_ptr;
          const auto local_source_point = source_entity.geometry().local(global_point);
          const auto local_source = source.local_function(source_entity);
          const auto source_value = local_source->evaluate(local_source_point);
          for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
            local_range_vector.set(kk, source_value[jj]);
        }
        else
          kk += dimRange;
      }
    } // walk the grid
  } // ... apply(...)

private:
  const GridPartType& grid_part_;
}; // class Generic


} // namespace ProlongationOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROLONGATIONS_HH
