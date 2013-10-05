// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PROJECTIONS_HH
#define DUNE_GDT_OPERATOR_PROJECTIONS_HH

#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>

namespace Dune {
namespace GDT {
namespace ProjectionOperator {


template <class GridPartType>
class L2
{
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  L2(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class R, int r, class SpaceType, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& source,
             DiscreteFunction<SpaceType, V>& range) const
  {
    // check
    static_assert(std::is_same<R, typename SpaceType::RangeFieldType>::value, "Types do not match!");
    static_assert(SpaceType::dimRange == r, "Dimensions do not match!");
    static_assert(SpaceType::dimRangeCols == 1, "Not implemented yet!");
    typedef typename SpaceType::BaseFunctionSetType::RangeType RangeType;
    // clear
    range.vector().backend() *= R(0);
    // walk the grid
    RangeType source_value(0);
    std::vector<RangeType> basis_values(range.space().mapper().maxNumDofs());
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      // prepare
      const auto& entity      = *entity_it;
      const auto local_basis  = range.space().baseFunctionSet(entity);
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      DynamicMatrix<R> local_matrix(local_basis.size(), local_basis.size(), R(0));
      DynamicVector<R> local_vector(local_basis.size(), R(0));
      // create quadrature
      const size_t quadrature_order = std::max(local_source->order(), local_range.order());
      const auto& quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * quadrature_order + 1);
      // loop over all quadrature points
      for (const auto& quadrature_point : quadrature) {
        const auto local_point         = quadrature_point.position();
        const auto quadrature_weight   = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        local_source->evaluate(local_point, source_value);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (source_value * basis_values[ii]);
          auto& local_matrix_row = local_matrix[ii];
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix_row[jj] += integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]);
          }
        }
      } // loop over all quadrature points
      // compute local DoFs
      DynamicVector<R> local_DoFs(local_basis.size(), 0);
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


template <class GridPartImp>
class Dirichlet
{
public:
  typedef GridPartImp GridPartType;
  typedef Stuff::GridboundaryInterface<typename GridPartType::IntersectionType> BoundaryInfoType;

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

public:
  Dirichlet(const GridPartType& grid_part, const BoundaryInfoType& boundary_info)
    : grid_part_(grid_part)
    , boundary_info_(boundary_info)
  {
  }

  template <class R, class GP, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, p, R, 1, 1>, V>& range) const
  {
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, p, R, 1, 1>, V> RangeType;
    // clear range
    range.vector().backend() *= 0.0;
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // only consider entities with boundary intersections
      if (entity.hasBoundaryIntersections()) {
        const auto source_local_function = source.local_function(entity);
        auto range_local_function        = range.local_discrete_function(entity);
        auto range_local_DoF_vector      = range_local_function.vector();
        const auto lagrangePointSet      = range.space().backend().lagrangePointSet(entity);
        // get the lagrange points' coordinates
        typedef typename RangeType::SpaceType::BackendType::LagrangePointSetType::CoordinateType
            LagrangePointCoordinateType;
        std::vector<LagrangePointCoordinateType> lagrangePoints(lagrangePointSet.nop(), LagrangePointCoordinateType(0));
        for (size_t ii = 0; ii < lagrangePointSet.nop(); ++ii)
          lagrangePoints[ii] = entity.geometry().global(lagrangePointSet.point(ii));
        // walk the intersections
        const auto intersection_it_end = grid_part_.iend(entity);
        for (auto intersection_it = grid_part_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          // only consider dirichlet boundary intersection
          if (boundary_info_.dirichlet(intersection)) {
            // loop over all lagrange points
            for (size_t ii = 0; ii < lagrangePointSet.nop(); ++ii) {
              // if dof lies on the boundary intersection
              if (Dune::Stuff::Grid::intersectionContains(intersection, lagrangePoints[ii])) {
                // set the corresponding target dof
                range_local_DoF_vector.set(ii, source_local_function->evaluate(lagrangePointSet.point(ii)));
              } // if dof lies on the boundary intersection
            } // loop over all lagrange points
          } // only consider dirichlet boundary intersection
        } // walk the intersections
      } // only consider entities with boundary intersection
    } // walk the grid
  } // ... apply(...) const

  template <class R, class GP, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, R, 1, 1>, V>& range) const
  {
    typedef ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, R, 1, 1> RangeSpaceType;
    // checks
    static_assert(dimDomain == 2, "This does not work for other dimensions!");
    static_assert(p == 1, "Not tested for higher polynomial orders!");
    FieldVector<R, 1> tmp_source_value(R(0));
    // clear range
    range.vector().backend() *= 0.0;
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // only work if we are at the boundary
      if (entity.hasBoundaryIntersections()) {
        // get the local quantities
        const auto local_source      = source.local_function(entity);
        auto local_range             = range.local_discrete_function(entity);
        auto& local_range_DoF_vector = local_range.vector();
        // get local finite elements
        // * we are a CG space
        const auto cg_finite_element      = range.space().backend().finiteElement(entity);
        const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
        // * but we also need a local DG finite element
        typedef Dune::DGLocalFiniteElement<typename RangeSpaceType::Traits::FiniteElementType> DgFiniteElementType;
        const DgFiniteElementType dg_finite_element(entity.type(), p);
        const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
        assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");

        // first we loop over all vertices of the entity
        std::vector<DomainType> global_vertices(entity.template count<dimDomain>(), DomainType(0));
        std::vector<size_t> local_DoF_ids(global_vertices.size(), 0);
        for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
          // get the vertex
          const auto vertexPtr             = entity.template subEntity<dimDomain>(local_vertex_id);
          const auto& vertex               = *vertexPtr;
          global_vertices[local_vertex_id] = vertex.geometry().center();
          // find the global DoF id to this vertex, therefore
          // loop over all local DoFs
          for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
            const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
            if (entity_cg_local_key.subEntity() == local_vertex_id) {
              // get the local DoF to this vertex
              const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
              assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
              local_DoF_ids[local_vertex_id] = entity_dg_local_key.index();
              // there must be one and only one for a polorder 1 lagrange basis
              break;
            }
          } // loop over all local DoFs
        } // loop over all vertices of the entity

        // then we walk the intersections
        const auto intersection_it_end = grid_part_.iend(entity);
        for (auto intersection_it = grid_part_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          if (boundary_info_.dirichlet(intersection)) {
            const auto& intersection_geometry = intersection.geometry();
            // and walk its corners (i.e. the vertices in 2d)
            for (size_t local_intersection_corner_id = 0;
                 int(local_intersection_corner_id) < intersection_geometry.corners();
                 ++local_intersection_corner_id) {
              const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
              // to check which vertex this corner is
              // loop over all vertices of the entity again
              for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
                // and check for equality
                if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id])) {
                  // this vertex is on the dirichlet boundary
                  // * so we evaluate the source
                  local_source->evaluate(global_intersection_corner, tmp_source_value);
                  // * and set the corresponding local DoF
                  local_range_DoF_vector.set(local_DoF_ids[local_vertex_id], tmp_source_value[0]);
                }
              } // loop over all vertices of the entity
            } // walk its corners
          }
        } // walk the intersections
      }
    } // walk the grid
  } // ... apply(...) const

private:
  const GridPartType& grid_part_;
  const BoundaryInfoType& boundary_info_;
}; // class Dirichlet


} // namespace ProjectionOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROJECTIONS_HH
