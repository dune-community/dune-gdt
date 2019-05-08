// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_HH
#define DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_HH

#include <vector>
#include <set>
#include <limits>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/boundaryinfo/types.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/dg/default.hh>
#include <dune/gdt/playground/spaces/block.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class OswaldInterpolationOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class OswaldInterpolationOperatorTraits
{
public:
  typedef OswaldInterpolationOperator<GridLayerImp, FieldImp> derived_type;
  typedef GridLayerImp GridLayerType;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class FieldImp>
class OswaldInterpolationOperator
  : public OperatorInterface<internal::OswaldInterpolationOperatorTraits<GridLayerImp, FieldImp>>
{
public:
  typedef internal::OswaldInterpolationOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename Traits::FieldType FieldType;
  static const size_t dimDomain = GridLayerType::dimension;

private:
  typedef XT::Grid::extract_entity_t<GridLayerType> E;
  typedef typename GridLayerType::ctype D;
  static const constexpr size_t d = dimDomain;

public:
  OswaldInterpolationOperator(
      const GridLayerType& grd_layr,
      const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info)
    : grid_layer_(grd_layr)
    , boundary_info_(boundary_info)
  {}

  template <class GL, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, FieldType, 1>& source,
             DiscreteFunction<DiscontinuousLagrangeSpace<GL, 1, FieldType>, V>& range) const
  {
    apply_p1_dg(source, range);
  }

  template <class GL, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, FieldType, 1>& source,
             DiscreteFunction<BlockSpace<DiscontinuousLagrangeSpace<GL, 1, FieldType>>, V>& range) const
  {
    apply_p1_dg(source, range);
  }

private:
  template <class SourceType, class RangeType>
  void apply_p1_dg(const SourceType& source, RangeType& range) const
  {
    range.vector() *= 0.;
    // data structures we need
    // * a map from a global vertex index to global DoF indices
    //   given a vertex, one obtains a set of all global DoF ids, which are associated with this vertex
    std::map<size_t, std::set<size_t>> global_vertex_id_to_global_DoF_id_map;
    // * a map from a global DoF index to the global index of its associated vertex
    std::vector<size_t> global_DoF_id_to_global_vertex_id_map(range.space().mapper().size());
    // * a set to hold the global id of all boundary vertices
    std::set<size_t> boundary_vertices;

    // walk the grid to create the maps explained above and to find the boundary vertices
    for (auto&& entity : elements(grid_layer_)) {
      const size_t num_vertices = entity.subEntities(dimDomain);
      assert(entity.dimension == dimDomain);
      const auto basis = range.space().base_function_set(entity);
      if (basis.size() != num_vertices)
        DUNE_THROW(XT::Common::Exceptions::internal_error, "basis.size() = " << basis.size());

      // loop over all vertices of the entitity, to find their associated global DoF indices
      for (auto local_vertex_id : XT::Common::value_range(num_vertices)) {
        const auto vertex = entity.template subEntity<dimDomain>(static_cast<int>(local_vertex_id));
        const auto global_vertex_id = grid_layer_.indexSet().index(vertex);
        const auto vertex_center = vertex.geometry().center();
        // find the local basis function which corresponds to this vertex
        const auto basis_values = basis.evaluate(entity.geometry().local(vertex_center));
        if (basis_values.size() != num_vertices)
          DUNE_THROW(XT::Common::Exceptions::internal_error, "basis_values.size() = " << basis_values.size());
        size_t ones = 0;
        size_t zeros = 0;
        size_t failures = 0;
        size_t local_DoF_index = 0;
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          if (std::abs(basis_values[ii][0] - 1.0) < 1e-14) {
            local_DoF_index = ii;
            ++ones;
          } else if (std::abs(basis_values[ii][0] - 0.0) < 1e-14)
            ++zeros;
          else
            ++failures;
        }
        if (ones != 1 || zeros != (basis.size() - 1) || failures > 0) {
          std::stringstream ss;
          ss << "ones = " << ones << ", zeros = " << zeros << ", failures = " << failures
             << ", num_vertices = " << num_vertices << ", entity " << grid_layer_.indexSet().index(entity)
             << ", vertex " << local_vertex_id << ": [ " << vertex_center << "], ";
          XT::Common::print(basis_values, "basis_values", ss);
          DUNE_THROW(XT::Common::Exceptions::internal_error, ss.str());
        }
        // now we know that the local DoF index of this vertex is ii
        const size_t global_DoF_index = range.space().mapper().mapToGlobal(entity, local_DoF_index);
        global_DoF_id_to_global_vertex_id_map[global_DoF_index] = global_vertex_id;
        global_vertex_id_to_global_DoF_id_map[global_vertex_id].insert(global_DoF_index);
      } // loop over all vertices

      static const constexpr auto dirichlet = XT::Grid::DirichletBoundary();
      // in order to determine the boundary vertices, we need to
      // loop over all intersections
      const auto intersectionEndIt = grid_layer_.iend(entity);
      for (auto intersectionIt = grid_layer_.ibegin(entity); intersectionIt != intersectionEndIt; ++intersectionIt) {
        const auto& intersection = *intersectionIt;
        if (boundary_info_.type(intersection) == dirichlet) {
          const auto& intersection_geometry = intersection.geometry();
          for (auto local_intersection_corner_id : XT::Common::value_range(intersection_geometry.corners())) {
            const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
            // now, we need to find the entity's vertex this intersection's corner point equals to, so we
            // loop over all vertices of the entity
            for (auto local_vertex_id : XT::Common::value_range(num_vertices)) {
              const auto vertex = entity.template subEntity<dimDomain>(static_cast<int>(local_vertex_id));
              const auto global_vertex_id = grid_layer_.indexSet().index(vertex);
              const auto vertex_center = vertex.geometry().center();
              if (XT::Common::FloatCmp::eq(global_intersection_corner, vertex_center))
                boundary_vertices.insert(global_vertex_id);
            } // loop over all vertices of the entity
          } // loop over all intersection corners
        } // if (intersection.boundary() && !intersection.neighbor())
      } // loop over all intersections
    } // walk the grid for the first time

    // walk the grid for the second time
    for (auto&& entity : elements(grid_layer_)) {
      const auto num_vertices = entity.subEntities(dimDomain);
      // get the local functions
      const auto local_source = source.local_function(entity);

      // * loop over all local DoFs
      for (auto local_vertex_id : XT::Common::value_range(num_vertices)) {
        const size_t global_DoF_index =
            range.space().mapper().mapToGlobal(entity, XT::Common::numeric_cast<size_t>(local_vertex_id));
        const size_t global_vertex_id = global_DoF_id_to_global_vertex_id_map[global_DoF_index];
        // if we are on the domain boundary
        if (boundary_vertices.count(global_vertex_id)) {
          // set the dof to zero (we have dirichlet zero)
          range.vector().set_entry(global_DoF_index, FieldType(0));
        } else {
          // do the oswald projection
          const size_t num_DoFS_per_vertex = global_vertex_id_to_global_DoF_id_map[global_vertex_id].size();
          // * get the source value
          const auto vertex = entity.template subEntity<dimDomain>(local_vertex_id).geometry().center();
          const FieldType source_value = local_source->evaluate(entity.geometry().local(vertex));
          // * and add it to all target DoFs
          for (size_t target_global_DoF_id : global_vertex_id_to_global_DoF_id_map[global_vertex_id])
            range.vector().add_to_entry(target_global_DoF_id, source_value / num_DoFS_per_vertex);
        } // if (boundary_vertices.find(global_vertex_id))
      } // loop over all local DoFs
    } // walk the grid for the second time
  } // ... apply_p1_dg(...)

  const GridLayerType& grid_layer_;
  const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info_;
}; // class OswaldInterpolationOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_HH
