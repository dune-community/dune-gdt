// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "../quadrature.hh"
#include "reconstructed_function.hh"
#include "slopelimiters.hh"

namespace Dune {
namespace GDT {

template <size_t dimDomain>
class StencilCoordinate
{
  StencilCoordinate(const FieldVector<int, dimDomain>& stencil_sizes)
    : stencil_sizes_(stencil_sizes)
    , offsets_(1)
  {
    for (size_t ii = 0; ii < dimDomain; ++ii)
      for (size_t jj = 0; jj < ii; ++jj)
        offsets_[ii] *= stencil_sizes_[jj];
  }

  FieldVector<int, dimDomain> get_zerobased_coords(const size_t vector_index) const
  {
    FieldVector<int, dimDomain> ret(vector_index);
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      for (size_t jj = 0; jj < ii; ++jj) {
        ret[ii] -= ret[jj] * offsets_[jj];
      }
      ret[ii] /= offsets_[ii];
      ret[ii] = ret[ii] % stencil_sizes_[ii];
    }
    return ret;
  }

  FieldVector<int, dimDomain> convert_to_relative_coords(FieldVector<int, dimDomain> zerobased_coords)
  {
    for (size_t ii = 0; ii < dimDomain; ++ii)
      zerobased_coords[ii] -= stencil_sizes_[ii] / 2;
    return zerobased_coords;
  }

  FieldVector<int, dimDomain> convert_to_zerobased_coords(FieldVector<int, dimDomain> relative_coords)
  {
    for (size_t ii = 0; ii < dimDomain; ++ii)
      relative_coords[ii] += stencil_sizes_[ii] / 2;
    return relative_coords;
  }

  size_t get_vector_index(const FieldVector<int, dimDomain>& relative_coords) const
  {
    size_t ret = 0;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret += (relative_coords[ii] + stencil_sizes_[ii] / 2) * offsets_[ii];
    return ret;
  }

  size_t calculate_num_values() const
  {
    size_t ret = 1;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret *= stencil_sizes_[ii];
    return ret;
  }

  size_t get_reordered_index(const size_t new_x_dir, const size_t vector_index)
  {
    auto coords = get_zerobased_coords(vector_index);
    auto new_stencil_sizes = stencil_sizes_;
    auto new_coords = coords;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      auto old_index = (new_x_dir + ii) % dimDomain;
      new_coords[ii] = coords[old_index];
      new_stencil_sizes[ii] = stencil_sizes_[old_index];
    }
    StencilCoordinate<dimDomain> new_stencil_coordinate(new_stencil_sizes);
    return new_stencil_coordinate.get_vector_index(new_coords);
  }

private:
  FieldVector<int, dimDomain> stencil_sizes_;
  FieldVector<int, dimDomain> offsets_;
}

template <class GridLayerType, class RangeType>
class Stencil
{
public:
  static const size_t dimDomain = GridLayerType::dimension;

  template <class EntityType, class BoundaryValuesFunctionType>
  Stencil(const std::vector<RangeType>& source_values,
          const FieldVector<int, dimDomain>& stencil_sizes,
          const GridLayerType& grid_layer,
          const BoundaryValuesFunctionType& boundary_values_function,
          const EntityType& entity)
    : source_values_(source_values)
    , stencil_sizes_(stencil_sizes)
    , grid_layer_(grid_layer)
    , boundary_values_()
    , stencil_values_(calculate_num_values(), nullptr)
    , valid_(get_stencil_values(entity, -1, FieldVector<int, dimDomain>(0), boundary_values_function))
  {
    static_assert(is_localizable_boundary_value<BoundaryValuesFunctionType>::value,
                  "BoundaryValuesFunctionType has to be derived from LocalizableBoundaryValueInterface");
#ifndef NDEBUG
    for (const auto& stencil_size : stencil_sizes_)
      if (stencil_size > 1)
        DUNE_THROW(NotImplemented, "Not implemented for stencil_sizes <= 1 in any direction!");
#endif
  }

  const RangeType& get_value(const FieldVector<int, dimDomain>& relative_coords) const
  {
    return *stencil_values_[get_vector_index[relative_coords]];
  }

  bool valid() const
  {
    return valid_;
  }

  size_t num_values() const
  {
    return stencil_values_.size();
  }

  const std::vector<const RangeType*>& values() const
  {
    return stencil_values_;
  }

private:
  template <class EntityType, class BoundaryValuesFunctionType>
  bool get_stencil_values(const EntityType& entity,
                          const int direction,
                          const FieldVector<int, dimDomain>& relative_coords,
                          const BoundaryValuesFunctionType& boundary_values_function)
  {
    bool ret = true;
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    stencil_values_[get_vector_index(relative_coords)] = source_values_[entity_index];
    std::vector<int> boundary_dirs;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity)) {
      const auto& intersection_index = intersection.indexInInside();
      if (!end_of_stencil(intersection_index, relative_coords)) {
        auto new_relative_coords = relative_coords;
        if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
          boundary_dirs.push_back(intersection_index);
          boundary_values_.push_back(boundary_values_function.local_function(entity)->evaluate(
              intersection, entity.geometry().local(intersection.geometry().center()), source_values_[entity_index]));
          while (!end_of_stencil(intersection_index, new_relative_coords)) {
            walk_in_direction(intersection_index, new_relative_coords);
            stencil_values_[get_vector_index(new_relative_coords)] = &boundary_values.back();
          }
        } else if (intersection.neighbor()
                   && direction_allowed(direction, intersection_index)) { // inner and periodic intersections
          const auto& outside = intersection.outside();
          walk_in_direction(intersection_index, new_relative_coords);
          ret = ret && get_stencil_values(outside, intersection_index, new_relative_coords, boundary_values_function);
        } else if (direction_allowed(direction, intersection_index) && !intersection.neighbor()
                   && !intersection.boundary()) { // processor boundary
          return false;
        }
      } // if (!end_of_stencil(...))
    } // intersections

    // TODO: improve multiple boundary handling, currently everything is filled with the boundary value in the first
    // direction
    assert(boundary_dirs.size() <= dimDomain);
    if (boundary_dirs.size() > 1) {
      walk_in_direction(boundary_dirs[0], relative_coords);
      const auto* boundary_value = stencil_values_[get_vector_index(relative_coords)];
      for (const RangeType* value : values)
        if (value == nullptr)
          value = boundary_value;
    }
    return ret;
  } // bool get_stencil_values(...)

  // walk in direction dir (increase or decrease relative_coords in that direction)
  static void walk_in_direction(const int dir, FieldVector<int, dimDomain>& relative_coords)
  {
    dir % 2 ? relative_coords[dir / 2]++ : relative_coords[dir / 2]--;
  }

  // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
  // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
  // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
  // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
  // without emitting new iterators).
  static constexpr bool direction_allowed(const int dir, const int new_dir)
  {
    return dir == -1 || new_dir == dir || new_dir / 2 > dir / 2;
  }

  static bool end_of_stencil(const int dir, const FieldVector<int, dimDomain>& relative_coords)
  {
    return (dir != -1 && std::abs(relative_coords[dir / 2]) >= stencil_sizes_[dir / 2] / 2);
  }

  const std::vector<RangeType>& source_values_;
  const FieldVector<int, dimDomain>& stencil_sizes_;
  const GridLayerType& grid_layer_;
  std::vector<RangeType> boundary_values_;
  std::vector<const RangeType*> stencil_values_;
  const bool valid_;
  FieldVector<int, dimDomain> offsets_;
}; // class Stencil


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
