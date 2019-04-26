// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH
#define DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH

#include <dune/common/dynvector.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/pattern.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


/**
 *  \brief Computes an element sparsity pattern, where the test space determines the rows (outer) and the ansatz space
 *         determines the columns (inner).
 */
template <class TGV, size_t t_r, size_t t_rC, class TR, class AGV, size_t a_r, size_t a_rC, class AR, class GV>
XT::LA::SparsityPatternDefault make_element_sparsity_pattern(const SpaceInterface<TGV, t_r, t_rC, TR>& test_space,
                                                             const SpaceInterface<AGV, a_r, a_rC, AR>& ansatz_space,
                                                             const GV& grid_view)
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    ansatz_space.mapper().global_indices(element, column_indices);
    for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
      for (size_t jj = 0; jj < ansatz_space.mapper().local_size(element); ++jj)
        pattern.insert(row_indices[ii], column_indices[jj]);
  }
  pattern.sort();
  return pattern;
} // ... make_element_sparsity_pattern(...)


template <class SGV, size_t r, size_t rC, class R, class GV>
XT::LA::SparsityPatternDefault make_element_sparsity_pattern(const SpaceInterface<SGV, r, rC, R>& space,
                                                             const GV& grid_view)
{
  return make_element_sparsity_pattern(space, space, grid_view);
}


template <class GV, size_t r, size_t rC, class R>
XT::LA::SparsityPatternDefault make_element_sparsity_pattern(const SpaceInterface<GV, r, rC, R>& space)
{
  return make_element_sparsity_pattern(space, space, space.grid_view());
}


/**
 *  \brief Computes a coupling sparsity pattern, where the test space determines the rows (outer) and the ansatz space
 *         determines the columns (inner).
 */
template <class TGV, size_t t_r, size_t t_rC, class TR, class AGV, size_t a_r, size_t a_rC, class AR, class GV>
XT::LA::SparsityPatternDefault
make_intersection_sparsity_pattern(const SpaceInterface<TGV, t_r, t_rC, TR>& test_space,
                                   const SpaceInterface<AGV, a_r, a_rC, AR>& ansatz_space,
                                   const GV& grid_view)
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    for (auto& intersection : intersections(grid_view, element)) {
      if (intersection.neighbor()) {
        const auto neighbour = intersection.outside();
        ansatz_space.mapper().global_indices(neighbour, column_indices);
        for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
          for (size_t jj = 0; jj < ansatz_space.mapper().local_size(neighbour); ++jj)
            pattern.insert(row_indices[ii], column_indices[jj]);
      }
    }
  }
  pattern.sort();
  return pattern;
} // ... make_intersection_sparsity_pattern(...)


template <class SGV, size_t r, size_t rC, class R, class GV>
XT::LA::SparsityPatternDefault make_intersection_sparsity_pattern(const SpaceInterface<SGV, r, rC, R>& space,
                                                                  const GV& grid_view)
{
  return make_intersection_sparsity_pattern(space, space, grid_view);
}


template <class GV, size_t r, size_t rC, class R>
XT::LA::SparsityPatternDefault make_intersection_sparsity_pattern(const SpaceInterface<GV, r, rC, R>& space)
{
  return make_intersection_sparsity_pattern(space, space, space.grid_view());
}


/**
 *  \brief Computes an element and coupling sparsity pattern, where the test space determines the rows (outer) and the
 *         ansatz space determines the columns (inner).
 */
template <class TGV, size_t t_r, size_t t_rC, class TR, class AGV, size_t a_r, size_t a_rC, class AR, class GV>
XT::LA::SparsityPatternDefault
make_element_and_intersection_sparsity_pattern(const SpaceInterface<TGV, t_r, t_rC, TR>& test_space,
                                               const SpaceInterface<AGV, a_r, a_rC, AR>& ansatz_space,
                                               const GV& grid_view)
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    ansatz_space.mapper().global_indices(element, column_indices);
    for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
      for (size_t jj = 0; jj < ansatz_space.mapper().local_size(element); ++jj)
        pattern.insert(row_indices[ii], column_indices[jj]);
    for (auto& intersection : intersections(grid_view, element)) {
      if (intersection.neighbor()) {
        const auto neighbour = intersection.outside();
        ansatz_space.mapper().global_indices(neighbour, column_indices);
        for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
          for (size_t jj = 0; jj < ansatz_space.mapper().local_size(neighbour); ++jj)
            pattern.insert(row_indices[ii], column_indices[jj]);
      }
    }
  }
  pattern.sort();
  return pattern;
} // ... make_element_and_intersection_sparsity_pattern(...)


template <class SGV, size_t r, size_t rC, class R, class GV>
XT::LA::SparsityPatternDefault
make_element_and_intersection_sparsity_pattern(const SpaceInterface<SGV, r, rC, R>& space, const GV& grid_view)
{
  return make_element_and_intersection_sparsity_pattern(space, space, grid_view);
}


template <class GV, size_t r, size_t rC, class R>
XT::LA::SparsityPatternDefault make_element_and_intersection_sparsity_pattern(const SpaceInterface<GV, r, rC, R>& space)
{
  return make_element_and_intersection_sparsity_pattern(space, space, space.grid_view());
}


template <class TGV, size_t t_r, size_t t_rC, class TR, class AGV, size_t a_r, size_t a_rC, class AR, class GV>
XT::LA::SparsityPatternDefault make_sparsity_pattern(const SpaceInterface<TGV, t_r, t_rC, TR>& test_space,
                                                     const SpaceInterface<AGV, a_r, a_rC, AR>& ansatz_space,
                                                     const GV& grid_view,
                                                     const Stencil stencil)
{
  if (stencil == Stencil::element)
    return make_element_sparsity_pattern(test_space, ansatz_space, grid_view);
  else if (stencil == Stencil::intersection)
    return make_intersection_sparsity_pattern(test_space, ansatz_space, grid_view);
  else if (stencil == Stencil::element_and_intersection)
    return make_element_and_intersection_sparsity_pattern(test_space, ansatz_space, grid_view);
  else
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "Unknown Stencil encountered (see below), add an appropriate method!"
                   << "\n   stencil = " << stencil);
} // ... make_sparsity_pattern(...)


template <class SGV, size_t r, size_t rC, class R, class GV>
XT::LA::SparsityPatternDefault
make_sparsity_pattern(const SpaceInterface<SGV, r, rC, R>& space, const GV& grid_view, const Stencil stencil)
{
  return make_sparsity_pattern(space, space, grid_view, stencil);
}


template <class GV, size_t r, size_t rC, class R>
XT::LA::SparsityPatternDefault make_sparsity_pattern(const SpaceInterface<GV, r, rC, R>& space, const Stencil stencil)
{
  return make_sparsity_pattern(space, space.grid_view(), stencil);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH
