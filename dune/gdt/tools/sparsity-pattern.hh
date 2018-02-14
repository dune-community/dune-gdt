// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH
#define DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH

#include <dune/common/dynvector.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/pattern.hh>

#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 *  \brief Computes an element sparsity pattern, where the tes space determines the rows (outer) and the ansatz space
 *         determines the columns (inner).
 */
template <class TS, class GV, size_t Tr, size_t TrC, class AS, size_t Ar, size_t ArC>
XT::LA::SparsityPatternDefault
make_element_sparsity_pattern(const SpaceInterface<TS, GridView<GV>::dimension, Tr, TrC>& test_space,
                              const SpaceInterface<AS, GridView<GV>::dimension, Ar, ArC>& ansatz_space,
                              const GridView<GV>& grid_view) const
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    ansatz_space.mapper().global_indices(element, column_indices);
    for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
      for (size_t jj = 0; jj < ansatz_space.local_size(element); ++jj)
        pattern.insert(row_indices[ii], column_indices[jj]);
  }
  pattern.sort();
  return pattern;
} // ... make_element_sparsity_pattern(...)


template <class S, class GV, size_t r, size_t rC>
XT::LA::SparsityPatternDefault
make_element_sparsity_pattern(const SpaceInterface<S, GridView<GV>::dimension, r, rC>& space,
                              const GridView<GV>& grid_view) const
{
  return make_element_sparsity_pattern(space, space, grid_view);
}


template <class S, size_t d, size_t r, size_t rC>
XT::LA::SparsityPatternDefault make_element_sparsity_pattern(const SpaceInterface<S, d, r, rC>& space) const
{
  return make_element_sparsity_pattern(space, space, space.grid_view());
}


/**
 *  \brief Computes a coupling sparsity pattern, where the tes space determines the rows (outer) and the ansatz space
 *         determines the columns (inner).
 */
template <class TS, class GV, size_t Tr, size_t TrC, class AS, size_t Ar, size_t ArC>
XT::LA::SparsityPatternDefault
make_coupling_sparsity_pattern(const SpaceInterface<TS, GridView<GV>::dimension, Tr, TrC>& test_space,
                               const SpaceInterface<AS, GridView<GV>::dimension, Ar, ArC>& ansatz_space,
                               const GridView<GV>& grid_view) const
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    for (auto& intersection : intersections(grid_view, element)) {
      if (intersection.neighbor() && !intersection.boundary()) {
        const auto neighbour = intersection.outside();
        ansatz_space.mapper().global_indices(neighbour, column_indices);
        for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
          for (size_t jj = 0; jj < ansatz_space.mapper().local_size(neighbor); ++jj)
            pattern.insert(row_indices[ii], column_indices[jj]);
      }
    }
  }
  pattern.sort();
  return pattern;
} // ... make_coupling_sparsity_pattern(...)


template <class S, class GV, size_t r, size_t rC>
XT::LA::SparsityPatternDefault
make_coupling_sparsity_pattern(const SpaceInterface<S, GridView<GV>::dimension, r, rC>& space,
                               const GridView<GV>& grid_view) const
{
  return make_coupling_sparsity_pattern(space, space, grid_view);
}


template <class S, size_t d, size_t r, size_t rC>
XT::LA::SparsityPatternDefault make_coupling_sparsity_pattern(const SpaceInterface<S, d, r, rC>& space) const
{
  return make_coupling_sparsity_pattern(space, space, space.grid_view());
}


/**
 *  \brief Computes an element and coupling sparsity pattern, where the tes space determines the rows (outer) and the
 *         ansatz space determines the columns (inner).
 */
template <class TS, class GV, size_t Tr, size_t TrC, class AS, size_t Ar, size_t ArC>
XT::LA::SparsityPatternDefault
make_element_and_coupling_sparsity_pattern(const SpaceInterface<TS, GridView<GV>::dimension, Tr, TrC>& test_space,
                                           const SpaceInterface<AS, GridView<GV>::dimension, Ar, ArC>& ansatz_space,
                                           const GridView<GV>& grid_view) const
{
  XT::LA::SparsityPatternDefault pattern(test_space.mapper().size());
  DynamicVector<size_t> row_indices(test_space.mapper().max_local_size(), 0);
  DynamicVector<size_t> column_indices(ansatz_space.mapper().max_local_size(), 0);
  for (auto&& element : elements(grid_view)) {
    test_space.mapper().global_indices(element, row_indices);
    ansatz_space.mapper().global_indices(element, column_indices);
    for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
      for (size_t jj = 0; jj < ansatz_space.local_size(element); ++jj)
        pattern.insert(row_indices[ii], column_indices[jj]);
    for (auto& intersection : intersections(grid_view, element)) {
      if (intersection.neighbor() && !intersection.boundary()) {
        const auto neighbour = intersection.outside();
        ansatz_space.mapper().global_indices(neighbour, column_indices);
        for (size_t ii = 0; ii < test_space.mapper().local_size(element); ++ii)
          for (size_t jj = 0; jj < ansatz_space.mapper().local_size(neighbor); ++jj)
            pattern.insert(row_indices[ii], column_indices[jj]);
      }
    }
  }
  pattern.sort();
  return pattern;
} // ... make_element_and_coupling_sparsity_pattern(...)


template <class S, class GV, size_t r, size_t rC>
XT::LA::SparsityPatternDefault
make_element_and_coupling_sparsity_pattern(const SpaceInterface<S, GridView<GV>::dimension, r, rC>& space,
                                           const GridView<GV>& grid_view) const
{
  return make_element_and_coupling_sparsity_pattern(space, space, grid_view);
}


template <class S, size_t d, size_t r, size_t rC>
XT::LA::SparsityPatternDefault
make_element_and_coupling_sparsity_pattern(const SpaceInterface<S, d, r, rC>& space) const
{
  return make_element_and_coupling_sparsity_pattern(space, space, space.grid_view());
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_SPARSITY_PATTERN_HH
