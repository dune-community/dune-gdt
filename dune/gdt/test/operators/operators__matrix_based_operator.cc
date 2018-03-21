// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

using namespace Dune;
using namespace Dune::GDT;


GTEST_TEST(matrix_operator, make_matrix_operator)
{
  using G = YASP_1D_EQUIDISTANT_OFFSET;
  auto grid = XT::Grid::make_cube_grid<G>(0, 1, 2);
  auto grid_view = grid.leaf_view();
  using GV = decltype(grid_view);
  const ContinuousLagrangeSpace<GV, 1> space(grid_view);

  using MatrixType = XT::LA::CommonDenseMatrix<double>;
  auto matrix = XT::LA::eye_matrix<MatrixType>(space.mapper().size());

  // given a matrix
  make_matrix_operator(grid_view, space, space, matrix);
  make_matrix_operator(grid_view, space, matrix);
  make_matrix_operator(space, matrix);
  // given a pattern, creates a matrix
  XT::LA::SparsityPatternDefault pattern(3);
  pattern.insert(0, 0);
  pattern.insert(1, 1);
  pattern.insert(2, 2);
  pattern.sort();
  make_matrix_operator<MatrixType>(grid_view, space, space, pattern);
  make_matrix_operator<MatrixType>(grid_view, space, pattern);
  make_matrix_operator<MatrixType>(space, pattern);
  // given a stencil, creates a pattern, creates a matrix
  make_matrix_operator<MatrixType>(grid_view, space, space, Stencil::intersection);
  make_matrix_operator<MatrixType>(grid_view, space, Stencil::intersection);
  make_matrix_operator<MatrixType>(space, Stencil::intersection);
  // given nothing, uses largest possible stencil, creates pattern, creates matrix
  make_matrix_operator<MatrixType>(grid_view, space, space);
  make_matrix_operator<MatrixType>(grid_view, space);
  make_matrix_operator<MatrixType>(space);
}
