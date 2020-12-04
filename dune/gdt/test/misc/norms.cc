// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#  define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
#  define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/print.hh>
#include <dune/gdt/norms.hh>

using namespace Dune;
using namespace Dune::GDT;

using G = YASP_2D_EQUIDISTANT_OFFSET;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using F = double;
using M = XT::LA::IstlRowMajorSparseMatrix<F>;
using V = XT::LA::IstlDenseVector<F>;

GTEST_TEST(norms, l2_norm)
{
  auto grid = XT::Grid::make_cube_grid<G>();
  auto grid_view = grid.leaf_view();

  XT::Functions::GridFunction<E> f(1);
  XT::Functions::GridFunction<E, 3> g({1, 1, 1});
  XT::Functions::GridFunction<E, 3, 3> h(1);

  EXPECT_DOUBLE_EQ(1, l2_norm(grid_view, 1.));
  EXPECT_DOUBLE_EQ(std::sqrt(2), l2_norm(grid_view, 1., 2.));
  EXPECT_DOUBLE_EQ(1, l2_norm(grid_view, 1., f));
  EXPECT_DOUBLE_EQ(1, l2_norm(grid_view, f));
  EXPECT_DOUBLE_EQ(std::sqrt(2), l2_norm(grid_view, f, 2.));
  EXPECT_DOUBLE_EQ(1, l2_norm(grid_view, f, f));
  EXPECT_DOUBLE_EQ(std::sqrt(3), l2_norm(grid_view, g));
  EXPECT_DOUBLE_EQ(std::sqrt(3), l2_norm(grid_view, g, {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}));
  EXPECT_DOUBLE_EQ(std::sqrt(3), l2_norm(grid_view, g, h));
}


GTEST_TEST(norms, h1_semi_norm)
{
  auto grid = XT::Grid::make_cube_grid<G>();
  auto grid_view = grid.leaf_view();

  XT::Functions::GridFunction<E> f(1);
  XT::Functions::GridFunction<E, 3> g({1, 1, 1});
  XT::Functions::GridFunction<E, G::dimension, G::dimension> hd(1);

  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, 1.));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, 1., 2.));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, 1., f));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, f));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, f, 2.));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, f, f));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, g));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, g, {{{1, 0}, {0, 1}}}));
  EXPECT_DOUBLE_EQ(0, h1_semi_norm(grid_view, g, hd));
}


GTEST_TEST(norms, h1_norm)
{
  auto grid = XT::Grid::make_cube_grid<G>();
  auto grid_view = grid.leaf_view();

  XT::Functions::GridFunction<E> f(1);
  XT::Functions::GridFunction<E, 3> g({1, 1, 1});
  XT::Functions::GridFunction<E, 3, 3> h(1);
  XT::Functions::GridFunction<E, G::dimension, G::dimension> hd(1);

  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, 1.));
  EXPECT_DOUBLE_EQ(std::sqrt(2), h1_norm(grid_view, 1., 2.));
  EXPECT_DOUBLE_EQ(std::sqrt(2), h1_norm(grid_view, 1., 2., 2.));
  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, 1., f));
  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, 1., f, f));
  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, f));
  EXPECT_DOUBLE_EQ(std::sqrt(2), h1_norm(grid_view, f, 2.));
  EXPECT_DOUBLE_EQ(std::sqrt(2), h1_norm(grid_view, f, 2., 2.));
  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, f, f));
  EXPECT_DOUBLE_EQ(1, h1_norm(grid_view, f, f, f));
  EXPECT_DOUBLE_EQ(std::sqrt(3), h1_norm(grid_view, g));
  EXPECT_DOUBLE_EQ(std::sqrt(3), h1_norm(grid_view, g, 1.));
  EXPECT_DOUBLE_EQ(std::sqrt(3), h1_norm(grid_view, g, 1., 1.));
  EXPECT_DOUBLE_EQ(std::sqrt(3), h1_norm(grid_view, g, {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}, {{{1, 0}, {0, 1}}}));
  EXPECT_DOUBLE_EQ(std::sqrt(3), h1_norm(grid_view, g, h, hd));
}
