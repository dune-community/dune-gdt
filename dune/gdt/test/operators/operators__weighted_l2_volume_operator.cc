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
#include <dune/xt/functions/constant.hh>

#include <dune/gdt/operators/weighted-l2.hh>
#include <dune/gdt/spaces/l2/discontinuous-galerkin.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

using namespace Dune;
using namespace Dune::GDT;


GTEST_TEST(matrix_operator, make_weighted_l2_volume_matrix_operator)
{
  using G = YASP_1D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0, 1, 2);
  auto grid_view = grid.leaf_view();
  using GV = decltype(grid_view);
  using E = typename GV::template Codim<0>::Entity;
  const ContinuousLagrangeSpace<GV, 1> space(grid_view);

  const XT::Functions::ConstantFunction<d> weight(1);
  const auto& localizable_weight = weight.as_localizable<E>();
  const int over_integrate = 0;
  const XT::Common::Parameter param;
  const XT::Grid::ApplyOn::AllElements<GV> filter;

  using MatrixType = XT::LA::CommonDenseMatrix<double>;
  auto matrix = XT::LA::eye_matrix<MatrixType>(space.mapper().size());

  // given a matrix, a grid view and two spaces
  make_weighted_l2_volume_matrix_operator(
      grid_view, space, space, matrix, localizable_weight, over_integrate, param, filter);
  make_weighted_l2_volume_matrix_operator(grid_view, space, space, matrix, localizable_weight, over_integrate, param);
  make_weighted_l2_volume_matrix_operator(grid_view, space, space, matrix, localizable_weight, over_integrate);
  make_weighted_l2_volume_matrix_operator(grid_view, space, space, matrix, localizable_weight);
  // given a matrix, a grid view and a space
  make_weighted_l2_volume_matrix_operator(grid_view, space, matrix, localizable_weight, over_integrate, param, filter);
  make_weighted_l2_volume_matrix_operator(grid_view, space, matrix, localizable_weight, over_integrate, param);
  make_weighted_l2_volume_matrix_operator(grid_view, space, matrix, localizable_weight, over_integrate);
  make_weighted_l2_volume_matrix_operator(grid_view, space, matrix, localizable_weight);
  /// \todo Add the other ones.
}
