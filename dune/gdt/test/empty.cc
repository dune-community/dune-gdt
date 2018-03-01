// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)

/**
  * This file is intended as a starting point for quick testing.
  */

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/constant.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/weighted-l2.hh>
#include <dune/gdt/spaces/l2/discontinuous-galerkin.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>
#include <dune/gdt/assembler/global.hh>

using namespace Dune;
using namespace Dune::GDT;


GTEST_TEST(empty, main)
{
  using G = YASP_1D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0, 1, 2);
  auto grid_view = grid.leaf_view();
  using GV = decltype(grid_view);
  using E = typename GV::template Codim<0>::Entity;
  ContinuousLagrangeSpace<GV, 1> space(grid_view);

  const XT::Functions::ConstantFunction<d> func(1);
  const LocalElementProductIntegrand<E> product_integrand(func.as_localizable<E>());
  const LocalElementIntegralOperator<E> local_op(product_integrand);

  auto op = make_matrix_operator<XT::LA::CommonDenseMatrix<double>>(grid_view, space);
  op.append(local_op);
  op.assemble();

  auto functional =
      make_l2_volume_vector_functional<XT::LA::CommonDenseVector<double>>(space, func.as_localizable<E>());
  functional.assemble();

  auto assembler = make_global_assembler(space);
  assembler.append(functional);
  assembler.assemble();

  std::cout << "vector = " << functional.vector() << std::endl;
  std::cout << "matrix = \n" << op.matrix() << std::endl;

  XT::LA::CommonDenseVector<double> vector({1, 2, 3});
  auto df = make_discrete_function(space, vector);
  df.visualize("df");
  auto& dofs = df.dofs();
  std::cout << "dofs.vector() = " << dofs.vector() << std::endl;

  auto local_dofs = dofs.localize();
  for (auto&& element : elements(grid_view))
    std::cout << local_dofs.bind(element) << std::endl;

  std::cout << "functional.apply(df) = " << functional.apply(df) << std::endl;
  std::cout << "op.apply(df).dofs().vector() = " << op.apply(df).dofs().vector() << std::endl;
  std::cout << "op.apply_inverse(df).dofs().vector() = " << op.apply_inverse(df).dofs().vector() << std::endl;
}
