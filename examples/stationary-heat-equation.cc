// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#include "config.h"

#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>


using namespace Dune;
using namespace Dune::GDT;


// some global defines
using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using M = XT::LA::IstlRowMajorSparseMatrix<double>;
using V = XT::LA::IstlDenseVector<double>;


int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    const double diffusion = 1;
    const XT::Functions::GenericFunction<d> source(3, [](const auto& x, const auto& /*param*/) {
      return M_PI_2 * M_PI * std::cos(M_PI_2 * x[0]) * std::cos(M_PI_2 * x[1]);
    });
    const XT::Functions::GridFunction<E> exact_solution(XT::Functions::GenericFunction<d>(
        3,
        /*evaluate=*/
        [](const auto& x, const auto& /*param*/) { return std::cos(M_PI_2 * x[0]) * std::cos(M_PI_2 * x[1]); },
        /*name=*/"exact_solution",
        /*parameter_type=*/{},
        /*jacobian=*/
        [](const auto& x, const auto& /*param*/) {
          const auto pre = -0.5 * M_PI;
          const auto x_arg = M_PI_2 * x[0];
          const auto y_arg = M_PI_2 * x[1];
          FieldMatrix<double, 1, d> result;
          result[0] = {pre * std::sin(x_arg) * std::cos(y_arg), pre * std::cos(x_arg) * std::sin(y_arg)};
          return result;
        }));

    auto grid = XT::Grid::make_cube_grid<G>(/*lower_left=*/-1., /*upper_right=*/1., /*num_elements=*/128);
    auto grid_view = grid.leaf_view();

    XT::Grid::AllDirichletBoundaryInfo<I> boundary_info;

    auto space = make_continuous_lagrange_space(grid_view, /*polorder=*/1);

    auto lhs_op = make_matrix_operator<M>(space, Stencil::element);
    lhs_op.append(LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(diffusion)));

    auto rhs_func = make_vector_functional<V>(space);
    rhs_func.append(LocalElementIntegralFunctional<E>(LocalProductIntegrand<E>().with_ansatz(source)));

    auto dirichlet_constraints = make_dirichlet_constraints(space, boundary_info);

    auto walker = XT::Grid::make_walker(grid_view);
    walker.append(lhs_op);
    walker.append(rhs_func);
    walker.append(dirichlet_constraints);
    walker.walk(/*thread_parallel=*/true);

    dirichlet_constraints.apply(lhs_op.matrix(), rhs_func.vector());
    auto solution = make_discrete_function<V>(space);
    XT::LA::make_solver(lhs_op.matrix()).apply(rhs_func.vector(), solution.dofs().vector());

    solution.visualize("solution");

    const auto error = solution - exact_solution;

    auto h1_prod = make_bilinear_form(grid_view, error, error);
    h1_prod += LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(diffusion));

    auto l2_prod = make_bilinear_form(grid_view, error, error);
    l2_prod += LocalElementIntegralBilinearForm<E>(LocalProductIntegrand<E>());

    walker.append(h1_prod);
    walker.append(l2_prod);
    walker.walk(/*thread_parallel=*/true);

    logger.info() << "error in H^1 semi-norm: " << std::sqrt(h1_prod.result()) << std::endl;
    logger.info() << "error in L^2 norm:      " << std::sqrt(l2_prod.result()) << std::endl;

  } catch (Exception& e) {
    std::cerr << "\nDUNE reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "\nstl reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occured!" << std::endl;
    return EXIT_FAILURE;
  } // try
  return EXIT_SUCCESS;
} // ... main(...)
