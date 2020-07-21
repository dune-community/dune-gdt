// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include "config.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <utility>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/laplace-ipdg.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/laplace-ipdg-flux-reconstruction.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/oswald-interpolation.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/test/stationary-heat-equation/ESV2007.hh>
#include <dune/gdt/tools/adaptation-helper.hh>
#include <dune/gdt/tools/doerfler-marking.hh>

using namespace Dune;
using namespace Dune::GDT;


// some global defines
using G = ALU_2D_SIMPLEX_CONFORMING;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using F = double;
using M = XT::LA::IstlRowMajorSparseMatrix<F>;
using V = XT::LA::IstlDenseVector<F>;
static const LocalEllipticIpdgIntegrands::Method ipdg = LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor;


std::pair<V, F> compute_local_indicators(const DiscreteFunction<V, GV>& u_h,
                                         const GV& grid_view,
                                         const XT::Functions::GridFunctionInterface<E, d, d>& dt,
                                         const XT::Functions::GridFunctionInterface<E>& f,
                                         const XT::Grid::BoundaryInfo<I>& boundary_info,
                                         const int over_integrate = 3)
{
  const auto& dg_space = u_h.space();
  const XT::Functions::ConstantFunction<d> one(1.);
  const auto& df = one.template as_grid_function<E>();
  auto diffusion = df * dt;
  // oswald interpolation
  auto oswald_interpolation_operator =
      make_oswald_interpolation_operator<M>(grid_view, dg_space, dg_space, boundary_info);
  oswald_interpolation_operator.assemble(DXTC_CONFIG_GET("parallel", true));
  auto u = oswald_interpolation_operator.apply(u_h);
  // flux reconstruction
  auto rt_space = make_raviart_thomas_space(grid_view, std::max(dg_space.max_polorder() - 1, 0));
  auto reconstruction_op = make_ipdg_flux_reconstruction_operator<M, ipdg>(grid_view, dg_space, rt_space, df, dt);
  auto t_h = reconstruction_op.apply(u_h);
  // the best index set for a grid view of arbitrary elements is a scalar FV space ...
  auto fv_space = GDT::make_finite_volume_space<1>(grid_view);
  // ... and the best thread-safe way to associate data with grid elements is a corresponding discrete function
  auto indicators = GDT::make_discrete_function<V>(fv_space);
  auto walker = XT::Grid::make_walker(grid_view);
  walker.append([](/*prepare nothing*/) {},
                [&](const auto& element) {
                  // prepare data functions
                  auto u_h_el = u_h.local_function();
                  auto u_el = u.local_function();
                  auto t_h_el = t_h.local_function();
                  auto div_t_h_el = XT::Functions::divergence(*t_h_el);
                  auto f_el = f.local_function();
                  auto d_el = diffusion.local_function();
                  u_h_el->bind(element);
                  u_el->bind(element);
                  t_h_el->bind(element);
                  div_t_h_el.bind(element);
                  f_el->bind(element);
                  d_el->bind(element);
                  // eta_NC
                  const double eta_NC_element_2 =
                      LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(dt), over_integrate)
                          .apply2(*u_h_el - *u_el, *u_h_el - *u_el)[0][0];
                  // eta_R
                  // - approximate minimum eigenvalue of the diffusion over the element (evaluate at some points)
                  double min_EV = std::numeric_limits<double>::max();
                  for (auto&& quadrature_point :
                       QuadratureRules<double, d>::rule(element.type(), d_el->order() + over_integrate)) {
                    auto diff = d_el->evaluate(quadrature_point.position());
                    auto eigen_solver =
                        XT::LA::make_eigen_solver(diff,
                                                  {{"type", XT::LA::EigenSolverOptions<decltype(diff)>::types().at(0)},
                                                   {"assert_positive_eigenvalues", "1e-15"}});
                    min_EV = std::min(min_EV, eigen_solver.min_eigenvalues(1).at(0));
                  }
                  DUNE_THROW_IF(!(min_EV > 0.),
                                Exceptions::integrand_error,
                                "The minimum eigenvalue of a positiv definite matrix must not be negative!"
                                    << "\n\nmin_EV = " << min_EV);
                  const auto C_P = 1. / (M_PIl * M_PIl); // Poincare constant (known for simplices/cubes)
                  const auto h = XT::Grid::diameter(element);
                  auto L2_norm_2 = LocalElementIntegralBilinearForm<E>(LocalProductIntegrand<E>(), over_integrate)
                                       .apply2(*f_el - div_t_h_el, *f_el - div_t_h_el)[0][0];
                  const double eta_R_element_2 = (C_P * h * h * L2_norm_2) / min_EV;
                  // eta_DF
                  const double eta_DF_element_2 = XT::Grid::element_integral(
                      element,
                      [&](const auto& xx) {
                        const auto diff = d_el->evaluate(xx);
                        const auto diff_inv = XT::LA::invert_matrix(d_el->evaluate(xx));
                        const auto pressure_grad = u_h_el->jacobian(xx)[0];
                        const auto t_val = t_h_el->evaluate(xx);
                        auto difference = diff * pressure_grad + t_val;
                        return (diff_inv * difference) * difference;
                      },
                      std::max(d_el->order() + std::max(u_h_el->order() - 1, 0), t_h_el->order()) + over_integrate);
                  // compute indicators and estimator
                  auto local_indicator = indicators.local_discrete_function(element);
                  local_indicator->dofs()[0] = std::sqrt(
                      eta_NC_element_2 + std::pow(std::sqrt(eta_R_element_2) + std::sqrt(eta_DF_element_2), 2));
                },
                [](/*finalize nothing*/) {});
  walker.walk(/*parallel=*/true);
  return {indicators.dofs().vector(), indicators.dofs().vector().l2_norm()};
} // ... compute_local_indicators(...)


void mark_elements(
    G& grid, const GV& grid_view, const V& local_indicators, const double& refine_factor, const double& coarsen_factor)
{
  // compute dörfler marking
  auto marked_elements = GDT::doerfler_marking(local_indicators, refine_factor, coarsen_factor);
  const auto& elements_to_be_coarsened = marked_elements.first;
  const auto& elements_to_be_refined = marked_elements.second;
  // mark elements, as above, use a scalar FV space as index mapper
  auto fv_space = GDT::make_finite_volume_space<1>(grid_view);
  size_t corsend_elements = 0;
  size_t refined_elements = 0;
  for (auto&& element : elements(grid_view)) {
    const size_t index = fv_space.mapper().global_indices(element)[0];
    bool coarsened = false;
    if (std::find(elements_to_be_refined.begin(), elements_to_be_refined.end(), index)
        != elements_to_be_refined.end()) {
      grid.mark(/*refine*/ 2, element);
      refined_elements += 1;
      coarsened = false;
    } else if (std::find(elements_to_be_coarsened.begin(), elements_to_be_coarsened.end(), index)
               != elements_to_be_coarsened.end()) {
      grid.mark(/*coarsen*/ -1, element);
      coarsened = true;
    }
    if (coarsened)
      ++corsend_elements;
  }
  auto logger = XT::Common::TimedLogger().get("mark_elements");
  logger.info() << "marked " << corsend_elements << "/" << fv_space.mapper().size() << " elements for coarsening and "
                << refined_elements << "/" << fv_space.mapper().size() << " elements for refinement" << std::endl;
} // ... mark_elements(...)


int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // problem
    const Test::ESV2007DiffusionProblem<GV> problem;

    // grid
    auto grid_ptr = problem.make_initial_grid().grid_ptr();
    auto& grid = *grid_ptr;
    auto grid_view = grid.leafGridView();

    // space
    auto dg_space = make_discontinuous_lagrange_space(grid_view, 1);
    auto current_solution = make_discrete_function<V>(dg_space);

    // the main adaptation loop
    auto helper = make_adaptation_helper(grid, dg_space, current_solution);
    const double tolerance = DXTC_CONFIG_GET("tolerance", 1e-1);
    size_t counter = 0;
    while (true) {
      logger.info() << "step " << counter << ", space has " << dg_space.mapper().size() << " DoFs" << std::endl;

      // assemble
      auto lhs_op = make_matrix_operator<M>(dg_space, Stencil::element_and_intersection);
      lhs_op.append(LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(problem.diffusion)));
      lhs_op.append(LocalCouplingIntersectionIntegralBilinearForm<I>(
                        LocalLaplaceIPDGIntegrands::InnerCoupling<I>(/*SIPDG=*/1, problem.diffusion)
                        + LocalIPDGIntegrands::InnerPenalty<I>(/*penalty=*/16, problem.diffusion)),
                    {},
                    XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
      lhs_op.append(
          LocalIntersectionIntegralBilinearForm<I>(
              LocalIPDGIntegrands::BoundaryPenalty<I>(/*penalty=*/16, problem.diffusion)
              + LocalLaplaceIPDGIntegrands::DirichletCoupling<I>(/*SIPDG=*/1, problem.diffusion)),
          {},
          XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(problem.boundary_info, new XT::Grid::DirichletBoundary()));
      auto rhs_func = make_vector_functional<V>(dg_space);
      rhs_func.append(LocalElementIntegralFunctional<E>(LocalProductIntegrand<E>().with_ansatz(problem.force)));
      // ... add Dirichlet here
      // (if we add something here, the oswald interpolation needs to be adapted accordingly!)
      // ... add Neumann here
      // assemble everything in one grid walk
      lhs_op.append(rhs_func);
      lhs_op.assemble(DXTC_CONFIG_GET("parallel", true));

      // solve, use previous solution as initial guess (if the solver supports it)
      XT::LA::make_solver(lhs_op.matrix()).apply(rhs_func.vector(), current_solution.dofs().vector());
      current_solution.visualize("solution_" + XT::Common::to_string(counter));

      // compute local indicators and estimate
      const auto estimates = compute_local_indicators(current_solution,
                                                      grid_view,
                                                      problem.diffusion.template as_grid_function<E>(),
                                                      problem.force.template as_grid_function<E>(),
                                                      problem.boundary_info);
      const auto indicators = std::move(estimates.first);
      const auto estimate = estimates.second;
      logger.info() << "  estimated error: " << estimate << std::endl;
      if (estimate < tolerance) {
        logger.info() << "target accuracy reached, terminating!" << std::endl;
        break;
      }

      // use these as indicators for grid refinement
      mark_elements(grid,
                    grid_view,
                    indicators,
                    DXTC_CONFIG_GET("mark.refine_factor", 0.25),
                    DXTC_CONFIG_GET("mark.coarsen_factor", 0.01));

      // adapt the grid (restricts/prolongs the solution)
      helper.pre_adapt();
      helper.adapt();
      helper.post_adapt();

      ++counter;
    }

    double min_h = std::numeric_limits<double>::max();
    for (auto&& element : elements(grid_view))
      min_h = std::min(min_h, XT::Grid::diameter(element));

    logger.info() << "\nA uniformly refined grid of similar accuracy would have roughly required a grid width of "
                  << min_h << std::endl;

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
