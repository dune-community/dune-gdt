// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include <config.h>

#include <cmath>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/timedlogging.hh>


#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>
#include <dune/xt/functions/sliced.hh>
#include <dune/xt/functions/transformed.hh>

#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/sap.hh>

using namespace Dune;
using namespace Dune::GDT;


using G = YASP_1D_EQUIDISTANT_OFFSET;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t m = d + 1;

using DomainType = XT::Common::FieldVector<D, d>;
using RangeType = XT::Common::FieldVector<D, m>;


int main(int argc, char** argv)
{
  try {
// init MPI, global config and logger
#if HAVE_DUNE_FEM
    Fem::MPIManager::initialize(argc, argv);
#else
    MPIHelper::instance(argc, argv);
#endif
    DXTC_CONFIG.read_command_line(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG.get("logger.max_info_level", 999),
                                     DXTC_CONFIG.get("logger.max_debug_level", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // define model
    const R p_0 = DXTC_CONFIG.get("problem.p_0", 1.01325e5);
    const R rho_0 = DXTC_CONFIG.get("problem.rho_0", 9.982e2);
    const R c_f = DXTC_CONFIG.get("problem.c_f", 1.482e3);
    const R nu_s = DXTC_CONFIG.get("problem.nu_s", 1.002e-3);
    const R nu_b = DXTC_CONFIG.get("problem.nu_b", 2.87e-3);
    const R Delta_p = DXTC_CONFIG.get("problem.Delta_p", 2.e1);
    const R Delta_rho_tilde = DXTC_CONFIG.get("problem.Delta_rho_tilde", Delta_p / (rho_0 * c_f * c_f));
    const R f = DXTC_CONFIG.get("problem.f", 1.e6);
    const R L = DXTC_CONFIG.get("problem.L", 1.e-6);
    const R bv_frequency_factor = DXTC_CONFIG.get("problem.f_times_L_over_c_f", f * L / c_f);
    SelfacoustophoreticParticleTools<d, R> tools(p_0, rho_0, c_f, nu_s, nu_b);

    // create grid
    auto grid = XT::Grid::make_cube_grid<G>(
        -0.5 * DXTC_CONFIG.get("problem.l_1", 5e-6) / L, 0., DXTC_CONFIG.get("grid.num_elements", 128));

    auto grid_layer = grid.leaf_view();
    logger.info() << "grid has " << grid_layer.indexSet().size(0) << " elements" << std::endl;

    // define boundary types
    using GL = std::decay_t<decltype(grid_layer)>;
    using I = XT::Grid::extract_intersection_t<GL>;
    XT::Grid::NormalBasedBoundaryInfo<I> boundary_info;
    boundary_info.register_new_normal(-1, new XT::Grid::InflowOutflowBoundary());
    boundary_info.register_new_normal(1, new XT::Grid::ImpermeableBoundary());
    visualize_grid(grid_layer, grid_layer.indexSet(), boundary_info, "grid");

    // define boundary values
    const XT::Functions::GlobalLambdaFunction<E, D, d, R, m> acoustic_wave(
        [&](const auto& /*xx*/, const auto& param) {
          const R t = param.get("t_").at(0);
          FieldVector<R, m> primitive_variables(0);
          // density
          primitive_variables[0] = 1. + Delta_rho_tilde * std::sin(2. * M_PI * bv_frequency_factor * t);
          // velocity: zero
          return tools.to_conservative(primitive_variables);
        },
        /*order=*/3,
        /*parameter_type=*/{"t_", 1});
    const auto& bv = acoustic_wave;

    const auto heuristic_inflow_outflow_treatment = [&](const auto& intersection,
                                                        const auto& x_intersection,
                                                        const auto& /*f*/,
                                                        const auto& u,
                                                        const auto& param = {}) {
      // determine boundary
      const auto n = intersection.centerUnitOuterNormal();
      if (n[0] > 0)
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "There should not be an inflow/outflow boundary at the particle wall!");
      // left boundary: inflow -> evaluate boundary values
      const auto entity = intersection.inside();
      const auto x_entity = intersection.geometryInInside().global(x_intersection);
      const auto bv_cons = bv.local_function(entity)->evaluate(x_entity, param);
      const auto bv_prim = tools.to_primitive(bv_cons);
      // use given density...
      RangeType v_prim;
      v_prim[0] = bv_prim[0];
      // ... but extrapolated velocity
      v_prim[1] = tools.to_primitive(u)[1];
      return tools.to_conservative(v_prim);
    };
    const auto& inflow_outflow_treatment = heuristic_inflow_outflow_treatment;
    XT::Grid::ApplyOn::CustomBoundaryIntersections<GL> inflow_filter(boundary_info,
                                                                     new XT::Grid::InflowOutflowBoundary());

    // impermeable wall treatment, see [DF2015, p. 415, (8.66 - 8.67)]
    const auto inviscid_mirror_impermeable_wall_treatment = [&](const auto& intersection,
                                                                const auto& x_intersection,
                                                                const auto& /*f*/,
                                                                const auto& u,
                                                                const auto& /*mu*/ = {}) {
      const auto normal = intersection.unitOuterNormal(x_intersection);
      if (normal[0] < 0)
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "There should not be an impermeable wall boundary at the inflow!");
      const auto rho = u[0];
      auto velocity = tools.velocity_from_conservative(u);
      velocity -= normal * 2. * (velocity * normal);
      return tools.to_conservative(XT::Common::hstack(rho, velocity));
    };
    const auto& impermeable_wall_treatment = inviscid_mirror_impermeable_wall_treatment;
    XT::Grid::ApplyOn::CustomBoundaryIntersections<GL> impermeable_wall_filter(boundary_info,
                                                                               new XT::Grid::ImpermeableBoundary());

    using S = FvSpace<GL, R, m>;
    S space(grid_layer);
    logger.info() << "space has " << space.mapper().size() << " DoFs" << std::endl;

    using V = XT::LA::EigenDenseVector<R>;
    using DF = DiscreteFunction<S, V>;
    DF initial_values(space, "solution");

    // define initial values
    const XT::Functions::GlobalLambdaFunction<E, D, d, R, m> constant_initial_values(
        [&](const auto& /*xx*/, const auto& /*param*/) {
          // velocity zero
          FieldVector<R, m> primitive_variables(0);
          // density
          primitive_variables[0] = 1;
          // velocity: zero
          return tools.to_conservative(primitive_variables);
        },
        /*order=*/0);
    const auto& u_0 = constant_initial_values;
    tools.visualize(u_0, grid_layer, "initial_values");

    project(u_0, initial_values);
    tools.visualize(initial_values, grid_layer, "projected_initial_values");

    using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
    const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> inviscid_flux(
        [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*param*/) {
          return tools.flux(conservative_variables);
        },
        {},
        "inviscid_flux",
        [](const auto& /*param*/) { return 3; },
        [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*param*/) {
          return tools.flux_jacobian(conservative_variables);
        });
    const auto& flux = inviscid_flux;

    const auto eigenvalue_decomposition_tolerance =
        DXTC_CONFIG.get("numerical_flux.eigenvalue_decomposition_tolerance", 1e-15);
    const auto vijayasundaram_flux = make_numerical_vijayasundaram_flux(flux, [&](const auto& w, const auto& n) {
      const auto eigenvalues = tools.eigenvalues_flux_jacobi_matrix(w, n);
      const auto eigenvectors = tools.eigenvectors_flux_jacobi_matrix(w, n);
      const auto eigenvectors_inv = tools.eigenvectors_inv_flux_jacobi_matrix(w, n);

#ifndef DUNE_GDT_DISABLE_CHECKS
      const auto identity = XT::LA::eye_matrix<XT::Common::FieldMatrix<R, m, m>>(m, m);
      if ((eigenvectors_inv * eigenvectors - identity).infinity_norm() > eigenvalue_decomposition_tolerance)
        DUNE_THROW(InvalidStateException,
                   "\n\neigenvectors:\n\n"
                       << std::setprecision(17)
                       << eigenvectors
                       << "\n\neigenvectors_inverse:\n\n"
                       << eigenvectors_inv
                       << "\n\n|| eigenvectors_inv * eigenvectors - identity ||_infty = "
                       << (eigenvectors_inv * eigenvectors - identity).infinity_norm());

      const auto eigenvaluematrix = tools.eigenvaluematrix_flux_jacobi_matrix(w, n);
      if (((eigenvectors_inv * (tools.flux_jacobi_matrix(w, n) * eigenvectors)) - eigenvaluematrix).infinity_norm()
          > eigenvalue_decomposition_tolerance)
        DUNE_THROW(InvalidStateException,
                   "\n\neigenvectors:\n\n"
                       << std::setprecision(17)
                       << eigenvectors
                       << "\n\neigenvectors_inverse:\n\n"
                       << eigenvectors_inv
                       << "\n\neigenvalues:"
                       << eigenvalues
                       << "\n\nP:\n\n"
                       << tools.flux_jacobi_matrix(w, n)
                       << "\n\neigenvectors_inv * (P * eigenvectors):\n\n"
                       << eigenvectors_inv * (tools.flux_jacobi_matrix(w, n) * eigenvectors)
                       << "\n\n|| eigenvectors_inv * (P * eigenvectors) - eigenvalues||_infty = "
                       << ((eigenvectors_inv * (tools.flux_jacobi_matrix(w, n) * eigenvectors)) - eigenvaluematrix)
                              .infinity_norm());
#endif // DUNE_GDT_DISABLE_CHECKS
      return std::make_tuple(eigenvalues, eigenvectors, eigenvectors_inv);
    });
    const auto& numerical_flux = vijayasundaram_flux;

    using OpType = GDT::AdvectionFvOperator<DF>;
    OpType advec_op(grid_layer,
                    numerical_flux,
                    /*periodicity_restriction=*/impermeable_wall_filter || inflow_filter);
    // non-periodic boundary treatment
    advec_op.append(inflow_outflow_treatment, inflow_filter.copy(), {"t_", 1});
    advec_op.append(impermeable_wall_treatment, impermeable_wall_filter.copy());

    auto dt = estimate_dt_for_hyperbolic_system(grid_layer, u_0, flux, {1. - Delta_rho_tilde, 1. + Delta_rho_tilde});
    logger.info() << "estimated dt via [CCL1995] is " << dt;
    const auto max_user_dt = DXTC_CONFIG.get("timestepping.max_dt", bv_frequency_factor / 10.);
    if (dt > max_user_dt) {
      dt = max_user_dt;
      logger.info() << ", but user requests " << dt;
    }
    logger.info() << std::endl;

    ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
        advec_op, initial_values, -1.);
    const auto test_dt =
        time_stepper.find_suitable_dt(dt,
                                      10,
                                      DXTC_CONFIG.get("timestepping.max_solution_variation", 1.1)
                                          * std::max(initial_values.vector().sup_norm(), 1. + Delta_rho_tilde),
                                      25);
    if (!test_dt.first)
      DUNE_THROW(InvalidStateException,
                 "Could not determine optimal dt (in particular, the dt computed to match the "
                 "CFL condition did not yield a stable scheme)!");
    if (test_dt.second < dt) {
      DUNE_THROW(InvalidStateException,
                 "The computed dt (to match the CFL condition) does not yield a stable scheme: "
                     << dt
                     << "\n   The following dt seems to work fine: "
                     << test_dt.second);
    } else {
      R t_max = 0;
      const auto max_timesteps = DXTC_CONFIG.get("timestepping.max_steps", -1);
      if (max_timesteps > 0)
        t_max = dt * max_timesteps;
      else
        t_max = DXTC_CONFIG.get("problem.t_max", 10 * bv_frequency_factor);
      time_stepper.solve(t_max,
                         dt,
                         DXTC_CONFIG.get("timestepping.num_saves", -1),
                         false,
                         true,
                         "solution",
                         /*visualizer=*/[&](const auto& u, const auto& filename_prefix, const auto& step) {
                           tools.visualize(u, grid_layer, filename_prefix, XT::Common::to_string(step));
                         });
    }
  } catch (const Dune::Exception& e) {
    std::cerr << "\nDUNE reported: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "\nSTL reported: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return EXIT_FAILURE;
  } // try
  return EXIT_SUCCESS;
} // ... main(...)
