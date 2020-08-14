// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/bisect.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/bochner.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/numerical-fluxes/engquist-osher.hh>
#include <dune/gdt/local/numerical-fluxes/generic.hh>
#include <dune/gdt/local/numerical-fluxes/upwind.hh>
#include <dune/gdt/local/numerical-fluxes/vijayasundaram.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/identity.hh>
#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/spaces/bochner.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/tools/euler.hh>
#include <dune/gdt/tools/hyperbolic.hh>

using namespace Dune;
using namespace Dune::GDT;

using V = XT::LA::IstlDenseVector<double>;
using M = XT::LA::IstlRowMajorSparseMatrix<double>;


std::vector<double> time_points_from_vector_array(const XT::LA::ListVectorArray<V>& va)
{
  std::vector<double> time_points;
  time_points.reserve(va.length());
  for (const auto& note : va.notes())
    time_points.emplace_back(note.get("_t").at(0));
  return time_points;
}


template <class GV, size_t m>
void visualize(const SpaceInterface<GV, m>& space,
               const double& T_end,
               XT::LA::ListVectorArray<V>& solution_vectors,
               const std::string& prefix)
{
  const BochnerSpace<GV> bochner_space(space, time_points_from_vector_array(solution_vectors));
  const auto timedomain_solution_function = make_discrete_bochner_function(bochner_space, solution_vectors);
  for (size_t ii = 0; ii < 100; ++ii) {
    const double time = ii * (T_end / 100);
    timedomain_solution_function.evaluate(time).visualize(prefix + "_" + XT::Common::to_string(ii));
  }
}


template <class M, class GV, size_t m>
double estimate_fixed_explicit_dt_for_some_steps(const OperatorInterface<M, GV, m>& op,
                                                 const DiscreteFunction<V, GV, m>& w_0,
                                                 const double& T_end,
                                                 const double max_overshoot = 1.25,
                                                 const int max_steps_to_try = 250)
{
  const auto max_sup_norm = max_overshoot * w_0.dofs().vector().sup_norm();
  return XT::Common::find_largest_by_bisection(
      /*min_dt=*/10 * std::numeric_limits<double>::epsilon(),
      /*max_dt=*/T_end,
      /*success=*/[&](const auto& dt_to_test) {
        try {
          auto w = w_0.dofs().vector();
          const double T_end = max_steps_to_try * dt_to_test;
          double time = 0.;
          // explicit euler
          while (time < T_end + dt_to_test) {
            w -= op.apply(w, {{"_t", {time}}, {"_dt", {dt_to_test}}}) * dt_to_test;
            time += dt_to_test;
            if (w.sup_norm() > max_sup_norm)
              return false;
          }
          return true;
        } catch (...) {
          return false;
        }
      });
} // ... estimate_fixed_explicit_dt_for_some_steps(...)


template <class M, class GV, size_t m>
double estimate_fixed_explicit_dt_to_T_end(const OperatorInterface<M, GV, m>& op,
                                           const DiscreteFunction<V, GV, m>& w_0,
                                           const double& min_dt,
                                           const double& T_end,
                                           const double max_overshoot = 1.25)
{
  const auto max_sup_norm = max_overshoot * w_0.dofs().vector().sup_norm();
  return XT::Common::find_largest_by_bisection(
      /*min_dt=*/min_dt,
      /*max_dt=*/T_end,
      /*success=*/[&](const auto& dt_to_test) {
        try {
          auto w = w_0.dofs().vector();
          double time = 0.;
          // explicit euler
          while (time < T_end + dt_to_test) {
            w -= op.apply(w, {{"_t", {time}}, {"_dt", {dt_to_test}}}) * dt_to_test;
            time += dt_to_test;
            if (w.sup_norm() > max_sup_norm)
              return false;
          }
          return true;
        } catch (...) {
          return false;
        }
      });
} // ... estimate_fixed_explicit_dt_to_T_end(...)


template <class V, class GV, size_t m, class M>
XT::LA::ListVectorArray<V> explicit_euler(const DiscreteFunction<V, GV, m>& initial_values,
                                          const OperatorInterface<M, GV, m>& spatial_op,
                                          const double T_end,
                                          const double dt)
{
  // initial values
  XT::LA::ListVectorArray<V> solution(
      spatial_op.source_space().mapper().size(), /*length=*/0, /*reserve=*/std::ceil(T_end / (dt)));
  solution.append(initial_values.dofs().vector(), {"_t", 0.});
  // timestepping
  double time = 0.;
  while (time < T_end + dt) {
    const auto& u_n = solution.back().vector();
    auto u_n_plus_one = u_n - spatial_op.apply(u_n, {{"_t", {time}}, {"_dt", {dt}}}) * dt;
    time += dt;
    solution.append(std::move(u_n_plus_one), {"_t", time});
  }
  return solution;
} // ... explicit_euler(...)


template <class V, class GV, size_t m, class M>
XT::LA::ListVectorArray<V> implicit_euler(const DiscreteFunction<V, GV, m>& initial_values,
                                          const OperatorInterface<M, GV, m>& spatial_op,
                                          const double T_end,
                                          const double dt)
{
  // some preparations
  auto id = make_identity_operator(spatial_op);
  V zero(spatial_op.range_space().mapper().size(), 0.);
  // initial values
  XT::LA::ListVectorArray<V> solution(
      spatial_op.source_space().mapper().size(), /*length=*/0, /*reserve=*/std::ceil(T_end / (dt)));
  solution.append(initial_values.dofs().vector(), {"_t", 0.});
  // timestepping
  double time = 0.;
  while (time < T_end + dt) {
    time += dt;
    const auto& u_n = solution.back().vector();
    auto residual_op = (id - u_n) / dt + spatial_op;
    auto u_n_plus_one = residual_op.apply_inverse(zero);
    solution.append(std::move(u_n_plus_one), {"_t", time});
  }
  return solution;
} // ... implicit_euler(...)


GTEST_TEST(MPI201902TalkExamples, instationary_heat_equation)
{
  using G = ONED_1D;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 1024);
  auto grid_view = grid.leaf_view();
  using GV = decltype(grid_view);
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  const XT::Grid::AllDirichletBoundaryInfo<I> boundary_info;

  auto cg_space = make_continuous_lagrange_space(grid_view, 1);
  auto dirichlet_constraints = make_dirichlet_constraints(cg_space, boundary_info);
  auto spatial_op = make_matrix_operator<M>(cg_space, Stencil::element);
  spatial_op.append(LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(1.)));
  spatial_op.append(dirichlet_constraints);
  auto l2_op = make_matrix_operator<M>(cg_space, Stencil::element);
  l2_op.append(LocalElementIntegralBilinearForm<E>(LocalProductIntegrand<E>(1.)));
  spatial_op.append(l2_op);
  spatial_op.assemble(/*SMP=*/true);
  dirichlet_constraints.apply(spatial_op.matrix());
  dirichlet_constraints.apply(l2_op.matrix());

  // generic function to do the time loop
  auto time_loop = [&](const double& T_end, const double& dt) {
    XT::LA::ListVectorArray<V> solution(cg_space.mapper().size(), /*length=*/0, /*reserve=*/std::ceil(T_end / (dt)));
    // initial values
    solution.append(default_interpolation<V>(
                        0,
                        [](const auto& xx, const auto& /*param*/) {
                          if (0.25 <= xx[0] && xx[0] <= 0.5)
                            return 1.;
                          else
                            return 0.;
                        },
                        cg_space)
                        .dofs()
                        .vector(),
                    {"_t", 0.});
    // timestepping
    auto timestepping_mat = spatial_op.matrix();
    timestepping_mat *= dt;
    timestepping_mat += l2_op.matrix();
    double time = 0.;
    while (time < T_end + dt) {
      time += dt;
      const auto& u_n = solution.back().vector();
      auto u_n_plus_one = u_n.copy();
      XT::LA::make_solver(timestepping_mat).apply(l2_op.matrix() * u_n, u_n_plus_one);
      solution.append(std::move(u_n_plus_one), {"_t", time});
    }
    // visualization
    const BochnerSpace<GV> bochner_space(cg_space, time_points_from_vector_array(solution));
    const auto timdomain_solution = make_discrete_bochner_function(bochner_space, solution);
    for (size_t ii = 0; ii < 100; ++ii) {
      time = ii * (T_end / 100);
      timdomain_solution.evaluate(time).visualize(XT::Common::Test::get_unique_test_name() + "__dt"
                                                  + XT::Common::to_string(dt) + "_solution_"
                                                  + XT::Common::to_string(ii));
    }
  }; // ... time_loop(...)

  // try some different dts
  time_loop(/*T_end=*/1.0, /*dt=*/0.1);
  time_loop(/*T_end=*/1.0, /*dt=*/0.5);
} // GTEST_TEST(MPI201902TalkExamples, instationary_heat_equation)


GTEST_TEST(MPI201902TalkExamples, linear_transport)
{
  using G = ONED_1D;
  static const size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 1024);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  auto V_h_0 = make_finite_volume_space(grid_view);

  const XT::Functions::GenericFunction<1, d, 1> f(
      1,
      [&](const auto& w, const auto& /*param*/) { return 1. * w; },
      "transport_to_the_right",
      {},
      [&](const auto& /*w*/, const auto& /*param*/) { return 1.; });
  const NumericalUpwindFlux<I, d, 1> g(f);
  auto L_h = make_advection_fv_operator<M>(grid_view, g, V_h_0, V_h_0);

  auto w_0 = default_interpolation<V>(
      0,
      [](const auto& xx, const auto& /*param*/) {
        if (0.25 <= xx[0] && xx[0] <= 0.5)
          return 1.;
        else
          return 0.;
      },
      V_h_0);

  // explicit, the right choice
  const double T_end = 1.;
  double dt = 1. / 1024; // <- dt = h, only true for linear transport with velocity 1
  auto w_h_reference = explicit_euler(w_0, L_h, T_end, dt);
  visualize(V_h_0, T_end, w_h_reference, XT::Common::Test::get_unique_test_name() + "__reference_solution");

  // explicit, the wrong choice
  auto w_h_explicit_large_dt = explicit_euler(w_0, L_h, T_end, 10 * dt);
  visualize(V_h_0, T_end, w_h_explicit_large_dt, XT::Common::Test::get_unique_test_name() + "__larger_dt_solution");

  // implicit
  auto w_h_implicit = implicit_euler(w_0, L_h, T_end, 10 * dt);
  visualize(V_h_0, T_end, w_h_implicit, XT::Common::Test::get_unique_test_name() + "__implicit_solution");
} // GTEST_TEST(MPI201902TalkExamples, linear_transport)


GTEST_TEST(MPI201902TalkExamples, burgers)
{
  using G = ONED_1D;
  static const size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 1024);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  auto V_h_0 = make_finite_volume_space(grid_view);

  const XT::Functions::GenericFunction<1, d, 1> f(
      2,
      [&](const auto& w, const auto& /*param*/) { return 0.5 * w * w; },
      "burgers",
      {},
      [&](const auto& w, const auto& /*param*/) { return w; });
  const NumericalUpwindFlux<I, d, 1> g(f);
  auto L_h = make_advection_fv_operator<M>(grid_view, g, V_h_0, V_h_0);

  auto w_0 = default_interpolation<V>(
      3,
      [&](const auto& xx, const auto& /*mu*/) {
        return std::exp(-std::pow(xx[0] - 0.33, 2) / (2 * std::pow(0.075, 2)));
      },
      V_h_0);

  // explicit
  const double T_end = 1.;
  const double dt = estimate_dt_for_hyperbolic_system(grid_view, w_0, f);
  auto w_h_reference = explicit_euler(w_0, L_h, T_end, dt);
  visualize(V_h_0, T_end, w_h_reference, XT::Common::Test::get_unique_test_name() + "__reference_solution");

  // implicit
  auto w_h_implicit = implicit_euler(w_0, L_h, T_end, 10 * dt);
  visualize(V_h_0, T_end, w_h_implicit, XT::Common::Test::get_unique_test_name() + "__implicit_solution");
} // GTEST_TEST(MPI201902TalkExamples, burgers)


GTEST_TEST(MPI201902TalkExamples, linear_transport__central_differences)
{
  using G = ONED_1D;
  static const size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 16);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  auto V_h_0 = make_finite_volume_space(grid_view);

  const XT::Functions::GenericFunction<1, d, 1> f(
      1,
      [&](const auto& w, const auto& /*param*/) { return 1. * w; },
      "transport_to_the_right",
      {},
      [&](const auto& /*w*/, const auto& /*param*/) { return 1.; });
  const GenericNumericalFlux<I, d, 1> g(
      f, [&](const auto&, const auto&, const auto& w_minus, const auto& w_plus, const auto& n, const auto& /*param*/) {
        return 0.5 * (f.evaluate(w_minus) + f.evaluate(w_plus)) * n;
      });
  auto L_h = make_advection_fv_operator<M>(grid_view, g, V_h_0, V_h_0);

  auto w_0 = default_interpolation<V>(
      0,
      [](const auto& xx, const auto& /*param*/) {
        if (0.25 <= xx[0] && xx[0] <= 0.5)
          return 1.;
        else
          return 0.;
      },
      V_h_0);

  // explicit, the right choice
  const double T_end = 1.;
  double dt = 1. / 16; // <- dt = h, only true for linear transport with velocity 1
  auto w_h = explicit_euler(w_0, L_h, T_end, dt);
  visualize(V_h_0, T_end, w_h, XT::Common::Test::get_unique_test_name() + "_solution");
} // GTEST_TEST(MPI201902TalkExamples, linear_transport__central_differences)


GTEST_TEST(MPI201902TalkExamples, 2d_euler)
{
  using G = YASP_2D_EQUIDISTANT_OFFSET;
  static const size_t d = G::dimension;
  using DomainType = FieldVector<double, d>;
  static const size_t m = EulerTools<d>::m;
  const EulerTools<d> euler_tools(/*gamma=*/1.4);
  const XT::Functions::GenericFunction<m, d, m> f(
      euler_tools.flux_order(),
      [&](const auto& w, const auto& /*param*/) { return euler_tools.flux(w); },
      "euler_flux",
      {},
      [&](const auto& w, const auto& /*param*/) { return euler_tools.flux_jacobian(w); });

  auto grid = XT::Grid::make_cube_grid<G>(-1., 1., 128);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  auto V_h_0 = make_finite_volume_space<m>(grid_view);

  const NumericalVijayasundaramFlux<I, d, m> g(
      f,
      /*flux_eigen_decomposition=*/[&](const auto& /*f*/, const auto& w, const auto& n, const auto&
                                       /*param*/) {
        return std::make_tuple(euler_tools.eigenvalues_flux_jacobian(w, n),
                               euler_tools.eigenvectors_flux_jacobian(w, n),
                               euler_tools.eigenvectors_inv_flux_jacobian(w, n));
      });
  auto L_h = make_advection_fv_operator<M>(grid_view, g, V_h_0, V_h_0);

  auto w_0 = default_interpolation<V>(
      0,
      [&](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
          return euler_tools.conservative(/*density=*/4., /*velocity=*/0., /*pressure=*/1.6);
        else
          return euler_tools.conservative(/*density=*/1., /*velocity=*/0., /*pressure=*/0.4);
      },
      V_h_0);

  const double T_end = 1.;
  const double dt = estimate_dt_for_hyperbolic_system(grid_view, w_0, f);
  auto w_h = explicit_euler(w_0, L_h, T_end, dt);
  // visualize
  const BochnerSpace<GV, m> bochner_space(V_h_0, time_points_from_vector_array(w_h));
  const auto timedomain_solution_function = make_discrete_bochner_function(bochner_space, w_h);
  for (size_t ii = 0; ii < 100; ++ii) {
    const double time = ii * (T_end / 100.);
    const auto u_t = timedomain_solution_function.evaluate(time);
    euler_tools.visualize(
        u_t, u_t.space().grid_view(), XT::Common::Test::get_unique_test_name() + "_", XT::Common::to_string(ii));
  }
} // GTEST_TEST(MPI201902TalkExamples, 2d_euler)


GTEST_TEST(MPI201902TalkExamples, burgers_p1_unstable)
{
  using G = ONED_1D;
  static const size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 128);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  const DiscontinuousLagrangeSpace<GV> V_h_1(grid_view, 1);

  const XT::Functions::GenericFunction<1, d, 1> f(
      2,
      [&](const auto& w, const auto& /*param*/) { return 0.5 * w * w; },
      "burgers",
      {},
      [&](const auto& w, const auto& /*param*/) { return w; });
  const NumericalUpwindFlux<I, d, 1> g(f);
  // unstable DG operator
  const AdvectionDgOperator<M, GV> L_h(grid_view, g, V_h_1, V_h_1, XT::Grid::ApplyOn::NoIntersections<GV>(), 0., 0.);

  auto w_0 = default_interpolation<V>(
      3,
      [&](const auto& xx, const auto& /*mu*/) {
        return std::exp(-std::pow(xx[0] - 0.33, 2) / (2 * std::pow(0.075, 2)));
      },
      V_h_1);

  const double T_end = 1.;
  const double fv_dt = estimate_dt_for_hyperbolic_system(grid_view, w_0, f);
  auto w_h = explicit_euler(w_0, L_h, T_end, 0.1 * fv_dt);
  visualize(V_h_1, T_end, w_h, XT::Common::Test::get_unique_test_name() + "_solution");
} // GTEST_TEST(MPI201902TalkExamples, burgers_p1_unstable)


GTEST_TEST(MPI201902TalkExamples, burgers_shock_capturing)
{
  using G = ONED_1D;
  static const size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0., 1., 16u);
  auto grid_view = XT::Grid::make_periodic_grid_view(grid.leaf_view());
  using GV = decltype(grid_view);
  using I = XT::Grid::extract_intersection_t<GV>;

  const XT::Functions::GenericFunction<1, d, 1> f(
      2,
      [&](const auto& w, const auto& /*param*/) { return 0.5 * w * w; },
      "burgers",
      {},
      [&](const auto& w, const auto& /*param*/) { return w; });
  const NumericalEngquistOsherFlux<I, d, 1> g(f);

  auto perform_simulation = [&](const int p, const std::string& prefix, const double& CFL_factor = 0.99) {
    const DiscontinuousLagrangeSpace<GV> V_h_p(grid_view, p);

    // stabilized DG operator
    const AdvectionDgOperator<M, GV> L_h(grid_view, g, V_h_p, V_h_p);

    auto w_0 = default_interpolation<V>(
        3,
        [&](const auto& xx, const auto& /*mu*/) {
          return std::exp(-std::pow(xx[0] - 0.33, 2) / (2 * std::pow(0.075, 2)));
        },
        V_h_p);

    const double T_end = 1.;
    const double fv_dt = estimate_dt_for_hyperbolic_system(grid_view, w_0, f);
    const double dt = estimate_fixed_explicit_dt_to_T_end(L_h, w_0, 1e-4, T_end);
    std::cout << "p" << p << ": CFL = " << dt / fv_dt << ", using " << CFL_factor * (dt / fv_dt) << " on request"
              << std::endl;
    auto w_h = explicit_euler(w_0, L_h, T_end, CFL_factor * dt);
    visualize(V_h_p, T_end, w_h, XT::Common::Test::get_unique_test_name() + "_" + prefix + "_solution");
  };

  perform_simulation(1, "p1");
  perform_simulation(2, "p2");
  perform_simulation(3, "p3");
} // GTEST_TEST(MPI201902TalkExamples, burgers_shock_capturing)
