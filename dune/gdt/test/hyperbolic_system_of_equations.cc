// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
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

#include <algorithm>
#include <cmath>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>
#include <dune/xt/functions/sliced.hh>
#include <dune/xt/functions/transformed.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/lambda.hh>
#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>

using namespace Dune;
using namespace Dune::GDT;


template <class V>
typename std::enable_if<XT::Common::is_vector<V>::value, void>::type check_values(const V& vec)
{
  for (size_t ii = 0; ii < vec.size(); ++ii)
    if (XT::Common::isnan(vec[ii]) || XT::Common::isinf(vec[ii]))
      DUNE_THROW(InvalidStateException, vec);
}

template <class M>
typename std::enable_if<XT::Common::is_matrix<M>::value, void>::type check_values(const M& mat)
{
  using MM = XT::Common::MatrixAbstraction<M>;
  for (size_t ii = 0; ii < MM::rows(mat); ++ii)
    for (size_t jj = 0; jj < MM::cols(mat); ++jj)
      if (XT::Common::isnan(MM::get_entry(mat, ii, jj)) || XT::Common::isinf(MM::get_entry(mat, ii, jj)))
        DUNE_THROW(InvalidStateException, mat);
}


// using G = YASP_1D_EQUIDISTANT_OFFSET;
// using G = YASP_2D_EQUIDISTANT_OFFSET;
using G = ALU_2D_SIMPLEX_CONFORMING;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t m = d + 2;

using DomainType = XT::Common::FieldVector<D, d>;


GTEST_TEST(empty, main)
{
  auto grid =
      XT::Grid::make_cube_grid<G>(DomainType(-1.), DomainType(1.), XT::Common::FieldVector<unsigned int, d>(128));
  grid.global_refine(1);

  auto leaf_layer = grid.leaf_view();
  std::cout << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;
  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  auto& grid_layer = periodic_leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;

  //  using I = XT::Grid::extract_intersection_t<GL>;
  //  XT::Grid::NormalBasedBoundaryInfo<I> boundary_info;
  //  boundary_info.register_new_normal({-1., 0.}, new XT::Grid::InflowOutflowBoundary());
  //  boundary_info.register_new_normal({1., 0.}, new XT::Grid::InflowOutflowBoundary());

  const double gamma = 1.4; // air or water at roughly 20 deg Cels.

  const auto to_primitive = [&](const FieldVector<R, m>& conservative_variables) {
    // extract
    const auto& rho = conservative_variables[0];
    DomainType v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[ii + 1] / rho;
    const auto& e = conservative_variables[m - 1];
    // convert
    FieldVector<R, m> primitive_variables;
    // * density
    primitive_variables[0] = rho;
    // * velocity
    for (size_t ii = 0; ii < d; ++ii)
      primitive_variables[ii + 1] = v[ii];
    // * pressure
    primitive_variables[m - 1] = (gamma - 1.) * (e - 0.5 * rho * v.two_norm2());
    return primitive_variables;
  };
  const auto to_conservative = [&](const FieldVector<R, m>& primitive_variables) {
    // extract
    const auto& rho = primitive_variables[0];
    DomainType v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = primitive_variables[ii + 1];
    const auto& p = primitive_variables[m - 1];
    // convert
    FieldVector<R, m> conservative_variables;
    // * density
    conservative_variables[0] = rho;
    // * density times velocity component
    for (size_t ii = 0; ii < d; ++ii)
      conservative_variables[1 + ii] = rho * v[ii];
    // * energy
    conservative_variables[m - 1] = p / (gamma - 1.) + 0.5 * rho * v.two_norm2();
    return conservative_variables;
  };

  const auto visualizer = [&](const auto& u_conservative, const std::string& filename_prefix, const auto step) {
    XT::Functions::make_sliced_function<1>(u_conservative, {0}, "density")
        .visualize(grid_layer, filename_prefix + "_density_" + XT::Common::to_string(step), /*subsampling=*/false);
    XT::Functions::make_sliced_function<d>(u_conservative, {1, 2} /*{1}*/, "density_times_velocity")
        .visualize(grid_layer,
                   filename_prefix + "_density_times_velocity_" + XT::Common::to_string(step),
                   /*subsampling=*/false);
    XT::Functions::make_sliced_function<1>(u_conservative, {3} /*{2}*/, "energy")
        .visualize(grid_layer, filename_prefix + "_energy_" + XT::Common::to_string(step), /*subsampling=*/false);
    const auto u_primitive = XT::Functions::make_transformed_function<m, 1, R>(u_conservative, to_primitive);
    XT::Functions::make_sliced_function<d>(u_primitive, {1, 2} /*{1}*/, "velocity")
        .visualize(grid_layer, filename_prefix + "_velocity_" + XT::Common::to_string(step), /*subsampling=*/false);
    XT::Functions::make_sliced_function<1>(u_primitive, {3} /*{2}*/, "pressure")
        .visualize(grid_layer, filename_prefix + "_pressure_" + XT::Common::to_string(step), /*subsampling=*/false);
  };

  using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
  const U0 indicator(
      [](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::ge(xx, DomainType(0.25)) && XT::Common::FloatCmp::le(xx, DomainType(0.5)))
          return FieldVector<R, m>(2.);
        else
          return FieldVector<R, m>(1.);
      },
      0,
      {},
      "indicator");
  const U0 initial_values_euler( // see [Kröner, 1997, p.394]
      [&](const auto& xx, const auto& /*mu*/) {
        FieldVector<R, m> primitive_variables(0.);
        // density
        if (xx[0] < 0)
          primitive_variables[0] = 4.;
        else
          primitive_variables[0] = 1.;
        // velocity
        for (size_t ii = 0; ii < d; ++ii)
          primitive_variables[1 + ii] = 0.;
        // pressure
        if (xx[0] < 0)
          primitive_variables[m - 1] = 1.6;
        else
          primitive_variables[m - 1] = 0.4;
        return to_conservative(primitive_variables);
      },
      /*order=*/0,
      /*parameter_type=*/{},
      /*name=*/"initial_values_euler");
  const U0 periodic_initial_values_euler(
      [&](const auto& xx, const auto& /*mu*/) {
        FieldVector<R, m> primitive_variables(0.);
        // density
        if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
          primitive_variables[0] = 4.;
        else
          primitive_variables[0] = 1.;
        // velocity
        for (size_t ii = 0; ii < d; ++ii)
          primitive_variables[1 + ii] = 0.;
        // pressure
        if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
          primitive_variables[m - 1] = 1.6;
        else
          primitive_variables[m - 1] = 0.4;
        return to_conservative(primitive_variables);
      },
      /*order=*/0,
      /*parameter_type=*/{},
      /*name=*/"periodic_initial_values_euler");
  const auto& u_0 = periodic_initial_values_euler;
  visualizer(u_0, "initial_values", "");

  using S = FvSpace<GL, R, m>;
  S space(grid_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs" << std::endl;
  std::cout << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");

  project(u_0, initial_values);
  visualizer(initial_values, "projected_initial_values", "");

  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> euler_1d(
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        if (d != 1)
          DUNE_THROW(NotImplemented, "Only for 1d!\nd = " << d);
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        const auto& v = primitive_variables[1];
        const auto& e = conservative_variables[2];
        const auto& p = primitive_variables[2];
        FieldMatrix<R, d, m> ret;
        auto& f = ret[0];
        f[0] = rho * v;
        f[1] = rho * v * v + p;
        f[2] = (e + p) * v;
        check_values(ret);
        return ret;
      },
      {},
      "euler_1d",
      [](const auto& /*mu*/) { return 3; },
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        if (d != 1)
          DUNE_THROW(NotImplemented, "Only for 1d!\nd = " << d);
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        const auto& v = primitive_variables[1];
        const auto& e = conservative_variables[2];
        FieldVector<FieldMatrix<R, m, m>, d> ret;
        auto& jacobian_f = ret[0];
        jacobian_f[0][0] = 0.;
        jacobian_f[0][1] = 1.;
        jacobian_f[0][2] = 0.;
        jacobian_f[1][0] = 0.5 * (gamma - 3.) * v * v;
        jacobian_f[1][1] = (3. - gamma) * v;
        jacobian_f[1][2] = gamma - 1.;
        jacobian_f[2][0] = v * ((gamma - 1.) * v * v - gamma * (e / rho));
        jacobian_f[2][1] = gamma * (e / rho) - (3. * (gamma - 1.) / 2.) * v * v;
        jacobian_f[2][2] = gamma * v;
        check_values(jacobian_f);
        return ret;
      });

  // See DF2016, 8.1.1, pp. 402 - 405
  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> euler_2d(
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        if (d != 2)
          DUNE_THROW(NotImplemented, "Only for 2d!\nd = " << d);
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        DomainType v;
        for (size_t ii = 0; ii < d; ++ii)
          v = primitive_variables[ii + 1];
        const auto& e = conservative_variables[m - 1];
        const auto& p = primitive_variables[m - 1];
        FieldMatrix<R, d, m> ret;
        for (size_t ss = 0; ss < d; ++ss) {
          auto& f_s = ret[ss];
          f_s[0] = rho * v[ss];
          for (size_t ii = 0; ii < d; ++ii)
            f_s[1 + ii] = rho * v[ii] * v[ss] + (ss == 1 ? 1 : 0) * p;
          f_s[m - 1] = (e + p) * v[ss];
        }
        check_values(ret);
        return ret;
      },
      {},
      "euler_2d",
      [](const auto& /*mu*/) { return 3; },
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        if (d != 2)
          DUNE_THROW(NotImplemented, "Only for 2d!\nd = " << d);
        check_values(conservative_variables);
        const auto primitive_variables = to_primitive(conservative_variables);
        check_values(primitive_variables);
        const auto& rho = conservative_variables[0];
        DomainType v;
        for (size_t ii = 0; ii < d; ++ii)
          v = primitive_variables[ii + 1];
        const auto& e = conservative_variables[m - 1];
        FieldVector<FieldMatrix<R, m, m>, d> ret;
        const auto gamma_1 = gamma - 1.;
        // f_0
        auto& jacobian_f_0 = ret[0];
        jacobian_f_0[0] = {0., 1., 0., 0.};
        jacobian_f_0[1][0] = 0.5 * gamma_1 * v.two_norm2() - v[0] * v[0];
        jacobian_f_0[1][1] = (3. - gamma) * v[0];
        jacobian_f_0[1][2] = -1. * gamma_1 * v[1];
        jacobian_f_0[1][3] = gamma_1;
        jacobian_f_0[2][0] = -1. * v[0] * v[1];
        jacobian_f_0[2][1] = v[1];
        jacobian_f_0[2][2] = v[0];
        jacobian_f_0[2][3] = 0.;
        jacobian_f_0[3][0] = v[0] * (gamma_1 * v.two_norm2() - (gamma * e) / rho);
        jacobian_f_0[3][1] = ((gamma * e) / rho) - gamma_1 * v[0] * v[0] - 0.5 * gamma_1 * v.two_norm2();
        jacobian_f_0[3][2] = -1. * gamma_1 * v[0] * v[1];
        jacobian_f_0[3][3] = gamma * v[0];
        check_values(jacobian_f_0);
        // f_1
        auto& jacobian_f_1 = ret[1];
        jacobian_f_1[0] = {0., 0., 1., 0.};
        jacobian_f_1[1][0] = -1. * v[0] * v[1];
        jacobian_f_1[1][1] = v[1];
        jacobian_f_1[1][2] = v[0];
        jacobian_f_1[1][3] = 0.;
        jacobian_f_1[2][0] = 0.5 * gamma_1 * v.two_norm2() - v[1] * v[1];
        jacobian_f_1[2][1] = -1. * gamma_1 * v[0];
        jacobian_f_1[2][2] = (3. - gamma) * v[1];
        jacobian_f_1[2][3] = gamma_1;
        jacobian_f_1[3][0] = v[1] * (gamma_1 * v.two_norm2() - ((gamma * e) / rho));
        jacobian_f_1[3][1] = -1. * gamma_1 * v[0] * v[1];
        jacobian_f_1[3][2] = ((gamma * e) / rho) - gamma_1 * v[1] * v[1] - 0.5 * gamma_1 * v.two_norm2();
        jacobian_f_1[3][3] = gamma * v[1];
        check_values(jacobian_f_1);
        return ret;
      });
  const auto& flux = euler_2d;

  auto numerical_flux = GDT::make_numerical_vijayasundaram_euler_flux(flux, gamma);
  using OpType = GDT::AdvectionFvOperator<DF>;
  OpType advec_op(grid_layer /*, flux*/, numerical_flux);

  // compute dt via Cockburn, Coquel, LeFloch, 1995
  // (in general, looking for the min/max should also include the boundary data)
  FieldVector<R, m> data_minimum(std::numeric_limits<R>::max());
  FieldVector<R, m> data_maximum(std::numeric_limits<R>::min());
  for (auto&& entity : elements(grid_layer)) {
    const auto u0_local = u_0.local_function(entity);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(entity.type(), u0_local->order())) {
      const auto value = u0_local->evaluate(quadrature_point.position());
      for (size_t ii = 0; ii < m; ++ii) {
        data_minimum[ii] = std::min(data_minimum[ii], value[ii]);
        data_maximum[ii] = std::max(data_maximum[ii], value[ii]);
      }
    }
  }
  R max_flux_derivative = std::numeric_limits<R>::min();
  if (XT::Common::FloatCmp::eq(data_minimum, data_maximum)) {
    const auto df = flux.partial_u({}, data_minimum);
    for (size_t ss = 0; ss < d; ++ss)
      max_flux_derivative = std::max(max_flux_derivative, df[ss].infinity_norm());
  } else {
    const auto max_flux_grid = XT::Grid::make_cube_grid<YaspGrid<m, EquidistantOffsetCoordinates<double, m>>>(
        data_minimum, data_maximum, XT::Common::FieldVector<unsigned int, m>(1));
    const auto max_flux_interval = *max_flux_grid.leaf_view().template begin<0>();
    for (const auto& quadrature_point : QuadratureRules<R, m>::rule(max_flux_interval.type(), flux.order())) {
      const auto df = flux.partial_u({}, max_flux_interval.geometry().global(quadrature_point.position()));
      for (size_t ss = 0; ss < d; ++ss)
        max_flux_derivative = std::max(max_flux_derivative, df[ss].infinity_norm());
    }
  }
  D perimeter_over_volume = std::numeric_limits<D>::min();
  for (auto&& entity : elements(grid_layer)) {
    D perimeter = 0;
    for (auto&& intersection : intersections(grid_layer, entity))
      perimeter += intersection.geometry().volume();
    perimeter_over_volume = std::max(perimeter_over_volume, perimeter / entity.geometry().volume());
  }
  const auto dt = 1. / (perimeter_over_volume * max_flux_derivative);

  const double T = 5.;
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
  const auto test_dt = time_stepper.find_suitable_dt(dt, 10, 1.1 * T * initial_values.vector().sup_norm(), 25);
  if (!test_dt.first)
    DUNE_THROW(InvalidStateException,
               "Could not determine optimal dt (in particular, the dt computed to match the CFL "
               "condition did not yield a stable scheme)!");
  if (test_dt.second < dt) {
    DUNE_THROW(InvalidStateException,
               "The computed dt (to match the CFL condition) does not yield a stable scheme: "
                   << dt
                   << "\n   The following dt seems to work fine: "
                   << test_dt.second);
  } else
    time_stepper.solve(T, dt, std::min(100, int(T / dt)), false, true, "solution", visualizer);
}
