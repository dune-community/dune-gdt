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

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/lambda.hh>
#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/euler.hh>

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
using RangeType = XT::Common::FieldVector<D, m>;


GTEST_TEST(empty, main)
{
  auto grid =
      XT::Grid::make_cube_grid<G>(DomainType(-1.), DomainType(1.), XT::Common::FieldVector<unsigned int, d>(128));
  grid.global_refine(1);

  auto leaf_layer = grid.leaf_view();
  std::cout << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;
  //  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  auto& grid_layer = /*periodic_*/ leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;

  const double gamma = 1.4; // air or water at roughly 20 deg Cels.
  EulerTools<d> euler_tools(gamma);

  using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
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
        return euler_tools.to_conservative(primitive_variables);
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
        return euler_tools.to_conservative(primitive_variables);
      },
      /*order=*/0,
      /*parameter_type=*/{},
      /*name=*/"periodic_initial_values_euler");
  const auto& u_0 = periodic_initial_values_euler;
  euler_tools.visualize(u_0, grid_layer, "initial_values");

  using I = XT::Grid::extract_intersection_t<GL>;
  XT::Grid::NormalBasedBoundaryInfo<I> boundary_info;
  boundary_info.register_new_normal({-1., 0.}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_normal({1., 0.}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_normal({0., -1.}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_normal({0., 1.}, new XT::Grid::ImpermeableBoundary());
  XT::Grid::ApplyOn::CustomBoundaryIntersections<GL> impermeable_wall_filter(boundary_info,
                                                                             new XT::Grid::ImpermeableBoundary());

  // see [DF2015, p. 414, (8.58)]
  const auto euler_impermeable_wall_treatment = [&](const auto& source, const auto& intersection, auto& local_range) {
    const auto& entity = local_range.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto normal = intersection.centerUnitOuterNormal();
    RangeType u_cons;
    if (u_cons.size() != u_inside->vector().size())
      DUNE_THROW(InvalidStateException, "");
    for (size_t ii = 0; ii < u_cons.size(); ++ii)
      u_cons[ii] = u_inside->vector().get(ii);
    const auto pressure = euler_tools.to_primitive(u_cons)[m - 1];
    RangeType g(0.);
    const auto tmp = normal * pressure;
    for (size_t ii = 0; ii < d; ++ii)
      g[ii + 1] = tmp[ii];
    const auto h = local_range.entity().geometry().volume();
    for (size_t ii = 0; ii < m; ++ii)
      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  };

  const auto& impermeable_wall_treatment = euler_impermeable_wall_treatment;

  using S = FvSpace<GL, R, m>;
  S space(grid_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs" << std::endl;
  std::cout << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");

  project(u_0, initial_values);
  euler_tools.visualize(initial_values, grid_layer, "projected_initial_values");

  const XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, m> euler_flux(
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        return euler_tools.flux(conservative_variables);
      },
      {},
      "euler_flux",
      [](const auto& /*mu*/) { return 3; },
      [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
        return euler_tools.flux_jacobian(conservative_variables);
      });
  const auto& flux = euler_flux;

  auto numerical_flux = GDT::make_numerical_vijayasundaram_euler_flux(flux, gamma);
  using OpType = GDT::AdvectionFvOperator<DF>;
  OpType advec_op(grid_layer, numerical_flux, /*periodicity_restriction=*/impermeable_wall_filter.copy());
  // impermeable wall
  advec_op.append(impermeable_wall_treatment, impermeable_wall_filter.copy());

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
    time_stepper.solve(T,
                       dt,
                       std::min(100, int(T / dt)),
                       false,
                       true,
                       "solution",
                       /*visualizer=*/[&](const auto& u, const auto& filename_prefix, const auto& step) {
                         euler_tools.visualize(u, grid_layer, filename_prefix, XT::Common::to_string(step));
                       });
}
