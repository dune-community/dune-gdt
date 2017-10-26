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

#include <cmath>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
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
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>


using namespace Dune;
using namespace Dune::GDT;

using G = YASP_1D_EQUIDISTANT_OFFSET;
// using G = ALU_2D_SIMPLEX_CONFORMING;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t r = 1;

using DomainType = XT::Common::FieldVector<D, d>;

using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;


GTEST_TEST(hyperbolic, scalar_equation)
{
  auto grid =
      XT::Grid::make_cube_grid<G>(DomainType(-1.), DomainType(1.), XT::Common::FieldVector<unsigned int, d>(128));
  grid.global_refine(1);

  auto leaf_layer = grid.leaf_view();
  //  auto leaf_layer = grid.leaf_part(); // for DG with dune-fem
  //  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  auto& grid_layer = leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;
  using I = XT::Grid::extract_intersection_t<GL>;
  std::cout << "grid has " << grid_layer.indexSet().size(0) << " elements" << std::endl;

  using S = FvSpace<GL, R, r>;
  //  using S = DuneFemDgSpaceWrapper<GL, 1, R, r>;
  S space(grid_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs\n" << std::endl;

  // initial values
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, r>;
  U0 inflow_from_the_left(
      [](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::le(xx, DomainType(-0.9))
            || (XT::Common::FloatCmp::ge(xx, DomainType(-0.25)) && XT::Common::FloatCmp::le(xx, DomainType(0.))))
          return 1.;
        else
          return 0.;
      },
      0,
      {},
      "inflow_from_the_left");
  U0 smooth_inflow_from_the_left(
      [](const auto& xx, const auto& /*mu*/) {
        return std::exp(-0.5 * (((xx - DomainType(-1.)) * (xx - DomainType(-1.))) / std::pow(0.1, d)));
      },
      0,
      {},
      "smooth_inflow_from_the_left");
  U0 indicator(
      [](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0.)))
          return 1.;
        else
          return 0.;
      },
      0,
      {},
      "indicator");
  const auto& u_0 = smooth_inflow_from_the_left;
  u_0.visualize(grid_layer, "initial_values");

  // boundary values (needs to match grid layer)
  XT::Grid::NormalBasedBoundaryInfo<I> boundary_info;
  boundary_info.register_new_normal({-1.}, new XT::Grid::InflowOutflowBoundary());
  boundary_info.register_new_normal({1.}, new XT::Grid::InflowOutflowBoundary());

  U0 zero([](const auto& /*xx*/, const auto& /*mu*/) { return 0.; }, 0, {}, "zero");
  const auto& bv = smooth_inflow_from_the_left;
  bv.visualize(grid_layer, "boundary_values");

  const auto inflow_outflow_boundary_treatment =
      [&](const auto& intersection, const auto& x_intersection, const auto& f, const auto& u, const auto& /*mu*/ = {}) {
        const auto entity = intersection.inside();
        const auto normal = intersection.unitOuterNormal(x_intersection);
        const auto df = f.partial_u({}, u);
        if ((normal * df) > 0) { // inflow
          const auto local_bv = bv.local_function(entity);
          return local_bv->evaluate(intersection.geometryInInside().global(x_intersection));
        } else { // outflow
          return u;
        }
      };

  // flux
  XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, r> linear_transport(
      [](const auto& /*x*/, const auto& u, const auto& /*mu*/) { return XT::Common::FieldVector<R, d>(u[0]); },
      {},
      "linear_transport",
      [](const XT::Common::Parameter& /*mu*/) { return 1; },
      [](const auto& /*x*/, const auto& /*u*/, const auto& /*mu*/) { return 1.; });

  XT::Functions::GlobalLambdaFluxFunction<U, 0, R, d, r> burgers(
      [](const auto& /*x*/, const auto& u, const auto& /*mu*/) {
        return XT::Common::FieldVector<R, d>(std::pow(u[0], 2.) / 2.);
      },
      {},
      "burgers_flux",
      [](const XT::Common::Parameter& /*mu*/) { return 2; },
      [](const auto& /*x*/, const auto& u, const auto& /*mu*/) { return u[0]; });
  const auto& flux = burgers;

  // inital projection
  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");
  project(u_0, initial_values);

  // spatial discretization
  using OpType = GDT::AdvectionFvOperator<DF>;
  //  using OpType = GDT::AdvectionDgOperator<DF>;
  auto numerical_flux = GDT::make_numerical_engquist_osher_flux(flux);
  OpType advec_op(grid_layer, numerical_flux);
  advec_op.append(
      inflow_outflow_boundary_treatment,
      new XT::Grid::ApplyOn::CustomBoundaryIntersections<GL>(boundary_info, new XT::Grid::InflowOutflowBoundary()));

  // temporal discretization
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
  const auto dt = GDT::estimate_dt_for_scalar_advection_equation(grid_layer, u_0, flux);
  const auto test_dt = time_stepper.find_suitable_dt(dt, 10, 1.1 * initial_values.vector().sup_norm(), 25);
  if (!test_dt.first)
    DUNE_THROW(InvalidStateException,
               "Could not determine optimal dt (in particular, the dt computed to match the CFL "
               "conditiond did not yield a stable scheme)!");
  if (test_dt.second < dt)
    DUNE_THROW(InvalidStateException,
               "The computed dt (to match the CFL condition) does not yield a stable scheme: "
                   << dt
                   << "\n   The following dt seems to work fine: "
                   << test_dt.second);

  const double T = 5.;
  time_stepper.solve(T, dt, std::min(100, int(T / dt)), false, true);
}
