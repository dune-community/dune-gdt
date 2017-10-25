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

// using G = YASP_1D_EQUIDISTANT_OFFSET;
using G = ALU_2D_SIMPLEX_CONFORMING;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t r = 1;

using DomainType = XT::Common::FieldVector<D, d>;

using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;
using F = XT::Functions::LocalizableFluxFunctionInterface<E, D, d, U, 0, R, d>;


GTEST_TEST(hyperbolic, scalar_equation)
{
  auto grid = XT::Grid::make_cube_grid<G>(DomainType(0.), DomainType(1.), XT::Common::FieldVector<unsigned int, d>(32));
  grid.global_refine(1);

  auto leaf_layer = grid.leaf_view();
  //  auto leaf_layer = grid.leaf_part();
  std::cout << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;
  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  using GL = decltype(periodic_leaf_layer);

  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, r>;
  U0 bump(
      [](const auto& xx, const auto& /*mu*/) {
        return std::exp(-0.5 * (((xx - DomainType(0.5)) * (xx - DomainType(0.5))) / std::pow(0.1, d)));
      },
      3,
      {},
      "bump");
  U0 indicator(
      [](const auto& xx, const auto& /*mu*/) {
        if (XT::Common::FloatCmp::ge(xx, DomainType(0.25)) && XT::Common::FloatCmp::le(xx, DomainType(0.5)))
          return 1.;
        else
          return 0.;
      },
      0,
      {},
      "indicator");
  const auto& u_0 = indicator;
  u_0.visualize(periodic_leaf_layer, "initial_values");

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

  using S = FvSpace<GL, R, r>;
  //  using S = DuneFemDgSpaceWrapper<GL, 1, R, r>;
  S space(periodic_leaf_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs" << std::endl;
  std::cout << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");
  project(u_0, initial_values);

  using OpType = GDT::AdvectionFvOperator<DF>;
  //  using OpType = GDT::AdvectionDgOperator<DF>;
  auto numerical_flux = GDT::make_numerical_engquist_osher_flux(flux);
  OpType advec_op(periodic_leaf_layer, numerical_flux);

  const auto dt = GDT::estimate_dt_for_scalar_advection_equation(periodic_leaf_layer, u_0, flux);
  const double T = 10.;

  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
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
  time_stepper.solve(T, dt, std::min(100, int(T / dt)), false, true);
}
