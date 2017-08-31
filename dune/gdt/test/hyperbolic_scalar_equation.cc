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


template <class GL, class S, class V>
class FvOperator
{
  using DF = DiscreteFunction<S, V>;
  using I = XT::Grid::extract_intersection_t<GL>;

public:
  FvOperator(const S& s, const F& f)
    : s_(s)
    , f_(f)
  {
  }

  void apply(const DF& source, DF& range, const XT::Common::Parameter& /*mu*/ = {}) const
  {
    range.vector() *= 0.;

    auto numerical_flux = [](const auto& flux, const auto& u, const auto& v, const auto& n) {
      const auto df = flux.partial_u({}, (u + v) / 2.);
      if ((n * df) > 0)
        return flux.evaluate({}, u) * n;
      else
        return flux.evaluate({}, v) * n;
    };

    for (auto&& entity : Dune::elements(s_.grid_layer())) {
      const auto h = entity.geometry().volume();
      const auto u_inside = source.local_discrete_function(entity);
      for (auto&& intersection : Dune::intersections(s_.grid_layer(), entity)) {
        const auto neighbor = intersection.outside();
        const auto normal = intersection.centerUnitOuterNormal();
        const auto u_outside = source.local_discrete_function(neighbor);
        const auto g =
            numerical_flux(*f_.local_function(entity), u_inside->vector().get(0), u_outside->vector().get(0), normal);
        range.local_discrete_function(entity)->vector().add(0, (g * intersection.geometry().volume()) / h);
      }
    }
  } // ... apply(...)

private:
  const S& s_;
  const F& f_;
}; // class FvOperator


template <class GL, class S, class V>
class DgOperator
{
  using DF = DiscreteFunction<S, V>;
  using I = XT::Grid::extract_intersection_t<GL>;

public:
  DgOperator(const S& s, const F& f)
    : s_(s)
    , f_(f)
  {
  }

  void apply(const DF& source, DF& range, const XT::Common::Parameter& /*mu*/ = {}) const
  {
    range.vector() *= 0.;

    auto numerical_flux = [](const auto& flux, const auto& u, const auto& v, const auto& n) {
      const auto df = flux.partial_u({}, (u + v) / 2.);
      if ((n * df) > 0)
        return flux.evaluate({}, u) * n;
      else
        return flux.evaluate({}, v) * n;
    };

    for (auto&& entity : Dune::elements(s_.grid_layer())) {
      const auto h = entity.geometry().volume();
      // prepare local functions
      const auto basis_en = s_.base_function_set(entity);
      const auto local_flux = f_.local_function(entity);
      const auto local_source_en = source.local_function(entity);
      const auto local_range = range.local_discrete_function(entity);
      // prepare local containers
      auto& local_range_vector = local_range->vector();
      if (local_range_vector.size() != basis_en.size())
        DUNE_THROW(InvalidStateException, "");
      XT::LA::CommonDenseMatrix<R> local_mass_matrix(basis_en.size(), basis_en.size(), 0.);
      XT::LA::CommonDenseVector<R> local_advection_op(basis_en.size(), 0.);
      // create volume quadrature
      const auto volume_integrand_order =
          (local_flux->order() * local_source_en->order()) + (std::max(ssize_t(0), ssize_t(basis_en.order()) - 1));
      // loop over all quadrature points
      for (const auto& quadrature_point :
           QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
        const auto xx = quadrature_point.position();
        // integration factors
        const auto integration_factor = entity.geometry().integrationElement(xx);
        const auto quadrature_weight = quadrature_point.weight();
        // evaluate everything
        const auto basis_values = basis_en.evaluate(xx);
        const auto basis_jacobians = basis_en.jacobian(xx);
        const auto source_value = local_source_en->evaluate(xx);
        const auto flux_value = local_flux->evaluate(xx, source_value);
        // compute mass matrix
        for (size_t ii = 0; ii < basis_en.size(); ++ii)
          for (size_t jj = 0; jj < basis_en.size(); ++jj)
            local_mass_matrix.add_to_entry(
                ii, jj, integration_factor * quadrature_weight * (basis_values[ii] * basis_values[jj]));
        // compute volume part of advection operator
        for (size_t jj = 0; jj < basis_en.size(); ++jj)
          local_advection_op[jj] +=
              integration_factor * quadrature_weight * -1. * (basis_jacobians[jj][0] * flux_value);
      }
      // walk the neighbors
      for (auto&& intersection : Dune::intersections(s_.grid_layer(), entity)) {
        const auto neighbor = intersection.outside();
        const auto local_source_ne = source.local_function(neighbor);
        const auto face_integrand_order =
            basis_en.order() + local_flux->order() * std::max(local_source_en->order(), local_source_ne->order());
        for (const auto& quadrature_point :
             QuadratureRules<D, d - 1>::rule(intersection.type(), boost::numeric_cast<int>(face_integrand_order))) {
          const auto x_intersection = quadrature_point.position();
          const auto integration_factor = intersection.geometry().integrationElement(x_intersection);
          const auto quadrature_weight = quadrature_point.weight();
          const auto normal = intersection.unitOuterNormal(x_intersection);
          const auto x_entity = intersection.geometryInInside().global(x_intersection);
          const auto x_neighbor = intersection.geometryInOutside().global(x_intersection);
          // evaluate everything
          const auto basis_values = basis_en.evaluate(x_entity);
          const auto g = numerical_flux(
              *local_flux, local_source_en->evaluate(x_entity), local_source_ne->evaluate(x_neighbor), normal);
          for (size_t jj = 0; jj < basis_en.size(); ++jj) {
            local_advection_op[jj] += integration_factor * quadrature_weight * (g * basis_values[jj]);
          }
        }
      }
      XT::LA::CommonDenseVector<R> local_DoF_vector(basis_en.size(), 0.);
      XT::LA::make_solver(local_mass_matrix).apply(local_advection_op, local_DoF_vector);
      for (size_t ii = 0; ii < basis_en.size(); ++ii)
        local_range_vector.set(ii, local_DoF_vector[ii]);
    }
  } // ... apply(...)

private:
  const S& s_;
  const F& f_;
}; // class DgOperator


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
  auto flux = burgers;

  using S = FvSpace<GL, R, r>;
  //  using S = DuneFemDgSpaceWrapper<GL, 1, R, r>;
  S space(periodic_leaf_layer);
  std::cout << "space has " << space.mapper().size() << " DoFs" << std::endl;
  std::cout << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");

  project(u_0, initial_values);

  //  using OpType = FvOperator<GL, S, V>;
  //  OpType advec_op(space, flux);
  //  using OpType = DgOperator<GL, S, V>;
  //  OpType advec_op(space, flux);
  using OpType = GDT::AdvectionFvOperator<DF>;
  //  using OpType = GDT::AdvectionDgOperator<DF>;
  auto upwinding = [&](const auto& /*f*/, const auto& u, const auto& v, const auto& n, const auto& /*mu*/) {
    const auto df = flux.partial_u({}, (u + v) / 2.);
    if ((n * df) > 0)
      return flux.evaluate({}, u) * n;
    else
      return flux.evaluate({}, v) * n;
  };
  auto lax_friedrichs = [&](const auto& /*f*/, const auto& u, const auto& v, const auto& n, const auto& /*mu*/) {
    const auto lambda = 1. / std::max(flux.partial_u({}, u).infinity_norm(), flux.partial_u({}, v).infinity_norm());
    return 0.5 * ((flux.evaluate({}, u) + flux.evaluate({}, v)) * n) + 0.5 * ((u - v) / lambda);
  };
  auto engquist_osher = [&](const auto& /*f*/, const auto& u, const auto& v, const auto& n, const auto& /*mu*/) {
    auto integrate_f = [&](const auto& s, const std::function<double(const double&, const double&)>& min_max) {
      if (XT::Common::FloatCmp::eq(s[0], 0.))
        return 0.;
      D ret = 0.;
      const OneDGrid state_grid(1, 0., s[0]);
      const auto state_interval = *state_grid.leafGridView().template begin<0>();
      for (const auto& quadrature_point : QuadratureRules<double, 1>::rule(state_interval.type(), flux.order())) {
        const auto local_uu = quadrature_point.position();
        const auto uu = state_interval.geometry().global(local_uu);
        const auto df = flux.partial_u({}, uu);
        ret += state_interval.geometry().integrationElement(local_uu) * quadrature_point.weight() * min_max(n * df, 0.);
      }
      return ret;
    };
    return (flux.evaluate({}, 0.) * n) + integrate_f(u, [](const double& a, const double& b) { return std::max(a, b); })
           + integrate_f(v, [](const double& a, const double& b) { return std::min(a, b); });
  };
  auto numerical_flux = engquist_osher;
  OpType advec_op(periodic_leaf_layer, flux, numerical_flux);

  // compute dt via Cockburn, Coquel, LeFloch, 1995
  // (in general, looking for the min/max should also include the boundary data)
  R data_minimum = std::numeric_limits<R>::max();
  R data_maximum = std::numeric_limits<R>::min();
  for (auto&& entity : elements(periodic_leaf_layer)) {
    const auto u0_local = u_0.local_function(entity);
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(entity.type(), u0_local->order())) {
      const auto value = u0_local->evaluate(quadrature_point.position());
      data_minimum = std::min(data_minimum, value[0]);
      data_maximum = std::max(data_maximum, value[0]);
    }
  }
  R max_flux_derivative = std::numeric_limits<R>::min();
  OneDGrid max_flux_grid(1, data_minimum, data_maximum);
  const auto max_flux_interval = *max_flux_grid.leafGridView().template begin<0>();
  for (const auto& quadrature_point : QuadratureRules<R, r>::rule(max_flux_interval.type(), flux.order())) {
    const auto df = flux.partial_u({}, max_flux_interval.geometry().global(quadrature_point.position()));
    max_flux_derivative = std::max(max_flux_derivative, df.infinity_norm());
  }
  D perimeter_over_volume = std::numeric_limits<D>::min();
  for (auto&& entity : elements(periodic_leaf_layer)) {
    D perimeter = 0;
    for (auto&& intersection : intersections(periodic_leaf_layer, entity))
      perimeter += intersection.geometry().volume();
    perimeter_over_volume = std::max(perimeter_over_volume, perimeter /= entity.geometry().volume());
  }
  const auto dt = 1. / (perimeter_over_volume * max_flux_derivative);

  const double T = 5.;
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_rungekutta_third_order_ssp> time_stepper(
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
