// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <algorithm>
#include <cmath>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/vector.hh>
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
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/euler.hh>

using namespace Dune;
using namespace Dune::GDT;

using G = YASP_1D_EQUIDISTANT_OFFSET;
using E = typename G::template Codim<0>::Entity;
using D = double;
static const constexpr size_t d = G::dimension;
using R = double;
static const constexpr size_t m = d + 2;

using DomainType = XT::Common::FieldVector<D, d>;
using RangeType = XT::Common::FieldVector<D, m>;


GTEST_TEST(inviscid_compressible_flow__euler_1d_explicit_fv, periodic_boundary_values)
{
  auto logger = XT::Common::TimedLogger().get("main");

  auto grid = XT::Grid::make_cube_grid<G>(-1., 1., 128u);
  auto leaf_layer = grid.leaf_view();
  logger.info() << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;

  auto periodic_leaf_layer = XT::Grid::make_periodic_grid_layer(leaf_layer);
  auto& grid_layer = periodic_leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;

  const double gamma = 1.4; // air or water at roughly 20 deg Cels.
  EulerTools<d> euler_tools(gamma);

  using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
  const U0 periodic_initial_values_euler( // compare [Kröner, 1997, p. 394]
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

  using S = FvSpace<GL, R, m>;
  S space(grid_layer);
  logger.info() << "space has " << space.mapper().size() << " DoFs\n" << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");

  project(u_0, initial_values);

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

  const auto numerical_flux = GDT::make_numerical_vijayasundaram_flux(flux);

  using OpType = GDT::AdvectionFvOperator<DF>;
  OpType advec_op(grid_layer, numerical_flux);

  const auto dt = estimate_dt_for_hyperbolic_system(grid_layer, initial_values, flux);
  EXPECT_DOUBLE_EQ(0.002708880865541605, dt);

  // do a hundred time steps for this test
  const double T = 100 * dt;
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
  const auto test_dt =
      time_stepper.find_suitable_dt(dt, 10, 1.1 * std::max(T, 1.) * initial_values.vector().sup_norm(), 25);
  ASSERT_TRUE(test_dt.first) << "Could not determine optimal dt (in particular, the dt computed to match the CFL "
                                "condition did not yield a stable scheme)!";
  ASSERT_LE(test_dt.second, dt) << "The computed dt (to match the CFL condition) does not yield a stable scheme: " << dt
                                << "\nThe following dt seems to work fine: " << test_dt.second;
  time_stepper.solve(T,
                     dt,
                     -1,
                     /*save_solution=*/true,
                     /*visualize=*/false,
                     "solution",
                     /*visualizer=*/[&](const auto& u, const auto& filename_prefix, const auto& step) {
                       euler_tools.visualize(u, grid_layer, filename_prefix, XT::Common::to_string(step));
                     });

  // check results
  const auto initial_l_infty_norm = initial_values.vector().sup_norm();
  const auto initial_l_2_norm = make_l2_operator(grid_layer)->induced_norm(initial_values);
  double infty = 0.;
  double l2 = 0.;
  size_t saved_timesteps = 0;
  for (const auto& t_and_u : time_stepper.solution()) {
    ++saved_timesteps;
    const auto& t = t_and_u.first;
    const auto& u_conservative = t_and_u.second;
    const auto relative_l_infty_error =
        std::abs(initial_l_infty_norm - u_conservative.vector().sup_norm()) / initial_l_infty_norm;
    EXPECT_LT(relative_l_infty_error, 0.065554) << "initial_l_infty_norm = " << initial_l_infty_norm
                                                << "\nl_infty_norm = " << u_conservative.vector().sup_norm();
    const auto relative_l_2_error =
        std::abs(initial_l_2_norm - make_l2_operator(grid_layer)->induced_norm(u_conservative)) / initial_l_2_norm;
    EXPECT_LT(relative_l_2_error, 0.10577)
        << "initial_l_2_norm = " << initial_l_2_norm
        << "\nl_2_norm = " << make_l2_operator(grid_layer)->induced_norm(u_conservative);
    infty = std::max(infty, relative_l_infty_error);
    l2 = std::max(l2, relative_l_2_error);
  }
  EXPECT_EQ(101, saved_timesteps);

  V expected_end_state(space.mapper().size());
  size_t counter = 0;
  for (const auto& element :
       {1.0000013466060014,      -1.007708936471061e-06,  1.0000018852500976,      1.0000044547440901,
        -3.3336389923735734e-06, 1.0000062366599656,      1.000014162794056,       -1.0598601470495452e-05,
        1.0000198280927333,      1.000043278400877,       -3.2387847239154085e-05, 1.0000605914221288,
        1.0001271363870938,      -9.5150841989372032e-05, 1.0001780050273525,      1.0003591157268059,
        -0.00026882151151833586, 1.0005028725006329,      1.0009754377211006,      -0.00073056045889799773,
        1.001366414575251,       1.0025466912916001,      -0.0019098633059258037,  1.0035707470786868,
        1.0063791334121701,      -0.0047990073260674075,  1.0089640383862377,      1.0152525334758986,
        -0.011556252795424875,   1.0215410194458294,      1.0344167689618096,      -0.026465478852591921,
        1.0491252572310434,      1.0717895217997457,      -0.056729189674759324,   1.1045402937983513,
        1.1345944672335124,      -0.11088822552498599,    1.2023203685383403,      1.2214962381923051,
        -0.19199643614110681,    1.3466923165714426,      1.318393207671563,       -0.28991309450150105,
        1.5191861689934714,      1.4064898984983827,      -0.38498980936562155,    1.6857360789318017,
        1.4740706738430702,      -0.46144563975373926,    1.8193503686036647,      1.519606001779207,
        -0.51450081921013058,    1.9119474298981303,      1.5478521082711401,      -0.54757063753091517,
        1.9693139550900038,      1.5658078733144107,      -0.56708662728700576,    2.0019170714815573,
        1.5811005497742998,      -0.57952996128806766,    2.0194077923147566,      1.6020567438188698,
        -0.59064351046306418,    2.0288955951068219,      1.637910468064709,       -0.60550700348793163,
        2.0350692009495228,      1.6976821954586885,      -0.6283864059322446,     2.0408008920090324,
        1.7869319023449324,      -0.66181003497671198,    2.0476112699953899,      1.9032430479136619,
        -0.70510282990341322,    2.0558660436391993,      2.0334797411187977,      -0.75349613653548186,
        2.0648975771168909,      2.1565336929162089,      -0.79921077698650322,    2.0733560945837159,
        2.2522404991379013,      -0.83479169848126478,    2.0798968043590746,      2.3114265912752696,
        -0.85684760016965011,    2.083899954809131,       2.3392396063595036,      -0.86729577601409358,
        2.0857112114673102,      2.3495693584200379,      -0.8708360702921687,     2.0868782867966438,
        2.3594878090231135,      -0.8700490995609067,     2.0935220152989653,      2.3858965507748144,
        -0.86320896878437114,    2.1172843206457626,      2.4395232929593429,      -0.84710304992662333,
        2.168019040837367,       2.5213650728113812,      -0.82007172518620353,    2.247679359487039,
        2.6259471479126888,      -0.78220350934058758,    2.3522148436628951,      2.7461027591121283,
        -0.73442560285831604,    2.4758328024512104,      2.875508269532367,       -0.67801868499946605,
        2.613128798835147,       3.0090513581218814,      -0.61445262326237748,    2.7593667089725793,
        3.1424675478253601,      -0.54534702557942139,    2.9101211473906226,      3.2719132514803162,
        -0.47243300608725741,    3.0608637997446988,      3.3936513815634881,      -0.39746595752458397,
        3.2066604292704208,      3.5038736821181184,      -0.32207774266836353,    3.3420261199208121,
        3.5986697395755836,      -0.24758527153533763,    3.4609838548574712,      3.6741637361057897,
        -0.17480020554354553,    3.5573815207499746,      3.726826687177506,       -0.10390461773147322,
        3.6254886508954995,      3.7539167274869008,      -0.034455164260433178,   3.6607941606870735,
        3.7539167274868928,      0.034455164260433054,    3.6607941606870731,      3.7268266871774998,
        0.10390461773147293,     3.6254886508954995,      3.6741637361058204,      0.17480020554354664,
        3.5573815207499742,      3.5986697395755729,      0.24758527153533685,     3.4609838548574707,
        3.5038736821181176,      0.32207774266836303,     3.3420261199208112,      3.393651381563485,
        0.39746595752458352,     3.2066604292704195,      3.2719132514803055,      0.47243300608725658,
        3.0608637997446966,      3.1424675478253561,      0.54534702557942083,     2.9101211473906226,
        3.0090513581218783,      0.6144526232623776,      2.7593667089725789,      2.8755082695323657,
        0.67801868499946627,     2.6131287988351448,      2.7461027591121265,      0.73442560285831648,
        2.47583280245121,        2.6259471479126861,      0.78220350934058758,     2.3522148436628947,
        2.5213650728113803,      0.82007172518620264,     2.247679359487039,       2.4395232929593424,
        0.84710304992662344,     2.1680190408373674,      2.3858965507748127,      0.86320896878437159,
        2.1172843206457626,      2.3594878090231139,      0.87004909956090648,     2.0935220152989662,
        2.3495693584200379,      0.87083607029216892,     2.0868782867966438,      2.3392396063595027,
        0.86729577601409458,     2.0857112114673111,      2.3114265912752687,      0.85684760016965045,
        2.0838999548091315,      2.2522404991378995,      0.83479169848126455,     2.0798968043590751,
        2.156533692916208,       0.79921077698650256,     2.0733560945837142,      2.0334797411187986,
        0.75349613653548275,     2.0648975771168909,      1.9032430479136633,      0.70510282990341322,
        2.0558660436391998,      1.7869319023449348,      0.66181003497671265,     2.0476112699953917,
        1.6976821954586898,      0.62838640593224493,     2.0408008920090333,      1.6379104680647092,
        0.60550700348793218,     2.0350692009495237,      1.6020567438188695,      0.59064351046306418,
        2.0288955951068219,      1.5811005497742998,      0.57952996128806766,     2.0194077923147566,
        1.5658078733144112,      0.56708662728700587,     2.0019170714815577,      1.5478521082711409,
        0.54757063753091517,     1.9693139550900034,      1.519606001779207,       0.51450081921013024,
        1.9119474298981289,      1.4740706738430702,      0.46144563975374009,     1.8193503686036647,
        1.4064898984983805,      0.38498980936562099,     1.6857360789318021,      1.3183932076715654,
        0.2899130945015016,      1.5191861689934716,      1.221496238192302,       0.19199643614110648,
        1.3466923165714419,      1.1345944672335149,      0.11088822552498603,     1.2023203685383399,
        1.0717895217997426,      0.056729189674758818,    1.1045402937983511,      1.0344167689618122,
        0.026465478852591956,    1.0491252572310437,      1.0152525334759077,      0.011556252795424923,
        1.0215410194458294,      1.0063791334121563,      0.0047990073260674084,   1.0089640383862379,
        1.0025466912916026,      0.0019098633059258323,   1.003570747078687,       1.0009754377211053,
        0.00073056045889800716,  1.0013664145752508,      1.0003591157268044,      0.00026882151151820115,
        1.0005028725006326,      1.0001271363870972,      9.515084198936241e-05,   1.0001780050273525,
        1.0000432784008766,      3.2387847239163707e-05,  1.0000605914221286,      1.0000141627940529,
        1.059860147056282e-05,   1.0000198280927326,      1.0000044547440927,      3.3336389926430426e-06,
        1.0000062366599654,      1.0000013466059992,      1.0077089366539149e-06,  1.0000018852500976,
        1.0000003911827935,      2.9273450724623062e-07,  1.0000005476560563,      1.0000001092050823,
        8.172160862392107e-08,   1.0000001528871254,      1.0000000292989173,      2.192530327153722e-08,
        1.0000000410184859,      1.0000000075551212,      5.6537353301225452e-09,  1.0000000105771703,
        1.0000000018726529,      1.4013644091063348e-09,  1.0000000026217131,      1.0000000004462188,
        3.3392012206851674e-10,  1.0000000006247076,      1.0000000001022271,      7.6500033828831142e-11,
        1.0000000001431189,      1.0000000000225202,      1.6852238215810858e-11,  1.0000000000315279,
        1.0000000000047693,      3.5702951711201839e-12,  1.0000000000066795,      1.0000000000009728,
        7.273727071948405e-13,   1.0000000000013611,      1.0000000000001921,      1.4246226759202467e-13,
        1.0000000000002669,      1.0000000000000335,      2.6869867669859681e-14,  1.0000000000000504,
        1.000000000000006,       4.9370494679935592e-15,  1.0000000000000093,      1.000000000000002,
        7.121669797885445e-16,   1.0000000000000018,      0.99999999999999845,     1.6360592778926019e-16,
        1.0000000000000004,      1.0000000000000000,      2.1172531831551323e-16,  1.0000000000000002,
        1.0000000000000002,      6.7367146736754204e-17,  1.0000000000000000,      1.0000000000000007,
        -1.8285368399976143e-16, 1.0000000000000004,      1.0000000000000027,      -8.5652515136730364e-16,
        1.0000000000000016,      1.0000000000000064,      -4.9081778336778064e-15, 1.0000000000000091,
        1.0000000000000351,      -2.6840996035543926e-14, 1.0000000000000502,      1.0000000000001898,
        -1.425103869825509e-13,  1.0000000000002667,      1.0000000000009721,      -7.2740157882915619e-13,
        1.0000000000013611,      1.0000000000047702,      -3.5701411890705e-12,    1.0000000000066795,
        1.0000000000225202,      -1.6852382573982434e-11, 1.0000000000315281,      1.0000000001022293,
        -7.650021668251514e-11,  1.0000000001431189,      1.0000000004462193,      -3.3392008357300429e-10,
        1.0000000006247074,      1.0000000018726496,      -1.4013647940614589e-09, 1.0000000026217135,
        1.0000000075551243,      -5.6537354841045952e-09, 1.0000000105771707,      1.0000000292989182,
        -2.1925303204170072e-08, 1.0000000410184864,      1.0000001092050819,      -8.172160844106737e-08,
        1.0000001528871256,      1.0000003911827926,      -2.9273450717886344e-07, 1.0000005476560563}) {
    expected_end_state[counter] = element;
    ++counter;
  }

  EXPECT_LT((expected_end_state - time_stepper.solution().rbegin()->second.vector()).sup_norm(), 1e-15);
} // ... (inviscid_compressible_flow__euler_1d_explicit_fv, periodic_boundary_values)


GTEST_TEST(inviscid_compressible_flow__euler_1d_explicit_fv, impermeable_walls_boundary_values__direct_euler_treatment)
{
  auto logger = XT::Common::TimedLogger().get("main");

  const double gamma = 1.4; // air or water at roughly 20 deg Cels.
  EulerTools<d> euler_tools(gamma);

  auto grid = XT::Grid::make_cube_grid<G>(-1., 1., 128u);
  auto leaf_layer = grid.leaf_view();
  logger.info() << "grid has " << leaf_layer.indexSet().size(0) << " elements" << std::endl;

  auto& grid_layer = leaf_layer;
  using GL = std::decay_t<decltype(grid_layer)>;

  using I = XT::Grid::extract_intersection_t<GL>;
  XT::Grid::NormalBasedBoundaryInfo<I> boundary_info;
  boundary_info.register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
  XT::Grid::ApplyOn::CustomBoundaryIntersections<GL> impermeable_wall_filter(boundary_info,
                                                                             new XT::Grid::ImpermeableBoundary());

  const auto impermeable_wall_treatment = [&](const auto& u, const auto& n, const auto& /*mu*/ = {}) {
    return euler_tools.flux_at_impermeable_walls(u, n);
  };

  using U = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U0 = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
  const U0 periodic_initial_values_euler( // compare [Kröner, 1997, p. 394]
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

  using S = FvSpace<GL, R, m>;
  S space(grid_layer);
  logger.info() << "space has " << space.mapper().size() << " DoFs\n" << std::endl;

  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  DF initial_values(space, "solution");
  project(u_0, initial_values);

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

  const auto numerical_flux = GDT::make_numerical_vijayasundaram_flux(flux);

  using OpType = GDT::AdvectionFvOperator<DF>;
  OpType advec_op(grid_layer, numerical_flux);
  advec_op.append(impermeable_wall_treatment, impermeable_wall_filter.copy(), {});

  const auto dt = estimate_dt_for_hyperbolic_system(grid_layer, initial_values, flux);
  EXPECT_DOUBLE_EQ(0.002708880865541605, dt);

  // do 250 time steps for this test, otherwise we do not see the boundary interaction
  const double T = 250 * dt;
  ExplicitRungeKuttaTimeStepper<OpType, DF, TimeStepperMethods::explicit_euler> time_stepper(
      advec_op, initial_values, -1.);
  const auto test_dt =
      time_stepper.find_suitable_dt(dt, 10, 1.1 * std::max(T, 1.) * initial_values.vector().sup_norm(), 25);
  ASSERT_TRUE(test_dt.first) << "Could not determine optimal dt (in particular, the dt computed to match the CFL "
                                "condition did not yield a stable scheme)!";
  ASSERT_LE(test_dt.second, dt) << "The computed dt (to match the CFL condition) does not yield a stable scheme: " << dt
                                << "\nThe following dt seems to work fine: " << test_dt.second;
  time_stepper.solve(T,
                     dt,
                     100,
                     /*save_solution=*/true,
                     /*visualize=*/false,
                     "solution",
                     /*visualizer=*/[&](const auto& u, const auto& filename_prefix, const auto& step) {
                       euler_tools.visualize(u, grid_layer, filename_prefix, XT::Common::to_string(step));
                     });

  // check results
  EXPECT_EQ(101, time_stepper.solution().size());
  V expected_end_state(space.mapper().size());
  size_t counter = 0;
  for (const auto& element :
       {2.4287195271302164,     -0.00032422691990605974, 3.4982911844252076,     2.4079264174279666,
        -0.0013236992919630915, 3.4964005887424094,      2.4020195519885088,     -0.005429704973102965,
        3.4894456975913135,     2.3910057771571176,      -0.015430514252899416,  3.4680874772877628,
        2.3671325956594558,     -0.03681284480133995,    3.4204281435472983,     2.3223027395944076,
        -0.076882808226762786,  3.3298244516393507,      2.247215332769346,      -0.14334466015756481,
        3.1777311047779047,     2.1358370166279212,      -0.23914173028958638,   2.9547368361970494,
        1.9959075430474098,     -0.35331902236289886,    2.6818329724329226,     1.8578208815503827,
        -0.45945860225814217,   2.4192491095573119,      1.7600203908629302,     -0.53559431913273403,
        2.227516684206253,      1.7183065469021335,      -0.58211857588721072,   2.1188750878231057,
        1.7228085692439847,     -0.61191739912392462,    2.0674581756801591,     1.7578944079137815,
        -0.63567274027409904,   2.0452841583958734,      1.8124288123169117,     -0.65861032280541387,
        2.0356769889792279,     1.8790140822246693,      -0.68220885360217198,   2.0305954405918891,
        1.9514766613139378,     -0.70587118878439536,    2.0262780422278581,     2.0237357660082504,
        -0.72795348690668049,   2.0207890549070373,      2.0899257202405517,     -0.74646992668323597,
        2.012971525667032,      2.1451267857450351,      -0.75964596201599388,   2.0020760214414119,
        2.1861330890686328,     -0.76631791424714935,    1.9876567119402357,     2.2118527116552751,
        -0.76611978464238084,   1.969541510710995,       2.2231669489282737,     -0.75943459575539407,
        1.9477951402236462,     2.2223287931608731,      -0.74716282056784367,   1.9226555596917831,
        2.2121602600364119,     -0.73041791429339975,    1.8944569067319204,     2.1953428197368718,
        -0.71026377431039522,   1.863561983608141,       2.1739923957917644,     -0.68756353662011638,
        1.8303201879658415,     2.1495493657717115,      -0.66294550697937571,   1.7950526439427796,
        2.1228882357616161,     -0.63684616339676747,    1.7580553303981474,     2.0945108979963489,
        -0.6095786784778543,    1.7196081946401147,      2.0647204031744333,     -0.58139088518987703,
        1.6799820332730566,     2.0337329743296411,      -0.55249995638745231,   1.6394406482425259,
        2.0017339449653608,     -0.52310689990551429,    1.5982405641151154,     1.9689055256099013,
        -0.49339739850011083,   1.5566354288068398,      1.935457990335868,      -0.46353054341056044,
        1.5148975049137092,     1.9016860415298071,      -0.43361267306640294,   1.4733669968451444,
        1.8680450233870272,     -0.40365871444303958,    1.4325205940132388,     1.8352042708543186,
        -0.37355805581217666,   1.393019126026023,       1.8040189975629839,     -0.34307205389445594,
        1.3556815756466531,     1.7753977975129291,      -0.31187942990633011,   1.3213666532565709,
        1.7501129792931043,     -0.27965689197196075,    1.2908063394169014,     1.7286445812562943,
        -0.24615980842683699,   1.2644722554710461,      1.7111260133424659,     -0.2112711793273882,
        1.2425321830672,        1.6973953566344131,      -0.17500977420015407,   1.2248969402433736,
        1.6871073841761828,     -0.13750886816931551,    1.2113169313722738,     1.6798538375311267,
        -0.098983531677978201,  1.2014835440670222,      1.6752596875355579,     -0.059699880581338682,
        1.1951090651490672,     1.6730465473358933,      -0.019952103235044864,  1.1919787001022739,
        1.6730465473358895,     0.019952103235044608,    1.1919787001022739,     1.675259687535553,
        0.059699880581338897,   1.1951090651490677,      1.6798538375311298,     0.098983531677978465,
        1.2014835440670224,     1.6871073841761857,      0.13750886816931607,    1.2113169313722747,
        1.6973953566344138,     0.17500977420015432,     1.2248969402433745,     1.7111260133424655,
        0.21127117932738818,    1.2425321830672009,      1.7286445812562949,     0.24615980842683713,
        1.2644722554710472,     1.7501129792931076,      0.27965689197196147,    1.2908063394169031,
        1.7753977975129329,     0.31187942990633077,     1.3213666532565722,     1.8040189975629859,
        0.34307205389445694,    1.3556815756466538,      1.8352042708543197,     0.37355805581217755,
        1.3930191260260223,     1.8680450233870236,      0.40365871444304202,    1.4325205940132335,
        1.9016860415297849,     0.4336126730664121,      1.473366996845124,      1.9354579903357749,
        0.46353054341059924,    1.5148975049136204,      1.9689055256095054,     0.49339739850026765,
        1.5566354288064705,     2.0017339449637337,      0.52310689990614434,    1.5982405641135942,
        2.0337329743230459,     0.55249995638995475,     1.6394406482363471,     2.0647204031480832,
        0.58139088519967408,    1.6799820332483066,      2.0945108978926545,     0.60957867851565462,
        1.7196081945424551,     2.1228882353600724,      0.6368461635403958,     1.7580553300189641,
        2.1495493642434669,     0.66294550751635406,     1.7950526424956172,     2.1739923900830673,
        0.6875635385936838,     1.8303201825430953,      2.1953427988424483,     0.7102637814337931,
        1.8635619636824288,     2.2121601852687687,      0.73041793951433931,    1.8944568350381046,
        2.2223285323385449,     0.74716290803581886,     1.9226553075270334,     2.2231660650343956,
        0.75943489233234107,    1.9477942749954515,      2.211849813093322,      0.76612076545101593,
        1.9695386214926607,     2.1861239257627192,      0.76632106887383389,    1.9876473467973301,
        2.1450989402555543,     0.75965580042300573,     2.0020466257069098,     2.0898444789940545,
        0.74649960246276781,    2.0128823341164122,      2.023507933485003,      0.7280399393010234,
        2.0205276154921386,     1.9508602733696871,      0.70611454078487368,    2.0255370297272588,
        1.8773963072668505,     0.68287237139756984,     2.0285596283523049,     1.8082846905287007,
        0.66036979659209594,    2.0302381468162172,      1.7474945080588224,     0.64022852658388707,
        2.0311237443757615,     1.6973217541912564,      0.62344371959497114,    2.0316260995995052,
        1.6583684918666544,     0.61036531370474667,     2.0320035484929773,     1.6298752221309873,
        0.60081647419738959,    2.0323847791259504,      1.6102287485388747,     0.59428198474366423,
        2.0328069817648231,     1.5974716090881289,      0.59009998633893201,    2.0332557496801353,
        1.5896968870698536,     0.58761206772364238,     2.0336968745365094,     1.5852845581175492,
        0.58625420696726294,    2.0340960929247465,      1.5829925836013352,     0.585592872052207,
        2.034427337542605,      1.5819455040901409,      0.58532191652256316,    2.0346721748722003,
        1.5815679853618057,     0.58523770158792721,     2.0348131238867224,     1.5815004648688928,
        0.58520585571231976,    2.034822222351,          1.5815187776147746,     0.58512704325752485,
        2.0346442681071437,     1.5814659878846526,      0.5849035286930615,     2.0341720308688842,
        1.5811952385729477,     0.58440415921797051,     2.033208509577066,      1.5805171336537609,
        0.58342259360747839,    2.0314089405271361,      1.5791427090243273,     0.58162176646410424,
        2.0281927729351548,     1.5766122241043838,      0.57845656176320692,    2.0226138102972846,
        1.5722006159975765,     0.5730673827288072,      2.0131773291137574,     1.564794309056396,
        0.56414300085195546,    1.9976021803466089,      1.552746785564848,      0.54976975670767492,
        1.9725576211698477,     1.5337541847315457,      0.52733183255660021,    1.9334865180719294,
        1.5048675107214002,     0.49362655873531291,     1.8747976907486499,     1.4628893135643182,
        0.44551095925032663,    1.7909728097300366,      1.4055295637036807,     0.3814732898383919,
        1.6792668399106168,     1.3335159347432795,      0.30407577224826471,    1.5439136297516218,
        1.2527940683764187,     0.22168569180494044,     1.3991485893485682,     1.1742542377091163,
        0.14645313476944921,    1.2659100110493404,      1.1090566212787396,     0.088111397461354629,
        1.1614224186453435,     1.062725248299383,       0.049143429181278217,   1.0907395715628008,
        1.0338394645033646,     0.025975136539660233,    1.0482321435664992,     1.0174715585005327,
        0.013250149964432939,   1.0246897182273387,      1.008760935216894,      0.006600566664965943,
        1.0123232752322964,     1.0043044379834345,      0.0032319460500384134,  1.0060402679496061,
        1.0020817646924807,     0.0015603889364410647,   1.0029177734288544,     1.0009931261867016,
        0.00074376823510747593, 1.0013911321108913,      1.0004676878282481,     0.00035011491520890698,
        1.0006549313856377,     1.0002174276763773,      0.00016273603349214043, 1.0003044353527186,
        1.0000997640913944,     7.4662410862154883e-05,  1.0001396774797513,     1.0000451622328106,
        3.3797135764305304e-05, 1.0000632287241606,      1.0000201629313166,     1.5087687516487832e-05,
        1.0000282284243927,     1.0000088755948804,      6.6389588405568395e-06, 1.0000124258953313,
        1.0000038538567255,     2.8761508570101123e-06,  1.0000053954112449,     1.0000016579572197,
        1.2204806728969253e-06, 1.0000023211422868,      1.0000007264076367,     4.921326823520219e-07,
        1.0000010169710787,     1.0000004246025982,      1.8885628667237185e-07, 1.0000005944437216}) {
    expected_end_state[counter] = element;
    ++counter;
  }
  EXPECT_LT((expected_end_state - time_stepper.solution().rbegin()->second.vector()).sup_norm(), 1e-15);
} // ... (inviscid_compressible_flow__euler_1d_explicit_fv, impermeable_walls_boundary_values__direct_euler_treatment)
