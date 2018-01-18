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
