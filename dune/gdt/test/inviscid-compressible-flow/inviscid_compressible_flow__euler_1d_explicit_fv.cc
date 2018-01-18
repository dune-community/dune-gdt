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

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

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


class InviscidCompressibleFlow1dExplicitEuler : public ::testing::Test
{
protected:
  using G = YASP_1D_EQUIDISTANT_OFFSET;
  using E = typename G::template Codim<0>::Entity;
  using D = double;
  static const constexpr size_t d = G::dimension;
  using R = double;
  static const constexpr size_t m = d + 2;

  using DomainType = XT::Common::FieldVector<D, d>;
  using RangeType = XT::Common::FieldVector<D, m>;

  using GridProvider = XT::Grid::GridProvider<G>;
  using LeafGL = typename GridProvider::LeafGridViewType;
  using GL = XT::Grid::PeriodicGridLayer<LeafGL>;

  using UI = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>;
  using U = XT::Functions::GlobalLambdaFunction<E, D, d, R, m>;
  using F = XT::Functions::GlobalLambdaFluxFunction<UI, 0, R, d, m>;

  using NF = NumericalVijayasundaramFlux<E, D, d, R, m>;

  using S = FvSpace<GL, R, m>;
  using V = XT::LA::EigenDenseVector<R>;
  using DF = DiscreteFunction<S, V>;
  using FvOperator = GDT::AdvectionFvOperator<DF>;

  using TimeStepper = ExplicitRungeKuttaTimeStepper<FvOperator, DF, TimeStepperMethods::explicit_euler>;

  InviscidCompressibleFlow1dExplicitEuler()
    : logger_(XT::Common::TimedLogger().get("main"))
    , euler_tools_(/*gamma=*/1.4) // air or water at roughly 20 deg Cels.
    , grid_(nullptr)
    , leaf_grid_layer_(nullptr)
    , grid_layer_(nullptr)
    , space_(nullptr)
    , initial_values_(nullptr)
    , flux_(nullptr)
    , numerical_flux_(nullptr)
    , fv_op_(nullptr)
  {
  }

  void SetUp() override final
  {
    grid_ = std::make_shared<GridProvider>(XT::Grid::make_cube_grid<G>(-1., 1., 128u));
    leaf_grid_layer_ = std::make_shared<LeafGL>(grid_->leaf_view());
    grid_layer_ = std::make_shared<GL>(XT::Grid::make_periodic_grid_layer(*leaf_grid_layer_));
    logger_.info() << "grid has " << leaf_grid_layer_->indexSet().size(0) << " elements" << std::endl;

    space_ = std::make_shared<S>(*grid_layer_);
    logger_.info() << "space has " << space_->mapper().size() << " DoFs\n" << std::endl;

    initial_values_ = std::make_shared<DF>(*space_, "solution");
    const U periodic_initial_values_euler( // compare [Kröner, 1997, p. 394]
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
          return euler_tools_.to_conservative(primitive_variables);
        },
        /*order=*/0,
        /*parameter_type=*/{},
        /*name=*/"periodic_initial_values_euler");
    project(periodic_initial_values_euler, *initial_values_);

    flux_ = std::shared_ptr<F>(new F([&](const auto& /*x*/,
                                         const auto& conservative_variables,
                                         const auto& /*mu*/) { return euler_tools_.flux(conservative_variables); },
                                     {},
                                     "euler_flux",
                                     [](const auto& /*mu*/) { return 3; },
                                     [&](const auto& /*x*/, const auto& conservative_variables, const auto& /*mu*/) {
                                       return euler_tools_.flux_jacobian(conservative_variables);
                                     }));
    numerical_flux_ = std::make_shared<NF>(*flux_);
  } // ... SetUp(...)

  void do_the_timestepping(
      const size_t num_timesteps = 100,
      const std::pair<XT::Common::FieldVector<R, m>, XT::Common::FieldVector<R, m>>& boundary_data_range = {
          XT::Common::FieldVector<R, m>(std::numeric_limits<R>::max()),
          XT::Common::FieldVector<R, m>(std::numeric_limits<R>::min())})
  {
    ASSERT_NE(fv_op_, nullptr) << "You have to create the fv_op_ before calling this method!";

    time_stepper_ = std::make_shared<TimeStepper>(*fv_op_, *initial_values_, -1.);

    // no idea why we need to provide the <GL, E, D, d, R, m> here, but the compiler could not figure it out without
    const auto dt = estimate_dt_for_hyperbolic_system<GL, E, D, d, R, m>(
        *grid_layer_, *initial_values_, *flux_, boundary_data_range);
    EXPECT_DOUBLE_EQ(0.002708880865541605, dt);
    const double T = num_timesteps * dt;

    // test for stability
    const auto test_dt =
        time_stepper_->find_suitable_dt(dt, 10, 1.1 * std::max(T, 1.) * initial_values_->vector().sup_norm(), 25);
    ASSERT_TRUE(test_dt.first) << "Could not determine optimal dt (in particular, the dt computed to match the CFL "
                                  "condition did not yield a stable scheme)!";
    ASSERT_LE(test_dt.second, dt) << "The computed dt (to match the CFL condition) does not yield a stable scheme: "
                                  << dt << "\nThe following dt seems to work fine: " << test_dt.second;

    // do the work
    time_stepper_->solve(T,
                         dt,
                         num_timesteps,
                         /*save_solution=*/true,
                         /*visualize=*/false,
                         "solution",
                         /*visualizer=*/[&](const auto& u, const auto& filename_prefix, const auto& step) {
                           euler_tools_.visualize(u, *grid_layer_, filename_prefix, XT::Common::to_string(step));
                         });
  } // ... do_the_timestepping(...)

  XT::Common::TimedLogManager logger_;
  EulerTools<d> euler_tools_;
  std::shared_ptr<GridProvider> grid_;
  std::shared_ptr<LeafGL> leaf_grid_layer_;
  std::shared_ptr<GL> grid_layer_;
  std::shared_ptr<S> space_;
  std::shared_ptr<DF> initial_values_;
  std::shared_ptr<F> flux_;
  std::shared_ptr<NF> numerical_flux_;
  std::shared_ptr<FvOperator> fv_op_;
  std::shared_ptr<TimeStepper> time_stepper_;
}; // class InviscidCompressibleFlow1dExplicitEuler


TEST_F(InviscidCompressibleFlow1dExplicitEuler, periodic_boundaries)
{
  ASSERT_NE(grid_layer_, nullptr);
  ASSERT_NE(numerical_flux_, nullptr);
  ASSERT_NE(initial_values_, nullptr);
  ASSERT_NE(space_, nullptr);

  // create operator, the layer is periodic and the operator includes handling of periodic boundaries
  fv_op_ = std::make_shared<FvOperator>(*grid_layer_, *numerical_flux_);

  // do timestepping
  const size_t num_timesteps = 100;
  this->do_the_timestepping(num_timesteps);
  ASSERT_NE(time_stepper_, nullptr);

  // check conservation principle
  const auto initial_l_infty_norm = initial_values_->vector().sup_norm();
  const auto initial_l_2_norm = make_l2_operator(*grid_layer_)->induced_norm(*initial_values_);
  double infty = 0.;
  double l2 = 0.;
  size_t saved_timesteps = 0;
  for (const auto& t_and_u : time_stepper_->solution()) {
    ++saved_timesteps;
    const auto& u_conservative = t_and_u.second;
    const auto relative_l_infty_error =
        std::abs(initial_l_infty_norm - u_conservative.vector().sup_norm()) / initial_l_infty_norm;
    EXPECT_LT(relative_l_infty_error, 0.065554) << "initial_l_infty_norm = " << initial_l_infty_norm
                                                << "\nl_infty_norm = " << u_conservative.vector().sup_norm();
    const auto relative_l_2_error =
        std::abs(initial_l_2_norm - make_l2_operator(*grid_layer_)->induced_norm(u_conservative)) / initial_l_2_norm;
    EXPECT_LT(relative_l_2_error, 0.10577)
        << "initial_l_2_norm = " << initial_l_2_norm
        << "\nl_2_norm = " << make_l2_operator(*grid_layer_)->induced_norm(u_conservative);
    infty = std::max(infty, relative_l_infty_error);
    l2 = std::max(l2, relative_l_2_error);
  }
  EXPECT_EQ(num_timesteps + 1, saved_timesteps);

  // check expected state at the end
  V expected_end_state(space_->mapper().size());
  size_t counter = 0;
  for (const auto& value :
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
    expected_end_state[counter] = value;
    ++counter;
  }

  const auto& actual_end_state = time_stepper_->solution().rbegin()->second.vector();
  EXPECT_LT((expected_end_state - actual_end_state).sup_norm(), 1e-15);

  if ((expected_end_state - actual_end_state).sup_norm() > 1e-15) {
    std::cout << "actual_end_state = {" << std::setprecision(17) << actual_end_state[0];
    for (size_t ii = 1; ii < actual_end_state.size(); ++ii)
      std::cout << ", " << std::setprecision(17) << actual_end_state[ii];
    std::cout << "}" << std::endl;
  }
} // ... (inviscid_compressible_flow__euler_1d_explicit_fv, periodic_boundaries)


TEST_F(InviscidCompressibleFlow1dExplicitEuler, impermeable_walls_by_direct_euler_treatment)
{
  ASSERT_NE(grid_layer_, nullptr);
  ASSERT_NE(numerical_flux_, nullptr);
  ASSERT_NE(space_, nullptr);

  // impermeable walls everywhere
  XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GL>> boundary_info;
  boundary_info.register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
  boundary_info.register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
  XT::Grid::ApplyOn::CustomBoundaryIntersections<GL> impermeable_wall_filter(boundary_info,
                                                                             new XT::Grid::ImpermeableBoundary());

  // create operator, the layer is periodic and the operator includes handling of periodic boundaries so we need to make
  // an exception for all non-periodic boundaries
  fv_op_ = std::make_shared<FvOperator>(
      *grid_layer_, *numerical_flux_, /*periodicity_restriction=*/impermeable_wall_filter.copy());

  // the actual handling of impermeable walls
  const auto euler_impermeable_wall_treatment = [&](const auto& u, const auto& n, const auto& /*mu*/ = {}) {
    return euler_tools_.flux_at_impermeable_walls(u, n);
  };
  fv_op_->append(euler_impermeable_wall_treatment, impermeable_wall_filter.copy(), {});

  // do timestepping
  this->do_the_timestepping(300); // we need enough for the wave to hit the boundary
  ASSERT_NE(time_stepper_, nullptr);

  // check expected state at the end
  V expected_end_state(space_->mapper().size());
  size_t counter = 0;
  for (const auto& value :
       {2.4333683449209049,      -0.00034510309372732661, 3.5070282077224091,      2.4129251526029525,
        -0.00048135403348335907, 3.5064729726900739,      2.4102795368359966,      -0.00077585929697480613,
        3.5061597418367745,      2.4094606046399374,      -0.00099329187530774587, 3.5055297521690312,
        2.408558773860435,       -0.0010958271868925979,  3.504367904647514,       2.4079702252657471,
        -0.0011053391961194706,  3.5022371280433693,      2.4089949999331002,      -0.0012221659306074822,
        3.4982136293279811,      2.4144451797505933,      -0.0020358236847386113,  3.4904241519559549,
        2.4288615756272685,      -0.004890008933088808,   3.475270693843457,       2.4570330768754882,
        -0.012550639873116673,   3.4461037207420779,      2.4995299909813693,      -0.03040983919076598,
        3.3910283704071111,      2.545553822030584,       -0.068129226165393839,   3.2902004237202367,
        2.5675244687087826,      -0.13964961033183751,    3.1164469198428,         2.5285593894752196,
        -0.25493861764348091,    2.8508391252368064,      2.4150802546676733,      -0.39909657415126193,
        2.5209846950326793,      2.2713606921342615,      -0.5262694610396732,     2.2151842769994352,
        2.1655306565307506,      -0.60216637594507949,    2.0076136726529508,      2.1168503403819887,
        -0.6326422784952892,     1.8933056527059249,      2.1011997694768874,      -0.63735295234516065,
        1.8297943004937052,      2.0945710613215969,      -0.62906248695577127,    1.7863293594080167,
        2.0857058880154495,      -0.61373043952860029,    1.7491365436277015,      2.0714322799315861,
        -0.59414717978247433,    1.7132625078929753,      2.0519435595307107,      -0.57184194172773084,
        1.6772066823679197,      2.0284510857141078,      -0.54779231362320702,    1.6406654059546013,
        2.0022197595551376,      -0.52266703481454146,    1.6037089656655776,      1.9742419813467598,
        -0.49692548003110443,    1.566501284652392,       1.9451926445596064,      -0.47087983805797512,
        1.5292123552223786,      1.915493545352245,       -0.44474440077580701,    1.4919961415081806,
        1.8854012268195584,      -0.41867335506936448,    1.4549886554763136,      1.8550795598099479,
        -0.39278622013675712,    1.4183107146625968,      1.8246466083646666,      -0.36718262807007956,
        1.3820707893382493,      1.7941997254488722,      -0.34194983393013972,    1.3463671297693427,
        1.7638271394429743,      -0.31716623881434702,    1.311289426249566,       1.7336131602499452,
        -0.29290308972833917,    1.2769204691284517,      1.703641709842816,       -0.26922519939177975,
        1.2433385572909916,      1.6740014583021761,      -0.24619035497768835,    1.210622083972672,
        1.6447955005133721,      -0.22384622586405092,    1.1788584407099003,      1.6161577304704366,
        -0.20222346463070379,    1.148159033676539,       1.5882752079434734,      -0.18132492809849177,
        1.1186796971263935,      1.5614108616406244,      -0.16111357127495607,    1.0906413574393254,
        1.5359163975388874,      -0.14150435455128185,    1.0643418654422354,      1.5122251942948925,
        -0.1223662545549494,     1.0401499993200831,      1.4908213721789294,      -0.10353762270221029,
        1.0184785176214364,      1.4721916932478558,      -0.084852500073558892,   0.99974244214126518,
        1.4567752628956356,      -0.066170219214407883,   0.98431582661532613,     1.4449266161929581,
        -0.047398911979264186,   0.97250039752200312,     1.4369003145741166,      -0.028506126355068349,
        0.96451269954371399,     1.4328546002007283,      -0.0095145744070903714,  0.9604873105835714,
        1.4328546002007325,      0.009514574407090675,    0.9604873105835714,      1.4369003145741146,
        0.028506126355068766,    0.96451269954371377,     1.4449266161929539,      0.047398911979264657,
        0.9725003975220019,      1.4567752628956325,      0.066170219214409576,    0.98431582661532335,
        1.4721916932478503,      0.084852500073564555,    0.99974244214125618,     1.490821372178899,
        0.10353762270222792,     1.0184785176214057,      1.512225194294782,       0.12236625455500889,
        1.0401499993199799,      1.5359163975385091,      0.14150435455148053,     1.0643418654418886,
        1.5614108616393474,      0.1611135712756131,      1.0906413574381595,      1.5882752079392146,
        0.18132492810064177,     1.1186796971225037,      1.6161577304563641,      0.20222346463767552,
        1.1481590336636784,      1.6447955004673245,      0.22384622588645015,     1.1788584406677711,
        1.6740014581529219,      0.24619035504896597,     1.2106220838359172,      1.7036417093635885,
        0.26922519961639152,     1.2433385568512121,      1.7336131587260069,      0.29290309042918505,
        1.2769204677276029,      1.7638271346440306,      0.31716624097946605,     1.3112894218302684,
        1.7941997104850789,      0.34194984055225808,     1.3463671159628441,      1.8246465621634915,
        0.36718264812361628,     1.3820707466233753,      1.8550794185474877,      0.39278628027346119,
        1.4183105837784973,      1.8854007989835124,      0.4186735337144683,      1.4549882581721825,
        1.9154922612340448,      0.44474492677993327,     1.4919949461572817,      1.9451888223846192,
        0.47088137443035638,     1.5292087880342873,      1.9742306885841578,      0.49692993674266756,
        1.5664907147753506,      2.0021866029676536,      0.52267989329160147,     1.6036778255569903,
        2.0283542225777103,      0.54782927871257614,     1.6405740400180842,      2.0516616845984972,
        0.57194802533474942,     1.6769392504001752,      2.0706146771479359,      0.59445162210373226,
        1.7124804303362862,      2.0833433030752762,      0.61460498209043968,     1.7468502116580378,
        2.0877908550522872,      0.6315751212443349,      1.7796588848745916,      2.0820477011888592,
        0.64454001106559577,     1.8105033846995473,      2.0647766964637757,      0.65284091835324076,
        1.8390115350580269,      2.0356239339363364,      0.65614373733121301,     1.8648933161555918,
        1.9954866654647552,      0.65456165075799555,     1.8879861019963664,      1.9465357482767609,
        0.6486953773896571,      1.9082804169921623,      1.8919591170249948,      0.63956944291929874,
        1.9259176447817434,      1.8354800777369022,      0.62847548695630706,     1.9411596011233374,
        1.7807742366561365,      0.61676318105153771,     1.9543384793000562,      1.7309339171600517,
        0.60563379128271044,     1.9658009274821626,      1.6881022776196775,      0.59598601418038855,
        1.9758601677545682,      1.653336646552704,       0.58834266097562049,     1.9847656986446407,
        1.6266897183133782,      0.58286008348576501,     1.9926935391637899,      1.6074440401041041,
        0.57940029251328307,     1.9997538856217738,      1.5944131392372038,      0.57763478274167401,
        2.0060093804271775,      1.5862305348362302,      0.57714981573575119,     2.0114964897594922,
        1.5815745779105104,      0.57753171576491036,     2.0162441042029728,      1.5793086405382393,
        0.57842238608037566,     2.0202861661488911,      1.5785418577822576,      0.57954535995969447,
        2.0236677039715785,      1.5786304030550515,      0.58070894346304147,     2.0264453977113899,
        1.5791433940799546,      0.58179515641368451,     2.0286845364711268,      1.5798142914671072,
        0.58274233387872365,     2.0304541522231663,      1.5804920797948037,      0.58352691590927153,
        2.0318215666389632,      1.5810997564173972,      0.58414737071329215,     2.0328468625312603,
        1.5816024027828757,      0.58461107275243229,     2.0335770575995453,      1.5819838188379787,
        0.5849234982603817,      2.0340390516754789,      1.5822290323571186,      0.58507818911100751,
        2.0342296791351311,      1.5823093280069866,      0.58504531697685225,     2.0341003079172473,
        1.5821661711523807,      0.58475608688342984,     2.033532239829178,       1.5816900793704474,
        0.58407943795580719,     2.0322975357876487,      1.5806898552438733,      0.58278640639662316,
        2.0299977239188691,      1.5788465796474089,      0.58049617000928955,     2.0259702948184266,
        1.5756456709022213,      0.57659670615040948,     2.019150830453921,       1.5702802619732577,
        0.57013383442315912,     2.0078799130690439,      1.5615232347177228,      0.5596695850393365,
        1.989656198863706,       1.5475808825588477,      0.54313439073904335,     1.9608775723847265,
        1.5259834472808385,      0.51775661011432761,     1.9167140540061662,      1.4936595833933037,
        0.48027107707924527,     1.8514600879017977,      1.4474906879906702,      0.42777099074502567,
        1.7599943997875562,      1.3857408046189987,      0.35957608961084353,     1.6409907446919481,
        1.3103950052215085,      0.27977329745970625,     1.5012969190500374,      1.2290385725359598,
        0.19825491546884244,     1.3577960828638163,      1.1533690491311199,      0.12726326486794629,
        1.231694946632331,       1.0934730078914803,      0.074729785223835637,    1.1372574100302058,
        1.0526935076165664,      0.040982819811075981,    1.0758204891008185,      1.0280828924623446,
        0.021459836039747677,    1.0399010326222926,      1.0144264454398266,      0.010910154434093764,
        1.0203518834346268,      1.0072412807218425,      0.00543914677486568,     1.010176875805129,
        1.0035853778734434,      0.0026687135362074691,   1.0050290806497879,      1.001772808907248,
        0.0012808311162860191,   1.0024842114992132,      1.0009088783081277,      0.00057737196828583782,
        1.0012729617122325,      1.0006280611243779,      0.00024873614355938508,  1.0008794416251166}) {
    expected_end_state[counter] = value;
    ++counter;
  }

  const auto& actual_end_state = time_stepper_->solution().rbegin()->second.vector();
  EXPECT_LT((expected_end_state - actual_end_state).sup_norm(), 1e-15);

  if ((expected_end_state - actual_end_state).sup_norm() > 1e-15) {
    std::cout << "actual_end_state = {" << std::setprecision(17) << actual_end_state[0];
    for (size_t ii = 1; ii < actual_end_state.size(); ++ii)
      std::cout << ", " << std::setprecision(17) << actual_end_state[ii];
    std::cout << "}" << std::endl;
  }
} // ... (inviscid_compressible_flow__euler_1d_explicit_fv, impermeable_walls_by_direct_euler_treatment)
