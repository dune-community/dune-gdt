// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH

#if HAVE_DUNE_XT_DATA

#  include <chrono>

#  include <dune/xt/common/parallel/threadmanager.hh>
#  include <dune/xt/common/string.hh>
#  include <dune/xt/test/gtest/gtest.h>

#  include <dune/xt/grid/information.hh>
#  include <dune/xt/grid/gridprovider.hh>

#  include <dune/xt/la/container.hh>

#  include <dune/gdt/discretefunction/default.hh>
#  include <dune/gdt/operators/advection-fv-entropybased.hh>
#  include <dune/gdt/operators/advection-fv.hh>
#  include <dune/gdt/interpolations/default.hh>
#  include <dune/gdt/test/momentmodels/entropysolver.hh>
#  include <dune/gdt/operators/rhs.hh>
#  include <dune/gdt/test/momentmodels/entropycalculator.hh>
#  include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#  include <dune/gdt/local/operators/advection-fv.hh>
#  include <dune/gdt/spaces/l2/finite-volume.hh>
#  include <dune/gdt/tools/timestepper/fractional-step.hh>
#  include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#  include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#  include <dune/gdt/test/momentmodels/kineticequation.hh>

#  include "pn-discretization.hh"

template <class GV, class ProblemType, class MapType>
class BoundaryFluxesFunctor : public Dune::XT::Grid::IntersectionFunctor<GV>
{
  using BaseType = typename Dune::XT::Grid::IntersectionFunctor<GV>;
  using IndexSetType = typename GV::IndexSet;

public:
  using typename BaseType::E;
  using typename BaseType::I;

  BoundaryFluxesFunctor(const ProblemType& problem, MapType& boundary_fluxes_map)
    : problem_(problem)
    , boundary_fluxes_map_(boundary_fluxes_map)
    , mutex_(std::make_shared<std::mutex>())
  {}

  BoundaryFluxesFunctor(const BoundaryFluxesFunctor& other)
    : BaseType(other)
    , problem_(other.problem_)
    , boundary_fluxes_map_(other.boundary_fluxes_map_)
    , mutex_(other.mutex_)
  {}

  Dune::XT::Grid::IntersectionFunctor<GV>* copy() override final
  {
    return new BoundaryFluxesFunctor(*this);
  }

  void apply_local(const I& intersection, const E& /*inside_element*/, const E& /*outside_element*/) override final
  {
    // store boundary fluxes
    const auto x = intersection.geometry().center();
    const auto dd = intersection.indexInInside() / 2;
    const auto n = intersection.centerUnitOuterNormal()[dd];
    auto boundary_flux = problem_.kinetic_boundary_flux(x, n, dd);
    mutex_->lock();
    boundary_fluxes_map_.insert(std::make_pair(x, boundary_flux));
    mutex_->unlock();
  }

private:
  const ProblemType& problem_;
  MapType& boundary_fluxes_map_;
  std::shared_ptr<std::mutex> mutex_;
};

template <class TestCaseType>
struct HyperbolicMnDiscretization
{
  // returns: (l1norm, l2norm, linfnorm, MPI rank)
  static std::pair<Dune::FieldVector<double, 3>, int> run(size_t num_save_steps = 1,
                                                          size_t num_output_steps = 0,
                                                          size_t quad_order = size_t(-1),
                                                          size_t quad_refinements = size_t(-1),
                                                          std::string grid_size = "",
                                                          size_t overlap_size = 2,
                                                          double t_end = 0.,
                                                          std::string filename = "",
                                                          bool disable_thread_cache = false)
  {
    using namespace Dune;
    using namespace Dune::GDT;

    //******************* get typedefs and constants from ProblemType **********************//
    using MomentBasis = typename TestCaseType::MomentBasis;
    using DiscreteFunctionType = typename TestCaseType::DiscreteFunctionType;
    using GridType = typename TestCaseType::GridType;
    using SpaceType = typename TestCaseType::SpaceType;
    using AdvectionSourceSpaceType = typename TestCaseType::AdvectionSourceSpaceType;
    using GV = typename TestCaseType::GridViewType;
    using I = XT::Grid::extract_intersection_t<GV>;
    using ProblemType = typename TestCaseType::ProblemType;
    using RangeFieldType = typename MomentBasis::RangeFieldType;
    using BoundaryValueType = typename ProblemType::BoundaryValueType;
    static constexpr size_t dimDomain = MomentBasis::dimDomain;
    static constexpr size_t dimRange = MomentBasis::dimRange;
    static constexpr auto la_backend = TestCaseType::la_backend;
    using DomainType = FieldVector<RangeFieldType, dimDomain>;
    using RangeType = FieldVector<RangeFieldType, dimRange>;
    using MatrixType = typename XT::LA::Container<RangeFieldType, la_backend>::MatrixType;
    using VectorType = typename XT::LA::Container<RangeFieldType, la_backend>::VectorType;
    using DynamicRangeType = DynamicVector<RangeFieldType>;

    //******************* create grid and FV space ***************************************
    auto grid_config = ProblemType::default_grid_cfg();
    if (!grid_size.empty())
      grid_config["num_elements"] = grid_size;
    grid_config["overlap_size"] = XT::Common::to_string(overlap_size);
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GV grid_view(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_view);
    const AdvectionSourceSpaceType advection_source_space(grid_view);

    //******************* create EquationType object ***************************************
    std::shared_ptr<const MomentBasis> basis_functions = std::make_shared<const MomentBasis>(
        quad_order == size_t(-1) ? MomentBasis::default_quad_order() : quad_order,
        quad_refinements == size_t(-1) ? MomentBasis::default_quad_refinements() : quad_refinements);
    const RangeFieldType psi_vac = DXTC_CONFIG_GET("psi_vac", 1e-8 / basis_functions->unit_ball_volume());
    const double tau = DXTC_CONFIG_GET("opt_tau", 1e-9);
    const std::unique_ptr<ProblemType> problem_ptr =
        std::make_unique<ProblemType>(*basis_functions, grid_view, psi_vac, grid_config, false, tau);
    const auto& problem = *problem_ptr;
    const auto initial_values = problem.initial_values();
    const auto boundary_values = problem.boundary_values();
    using AnalyticalFluxType = typename ProblemType::FluxType;
    using EntropyFluxType = typename ProblemType::ActualFluxType;
    auto analytical_flux = problem.flux();
    auto* entropy_flux = dynamic_cast<EntropyFluxType*>(analytical_flux.get());
    // for Legendre polynomials and real spherical harmonics, the results are sensitive to the initial guess in the
    // Newton algorithm. If the thread cache is enabled, the guess is different dependent on how many threads we are
    // using, so for the tests we disable this cache to get reproducible results.
    if (disable_thread_cache)
      dynamic_cast<EntropyFluxType*>(analytical_flux.get())->disable_thread_cache();
    // const RangeFieldType CFL = DXTC_CONFIG.get("timestepper.CFL", problem.CFL());

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    VectorType u_vec(fv_space.mapper().size(), 0., DXTC_CONFIG_GET("num_u_mutexes", 1));
    DiscreteFunctionType u(fv_space, u_vec, "rho");
    // project initial values
    default_interpolation(*initial_values, u, grid_view);

    // ******************** choose flux and rhs operator and timestepper ******************************************

    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using EigenvectorWrapperType = typename EigenvectorWrapperChooser<MomentBasis, AnalyticalFluxType>::type;
    using EntropySolverType = EntropySolver<MomentBasis, SpaceType, MatrixType>;
    using EntropyCalculatorType = EntropyCalculator<MomentBasis, SpaceType, MatrixType>;
    using ReconstructionOperatorType = PointwiseLinearReconstructionOperator<AnalyticalFluxType,
                                                                             BoundaryValueType,
                                                                             GV,
                                                                             VectorType,
                                                                             EigenvectorWrapperType>;
    using ReconstructionAdvectionOperatorType =
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
    // using FvOperatorType = EntropyBasedMomentFvOperator<
    //     std::conditional_t<TestCaseType::reconstruction,
    //                        ReconstructionAdvectionOperatorType,
    //                        AdvectionOperatorType>,
    //     EntropySolverType,
    //     EntropyCalculatorType>;
    constexpr TimeStepperMethods time_stepper_type = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
    using RhsOperatorType = RhsOperator<MomentBasis, SpaceType, MatrixType>;
    using FvOperatorType = EntropyBasedMomentCombinedFvOperator<
        std::conditional_t<TestCaseType::reconstruction, ReconstructionAdvectionOperatorType, AdvectionOperatorType>,
        EntropySolverType,
        EntropyCalculatorType,
        RhsOperatorType>;
    using TimeStepperType = ExplicitRungeKuttaTimeStepper<FvOperatorType, DiscreteFunctionType, time_stepper_type>;
    // using OperatorTimeStepperType =
    //     ExplicitRungeKuttaTimeStepper<FvOperatorType, DiscreteFunctionType, time_stepper_type>;
    // using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, MomentBasis>;
    // using TimeStepperType = StrangSplittingTimeStepper<RhsTimeStepperType, OperatorTimeStepperType>;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if constexpr (dimDomain == 2)
      dx /= std::sqrt(2);
    else if constexpr (dimDomain == 3)
      dx /= std::sqrt(3);
    // RangeFieldType dt = CFL * dx;
    const RangeFieldType dt = DXTC_CONFIG.get("timestepper.dt", problem.CFL() * dx);

    // *********************** create operators and timesteppers ************************************
    NumericalKineticFlux<GV, MomentBasis> numerical_flux(*analytical_flux, *basis_functions);
    AdvectionOperatorType advection_operator(grid_view, numerical_flux, advection_source_space, fv_space, true);

    // boundary treatment
    using BoundaryOperator =
        // LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
        LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator<I, VectorType, GV, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;

    // store boundary fluxes
    using BoundaryFluxesMapType = std::map<DomainType, DynamicRangeType, XT::Common::FieldVectorFloatLess>;
    BoundaryFluxesMapType boundary_fluxes;
    BoundaryFluxesFunctor<GV, ProblemType, BoundaryFluxesMapType> boundary_flux_functor(problem, boundary_fluxes);
    auto walker = XT::Grid::Walker<GV>(fv_space.grid_view());
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> boundary_intersection_filter;
    walker.append(boundary_flux_functor, boundary_intersection_filter);
    walker.walk(true);
    using GenericFunctionType = XT::Functions::GenericFunction<dimDomain, dimRange, 1, RangeFieldType>;
    GenericFunctionType boundary_kinetic_fluxes(
        1, [&](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) { ret = boundary_fluxes[x]; });
    LambdaType boundary_lambda =
        [&boundary_kinetic_fluxes,
         &entropy_flux](const I& intersection,
                        const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
                        const DynamicRangeType& u_in,
                        DynamicRangeType& g,
                        const XT::Common::Parameter& param) {
          // influx
          // boundary_kinetic flux is caclulated as <(v*n_outside) b psi>, so it should have a positive sign in the
          // timestepper update. As the advection operator is added to the timestepper with negative sign, we need to
          // multiply by -1 here.
          boundary_kinetic_fluxes.evaluate(
              intersection.geometry().global(xx_in_reference_intersection_coordinates), g, param);
          g *= -1.;
          // outflux
          const auto& entity = intersection.inside();
          const auto dd = intersection.indexInInside() / 2;
          const auto alpha_entity = entropy_flux->get_alpha(entity, u_in, true)->first;
          const auto outflux =
              entropy_flux->evaluate_kinetic_outflow(alpha_entity, intersection.centerUnitOuterNormal(), dd);
          // outflux is calculated as <(v*n) b psi>, so it should have a negative sign in the timestepper update, which
          // is why it has a positive sign here (see previous comment).
          g += outflux;
        };
    advection_operator.append(boundary_lambda, {}, boundary_intersection_filter);

    constexpr double epsilon = 1e-11;
    auto slope = TestCaseType::RealizabilityLimiterChooserType::template make_slope<EigenvectorWrapperType>(
        *entropy_flux, *basis_functions, epsilon);
    ReconstructionOperatorType reconstruction_operator(*analytical_flux, *boundary_values, fv_space, *slope, false);
    ReconstructionAdvectionOperatorType reconstruction_advection_operator(advection_operator, reconstruction_operator);

    if (XT::Common::is_zero(t_end))
      t_end = problem.t_end();

    if (!filename.empty())
      filename += "_";
    filename += ProblemType::static_id();
    filename +=
        (time_stepper_type == TimeStepperMethods::explicit_rungekutta_second_order_ssp)
            ? "_ssp2"
            : (time_stepper_type == TimeStepperMethods::explicit_rungekutta_third_order_ssp ? "_ssp3" : "_unknown");
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_dt_" + XT::Common::to_string(dt);
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad_" + XT::Common::to_string(quad_refinements) + "x" + XT::Common::to_string(quad_order);
    filename += "_threads_" + DXTC_CONFIG.get("threading.max_count", "1") + "x"
                + DXTC_CONFIG.get("threading.partition_factor", "1");
    filename += MomentBasis::entropy == EntropyType::MaxwellBoltzmann ? "_MaxwellBoltzmann" : "_BoseEinstein";
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->mn_name();

    EntropySolverType entropy_solver(*(dynamic_cast<EntropyFluxType*>(analytical_flux.get())),
                                     fv_space,
                                     problem.psi_vac() * basis_functions->unit_ball_volume() / 10,
                                     filename);
    EntropyCalculatorType entropy_calculator(*(dynamic_cast<EntropyFluxType*>(analytical_flux.get())), fv_space);

    // static const RangeType u_iso = basis_functions->u_iso();
    // static const RangeType basis_integrated = basis_functions->integrated();
    const auto sigma_a = problem.sigma_a();
    const auto sigma_s = problem.sigma_s();
    const auto Q = problem.Q();
    // auto rhs_func = [&](const auto& source,
    //                     const auto& local_source,
    //                     auto& local_range,
    //                     const Dune::XT::Common::Parameter& /*param*/) {
    //   const auto& element = local_range.element();
    //   local_source[0]->bind(element);
    //   const auto center = element.geometry().center();
    //   // const auto& u_elem = analytical_flux->get_precomputed_u(fv_space.grid_view().indexSet().index(element));
    //   const auto u_elem = local_source[0]->evaluate(center);
    //   const auto sigma_a_value = sigma_a->evaluate(center)[0];
    //   const auto sigma_s_value = sigma_s->evaluate(center)[0];
    //   const auto sigma_t_value = sigma_a_value + sigma_s_value;
    //   const auto Q_value = Q->evaluate(center)[0];
    //   auto ret = u_elem;
    //   ret *= -sigma_t_value;
    //   ret.axpy(basis_functions->density(u_elem) * sigma_s_value, u_iso);
    //   ret.axpy(Q_value, basis_integrated);
    //   auto& range_dofs = local_range.dofs();
    //   for (size_t ii = 0; ii < dimRange; ++ii)
    //     range_dofs.add_to_entry(ii, ret[ii]);
    // };
    RhsOperatorType rhs_operator(
        *(dynamic_cast<EntropyFluxType*>(analytical_flux.get())), fv_space, *sigma_a, *sigma_s, *Q);

    // rhs_operator.append(GenericLocalElementOperator<VectorType, GV, dimRange>(rhs_func, 1));
    FvOperatorType fv_operator(
        FvOperatorChooser<TestCaseType::reconstruction>::choose(advection_operator, reconstruction_advection_operator),
        entropy_solver,
        entropy_calculator,
        rhs_operator);

    // ******************************** do the time steps ***********************************************************
    // const auto sigma_a = problem.sigma_a();
    // const auto sigma_s = problem.sigma_s();
    // const auto Q = problem.Q();
    // OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    TimeStepperType timestepper(fv_operator, u);
    // RhsTimeStepperType timestepper_rhs(*basis_functions, u, *sigma_a, *sigma_s, *Q);
    // TimeStepperType timestepper(timestepper_rhs, timestepper_op);

    auto begin_time = std::chrono::steady_clock::now();
    timestepper.solve(t_end,
                      dt,
                      num_save_steps,
                      num_output_steps,
                      false,
                      DXTC_CONFIG_GET("visualize", true),
                      DXTC_CONFIG_GET("write_txt", true),
                      false,
                      true,
                      filename,
                      *basis_functions->visualizer(),
                      basis_functions->stringifier());
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (grid_view.comm().rank() == 0)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;
    timestepper.write_timings(filename);

    auto ret = std::make_pair(FieldVector<double, 3>(0.), int(0));
    double& l1norm = ret.first[0];
    double& l2norm = ret.first[1];
    double& linfnorm = ret.first[2];
    ret.second = grid_view.comm().rank();
    const auto& current_sol = timestepper.current_solution();
    const auto local_sol = current_sol.local_function();
    for (const auto& entity : elements(grid_view, Dune::Partitions::interior)) {
      local_sol->bind(entity);
      const auto val = local_sol->evaluate(entity.geometry().local(entity.geometry().center()));
      RangeFieldType psi = basis_functions->density(val);
      l1norm += std::abs(psi) * entity.geometry().volume();
      l2norm += std::pow(psi, 2) * entity.geometry().volume();
      linfnorm = std::max(std::abs(psi), linfnorm);
    }
    l1norm = grid_view.comm().sum(l1norm);
    l2norm = grid_view.comm().sum(l2norm);
    linfnorm = grid_view.comm().max(linfnorm);
    l2norm = std::sqrt(l2norm);
    return ret;
  }
};

template <class TestCaseType>
struct HyperbolicMnTest
  : public HyperbolicMnDiscretization<TestCaseType>
  , public ::testing::Test
{
  void run(const double tol = TestCaseType::ExpectedResultsType::tol)
  {
    auto norms = HyperbolicMnDiscretization<TestCaseType>::run(
                     DXTC_CONFIG.get("num_save_steps", 1),
                     DXTC_CONFIG.get("num_output_steps", 0),
                     DXTC_CONFIG.get("quad_order", TestCaseType::quad_order),
                     DXTC_CONFIG.get("quad_refinements", TestCaseType::quad_refinements),
                     DXTC_CONFIG.get("grid_size", ""),
                     DXTC_CONFIG.get("overlap_size", 2),
                     DXTC_CONFIG.get("t_end", TestCaseType::t_end),
                     DXTC_CONFIG.get("filename", "timings"),
                     DXTC_CONFIG.get("disable_thread_cache",
                                     Dune::GDT::is_full_moment_basis<typename TestCaseType::MomentBasis>::value))
                     .first;
    const double l1norm = norms[0];
    const double l2norm = norms[1];
    const double linfnorm = norms[2];
    using ResultsType = typename TestCaseType::ExpectedResultsType;
    EXPECT_NEAR(ResultsType::l1norm, l1norm, ResultsType::l1norm * tol);
    EXPECT_NEAR(ResultsType::l2norm, l2norm, ResultsType::l2norm * tol);
    EXPECT_NEAR(ResultsType::linfnorm, linfnorm, ResultsType::linfnorm * tol);
  }
};

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
