// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH

#include <chrono>

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/advection-fv-entropybased.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

#include "pn-discretization.hh"

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
    static const auto la_backend = TestCaseType::la_backend;
    using MatrixType = typename XT::LA::Container<RangeFieldType, la_backend>::MatrixType;
    using VectorType = typename XT::LA::Container<RangeFieldType, la_backend>::VectorType;

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
    const std::unique_ptr<ProblemType> problem_ptr =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_view, grid_config);
    const auto& problem = *problem_ptr;
    const auto initial_values = problem.initial_values();
    const auto boundary_values = problem.boundary_values();
    using AnalyticalFluxType = typename ProblemType::FluxType;
    using EntropyFluxType = typename ProblemType::ActualFluxType;
    auto analytical_flux = problem.flux();
    // for Legendre polynomials and real spherical harmonics, the results are sensitive to the initial guess in the
    // Newton algorithm. If the thread cache is enabled, the guess is different dependent on how many threads we are
    // using, so for the tests we disable this cache to get reproducible results.
    if (disable_thread_cache)
      dynamic_cast<EntropyFluxType*>(analytical_flux.get())->disable_thread_cache();
    const RangeFieldType CFL = problem.CFL();

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "solution");
    // project initial values
    default_interpolation(*initial_values, u, grid_view);

    // ******************** choose flux and rhs operator and timestepper ******************************************

    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using EigenvectorWrapperType = typename EigenvectorWrapperChooser<MomentBasis, AnalyticalFluxType>::type;
    using EntropySolverType = EntropySolver<MomentBasis, SpaceType, MatrixType>;
    //    using ReconstructionOperatorType =
    //        LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, GV, MatrixType,
    //        EigenvectorWrapperType>;
    using ReconstructionOperatorType = PointwiseLinearReconstructionOperator<AnalyticalFluxType,
                                                                             BoundaryValueType,
                                                                             GV,
                                                                             VectorType,
                                                                             EigenvectorWrapperType>;
    using ReconstructionAdvectionOperatorType =
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
    using FvOperatorType = EntropyBasedMomentFvOperator<
        std::conditional_t<TestCaseType::reconstruction, ReconstructionAdvectionOperatorType, AdvectionOperatorType>,
        EntropySolverType>;
    using OperatorTimeStepperType =
        ExplicitRungeKuttaTimeStepper<FvOperatorType,
                                      DiscreteFunctionType,
                                      TimeStepperMethods::explicit_rungekutta_second_order_ssp>;
    using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, MomentBasis>;
    using TimeStepperType = StrangSplittingTimeStepper<RhsTimeStepperType, OperatorTimeStepperType>;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;

    // *********************** create operators and timesteppers ************************************
    NumericalKineticFlux<GV, MomentBasis> numerical_flux(*analytical_flux, *basis_functions);
    AdvectionOperatorType advection_operator(grid_view, numerical_flux, advection_source_space, fv_space);

    // boundary treatment
    using BoundaryOperator =
        LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;
    using StateType = typename BoundaryOperator::StateType;
    LambdaType boundary_lambda =
        [&boundary_values](const I& intersection,
                           const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
                           const AnalyticalFluxType& /*flux*/,
                           const StateType& /*u*/,
                           const XT::Common::Parameter& /*param*/) {
          return boundary_values->evaluate(intersection.geometry().global(xx_in_reference_intersection_coordinates));
        };
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> filter;
    advection_operator.append(boundary_lambda, {}, filter);

    constexpr double epsilon = 1e-11;
    auto slope = TestCaseType::RealizabilityLimiterChooserType::template make_slope<EigenvectorWrapperType>(
        *dynamic_cast<EntropyFluxType*>(analytical_flux.get()), *basis_functions, epsilon);
    ReconstructionOperatorType reconstruction_operator(*analytical_flux, *boundary_values, fv_space, *slope, false);
    ReconstructionAdvectionOperatorType reconstruction_advection_operator(advection_operator, reconstruction_operator);

    if (XT::Common::is_zero(t_end))
      t_end = problem.t_end();

    if (!filename.empty())
      filename += "_";
    filename += ProblemType::static_id();
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad_" + XT::Common::to_string(quad_order);
    filename += MomentBasis::entropy == EntropyType::MaxwellBoltzmann ? "_MaxwellBoltzmann_" : "_BoseEinstein_";
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->mn_name();

    EntropySolverType entropy_solver(*(dynamic_cast<EntropyFluxType*>(analytical_flux.get())),
                                     fv_space,
                                     problem.psi_vac() * basis_functions->unit_ball_volume() / 10,
                                     filename);
    FvOperatorType fv_operator(
        FvOperatorChooser<TestCaseType::reconstruction>::choose(advection_operator, reconstruction_advection_operator),
        entropy_solver);

    // ******************************** do the time steps ***********************************************************
    const auto sigma_a = problem.sigma_a();
    const auto sigma_s = problem.sigma_s();
    const auto Q = problem.Q();
    OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    RhsTimeStepperType timestepper_rhs(*basis_functions, u, *sigma_a, *sigma_s, *Q);
    TimeStepperType timestepper(timestepper_rhs, timestepper_op);

    auto begin_time = std::chrono::steady_clock::now();
    timestepper.solve(t_end,
                      dt,
                      num_save_steps,
                      num_output_steps,
                      false,
                      true,
                      true,
                      false,
                      filename,
                      *basis_functions->visualizer(),
                      basis_functions->stringifier());
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (grid_view.comm().rank() == 0)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;

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
                     1,
                     0,
                     TestCaseType::quad_order,
                     TestCaseType::quad_refinements,
                     "",
                     2,
                     TestCaseType::t_end,
                     "test",
                     Dune::GDT::is_full_moment_basis<typename TestCaseType::MomentBasis>::value)
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

#endif // DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
