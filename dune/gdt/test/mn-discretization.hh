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
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/entropybasedmoments_fv.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/momentmodels/entropysolver.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/timestepper/fractional-step.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

#include "pn-discretization.hh"

template <class TestCaseType>
struct HyperbolicMnDiscretization
{
  // returns: (l1norm, l2norm, linfnorm, MPI rank)
  static std::pair<Dune::FieldVector<double, 3>, int>
  run(size_t num_save_steps = 1,
      size_t num_output_steps = 0,
      size_t num_quad_refinements = TestCaseType::RealizabilityLimiterChooserType::num_quad_refinements,
      size_t quad_order = TestCaseType::RealizabilityLimiterChooserType::quad_order,
      std::string grid_size = "",
      std::string overlap_size = "2",
      double t_end = TestCaseType::t_end,
      std::string filename = "")
  {
    using namespace Dune;
    using namespace Dune::GDT;

    //******************* get typedefs and constants from ProblemType **********************//
    using BasisfunctionType = typename TestCaseType::BasisfunctionType;
    using DiscreteFunctionType = typename TestCaseType::DiscreteFunctionType;
    using GridType = typename TestCaseType::GridType;
    using SpaceType = typename TestCaseType::SpaceType;
    using AdvectionSourceSpaceType = typename TestCaseType::AdvectionSourceSpaceType;
    using GridViewType = typename TestCaseType::GridViewType;
    using I = typename GridViewType::Intersection;
    using ProblemType = typename TestCaseType::ProblemType;
    using RangeFieldType = typename BasisfunctionType::RangeFieldType;
    using BoundaryValueType = typename ProblemType::BoundaryValueType;
    static constexpr size_t dimDomain = BasisfunctionType::dimDomain;
    static constexpr size_t dimRange = BasisfunctionType::dimRange;
    using MatrixType = typename XT::LA::Container<RangeFieldType>::MatrixType;
    using VectorType = typename XT::LA::Container<RangeFieldType>::VectorType;

    //******************* create grid and FV space ***************************************
    auto grid_config = ProblemType::default_grid_cfg();
    if (!grid_size.empty())
      grid_config["num_elements"] = grid_size;
    grid_config["overlap_size"] = overlap_size;
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GridViewType grid_view(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_view);
    const AdvectionSourceSpaceType advection_source_space(grid_view, 1);

    //******************* create EquationType object ***************************************
    std::shared_ptr<const BasisfunctionType> basis_functions =
        std::make_shared<const BasisfunctionType>(quad_order, num_quad_refinements);
    const std::unique_ptr<ProblemType> problem_ptr =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_config);
    const auto& problem = *problem_ptr;
    const auto initial_values = problem.initial_values();
    const auto boundary_values = problem.boundary_values();
    using AnalyticalFluxType = typename ProblemType::FluxType;
    using EntropyFluxType = typename ProblemType::ActualFluxType;
    const auto analytical_flux = problem.flux();
    const RangeFieldType CFL = problem.CFL();

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "solution");
    // project initial values
    default_interpolation(*initial_values, u, grid_view);

    // ******************** choose flux and rhs operator and timestepper ******************************************

    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GridViewType, dimRange>;
    using EigenvectorWrapperType = typename JacobianChooser<BasisfunctionType, AnalyticalFluxType>::type;
    using EntropySolverType = EntropySolver<BasisfunctionType, SpaceType>;
    using ReconstructionOperatorType = LinearReconstructionOperator<AnalyticalFluxType,
                                                                    BoundaryValueType,
                                                                    GridViewType,
                                                                    MatrixType,
                                                                    EigenvectorWrapperType>;
    using ReconstructionFvOperatorType = EntropyBasedMomentFvOperator<AdvectionOperatorType,
                                                                      ReconstructionOperatorType,
                                                                      RealizabilityLimiterType,
                                                                      RegularizationOperatorType>;
    using NoReconstructionFvOperatorType = EntropyBasedMomentFvOperator<AdvectionOperatorType, EntropySolverType>;
    using FvOperatorType =
        std::conditional_t<TestCaseType::reconstruction, ReconstructionFvOperatorType, NoReconstructionFvOperatorType>;
    using OperatorTimeStepperType =
        ExplicitRungeKuttaTimeStepper<FvOperatorType,
                                      DiscreteFunctionType,
                                      TimeStepperMethods::explicit_rungekutta_second_order_ssp>;
    using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, BasisfunctionType>;
    using TimeStepperType = StrangSplittingTimeStepper<RhsTimeStepperType, OperatorTimeStepperType>;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GridViewType> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;

    // *********************** create operators and timesteppers ************************************
    NumericalKineticFlux<BasisfunctionType> numerical_flux(*analytical_flux, *basis_functions);
    AdvectionOperatorType advection_operator(grid_view, numerical_flux, fv_space, fv_space);
    // boundary treatment
    using BoundaryOperator =
        LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GridViewType, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;
    using StateDofsType = typename BoundaryOperator::StateDofsType;
    LambdaType boundary_lambda =
        [&boundary_values](const I& intersection,
                           const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
                           const AnalyticalFluxType& /*flux*/,
                           const StateDofsType& /*u*/,
                           const XT::Common::Parameter& /*param*/) {
          return boundary_values->evaluate(intersection.geometry().global(xx_in_reference_intersection_coordinates));
        };
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridViewType> filter;
    advection_operator.append(boundary_lambda, {}, filter);

    // constexpr double epsilon = 1e-11;
    // auto slope = TestCaseType::RealizabilityLimiterChooserType::template make_slope<
    //     typename ReconstructionOperatorType::MatrixType>(*basis_functions, epsilon);
    // ReconstructionOperatorType reconstruction_operator(
    //     *analytical_flux, boundary_values, *slope, Dune::GDT::default_1d_quadrature<double>(1));

    filename += "_" + ProblemType::static_id();
    filename += "_grid_" + grid_size;
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad" + XT::Common::to_string(num_quad_refinements) + "x" + XT::Common::to_string(quad_order);
    //    filename += std::is_same<FvOperatorType, ReconstructionFvOperatorType>::value ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->short_id();
    filename += "_m" + Dune::XT::Common::to_string(dimRange);

    std::vector<typename EntropySolverType::LocalVectorType> alpha_vector(grid_view.size(0));
    EntropySolverType entropy_solver(*(dynamic_cast<const EntropyFluxType*>(analytical_flux.get())),
                                     fv_space,
                                     alpha_vector,
                                     problem.psi_vac() * basis_functions->unit_ball_volume() / 10,
                                     filename);

    // RealizabilityLimiterType realizability_limiter(*analytical_flux, *basis_functions, epsilon);
    // ReconstructionFvOperatorType reconstruction_fv_operator(
    //     advection_operator, reconstruction_operator, regularization_operator, realizability_limiter);
    NoReconstructionFvOperatorType no_reconstruction_fv_operator(advection_operator, entropy_solver);
    FvOperatorType& fv_operator = no_reconstruction_fv_operator;
    //    FvOperatorType& fv_operator =
    //    FvOperatorChooser<TestCaseType::reconstruction>::choose(no_reconstruction_fv_operator,
    //                                                                                          reconstruction_fv_operator);

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
  void run()
  {
    auto norms = HyperbolicMnDiscretization<TestCaseType>::run().first;
    const double l1norm = norms[0];
    const double l2norm = norms[1];
    const double linfnorm = norms[2];
    using ResultsType = typename TestCaseType::ExpectedResultsType;
    EXPECT_NEAR(ResultsType::l1norm, l1norm, ResultsType::l1norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::l2norm, l2norm, ResultsType::l2norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::linfnorm, linfnorm, ResultsType::linfnorm * ResultsType::tol);
  }
};

#endif // DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
