// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_MN_DISCRETIZATION_HH

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/operators/fv/quadrature.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/kineticequation.hh>

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
    using GridLayerType = typename TestCaseType::GridLayerType;
    using ProblemType = typename TestCaseType::ProblemType;
    using EquationType = Hyperbolic::Problems::KineticEquation<ProblemType>;
    using DomainFieldType = typename EquationType::DomainFieldType;
    using RangeFieldType = typename EquationType::RangeFieldType;
    using RhsType = typename EquationType::RhsType;
    using InitialValueType = typename EquationType::InitialValueType;
    static constexpr size_t dimDomain = BasisfunctionType::dimDomain;
    static constexpr size_t dimRange = BasisfunctionType::dimRange;

    //******************* create grid and FV space ***************************************
    auto grid_config = EquationType::default_grid_cfg();
    if (!grid_size.empty())
      grid_config["num_elements"] = grid_size;
    grid_config["overlap_size"] = overlap_size;
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GridLayerType grid_layer(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_layer);

    //******************* create EquationType object ***************************************
    std::shared_ptr<const BasisfunctionType> basis_functions =
        std::make_shared<const BasisfunctionType>(quad_order, num_quad_refinements);
    const std::unique_ptr<ProblemType> problem_imp =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_layer, grid_config);
    const EquationType problem(*problem_imp);
    const InitialValueType& initial_values = problem.initial_values();
    using BoundaryValueType = typename ProblemType::BoundaryValueType;
    const BoundaryValueType& boundary_values = problem.boundary_values();
    const RhsType& rhs = problem.rhs();
    const RangeFieldType CFL = problem.CFL();

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "solution");
    // project initial values
    project_l2(initial_values, u, 0, {}, true);

    // ************************* create analytical flux object ***************************************
    using AnalyticalFluxType = typename EquationType::FluxType;
    const AnalyticalFluxType& analytical_flux = problem.flux();

    //***************** choose RealizabilityLimiter ********************************************
    using LocalRealizabilityLimiterType =
        typename TestCaseType::RealizabilityLimiterChooserType::LocalRealizabilityLimiterType;
    using RealizabilityLimiterType = RealizabilityLimiter<LocalRealizabilityLimiterType>;
    using JacobianWrapperType = typename JacobianChooser<BasisfunctionType, AnalyticalFluxType>::type;

    // ******************** choose flux and rhs operator and timestepper ******************************************
    using AdvectionOperatorType =
        AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType, BasisfunctionType, GridLayerType>;
    using ReconstructionOperatorType =
        LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, JacobianWrapperType>;
    using RegularizationOperatorType = MomentRegularizer<AnalyticalFluxType, BasisfunctionType>;
    using ReconstructionFvOperatorType = EntropyBasedMomentFvOperator<AdvectionOperatorType,
                                                                      ReconstructionOperatorType,
                                                                      RealizabilityLimiterType,
                                                                      RegularizationOperatorType>;
    using NoReconstructionFvOperatorType =
        EntropyBasedMomentFvOperatorNoReconstruction<AdvectionOperatorType, RegularizationOperatorType>;
    using FvOperatorType =
        std::conditional_t<TestCaseType::reconstruction, ReconstructionFvOperatorType, NoReconstructionFvOperatorType>;
    using OperatorTimeStepperType = typename TimeStepperFactory<FvOperatorType,
                                                                DiscreteFunctionType,
                                                                TestCaseType::time_stepper_method>::TimeStepperType;
    using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, BasisfunctionType>;
    using TimeStepperType =
        typename Dune::GDT::TimeStepperSplittingFactory<RhsTimeStepperType,
                                                        OperatorTimeStepperType,
                                                        TestCaseType::time_stepper_splitting_method>::TimeStepperType;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GridLayerType> dimensions(grid_layer);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;

    // *********************** create operators and timesteppers ************************************
    AdvectionOperatorType advection_operator(analytical_flux, boundary_values, *basis_functions);

    constexpr double epsilon = 1e-11;
    auto slope = TestCaseType::RealizabilityLimiterChooserType::template make_slope<
        typename ReconstructionOperatorType::MatrixType>(*basis_functions, epsilon);
    ReconstructionOperatorType reconstruction_operator(
        analytical_flux, boundary_values, *slope, Dune::GDT::default_1d_quadrature<double>(1));

    filename += "_" + ProblemType::static_id();
    filename += "_grid_" + grid_size;
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad" + XT::Common::to_string(num_quad_refinements) + "x" + XT::Common::to_string(quad_order);
    filename += std::is_same<FvOperatorType, ReconstructionFvOperatorType>::value ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->short_id();
    filename += "_m" + Dune::XT::Common::to_string(dimRange);

    RegularizationOperatorType regularization_operator(
        analytical_flux, problem_imp->psi_vac() * basis_functions->unit_ball_volume() / 10, filename);

    RealizabilityLimiterType realizability_limiter(analytical_flux, *basis_functions, epsilon);
    ReconstructionFvOperatorType reconstruction_fv_operator(
        advection_operator, reconstruction_operator, regularization_operator, realizability_limiter);
    NoReconstructionFvOperatorType no_reconstruction_fv_operator(advection_operator, regularization_operator);
    FvOperatorType& fv_operator = FvOperatorChooser<TestCaseType::reconstruction>::choose(no_reconstruction_fv_operator,
                                                                                          reconstruction_fv_operator);

    // ******************************** do the time steps ***********************************************************
    OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    RhsTimeStepperType timestepper_rhs(
        *basis_functions, u, problem_imp->get_sigma_a(), problem_imp->get_sigma_s(), problem_imp->get_Q());
    TimeStepperType timestepper(timestepper_rhs, timestepper_op);
    auto visualizer = basis_functions->template visualizer<DiscreteFunctionType>();
    timestepper.solve(t_end,
                      dt,
                      num_save_steps,
                      num_output_steps,
                      false,
                      true,
                      true,
                      false,
                      filename,
                      visualizer,
                      basis_functions->stringifier());

    auto ret = std::make_pair(FieldVector<double, 3>(0.), int(0));
    double& l1norm = ret.first[0];
    double& l2norm = ret.first[1];
    double& linfnorm = ret.first[2];
    ret.second = grid_layer.comm().rank();
    const auto& current_sol = timestepper.current_solution();
    for (const auto& entity : elements(grid_layer, Dune::Partitions::interior)) {
      const auto local_sol = current_sol.local_function(entity);
      const auto val = local_sol->evaluate(entity.geometry().local(entity.geometry().center()));
      RangeFieldType psi = basis_functions->density(val);
      l1norm += std::abs(psi) * entity.geometry().volume();
      l2norm += std::pow(psi, 2) * entity.geometry().volume();
      linfnorm = std::max(std::abs(psi), linfnorm);
    }
    l1norm = grid_layer.comm().sum(l1norm);
    l2norm = grid_layer.comm().sum(l2norm);
    linfnorm = grid_layer.comm().max(linfnorm);
    l2norm = std::sqrt(l2norm);
    return ret;
  }
};

template <class TestCaseType>
struct HyperbolicMnTest : public HyperbolicMnDiscretization<TestCaseType>, public ::testing::Test
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
