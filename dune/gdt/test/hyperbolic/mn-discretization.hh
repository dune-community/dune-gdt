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
  Dune::FieldVector<double, 3>
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
    using IntersectionType = typename GridLayerType::Intersection;
    using EquationType = Hyperbolic::Problems::KineticEquation<ProblemType>;
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
    if ((num_quad_refinements == size_t(-1) || quad_order == size_t(-1)) && (num_quad_refinements != quad_order))
      std::cerr << "You specified either num_quad_refinements or quad_order, please also specify the other one!"
                << std::endl;
    std::shared_ptr<const BasisfunctionType> basis_functions =
        (num_quad_refinements == size_t(-1) || quad_order == size_t(-1))
            ? std::make_shared<const BasisfunctionType>()
            : std::make_shared<const BasisfunctionType>(quad_order, num_quad_refinements);
    const std::unique_ptr<ProblemType> problem_imp =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_layer, grid_config);
    const EquationType problem(*problem_imp);
    const InitialValueType& initial_values = problem.initial_values();
    using BoundaryValueType =
        MomentModelBoundaryValue<GridLayerType, BasisfunctionType, typename EquationType::BoundaryValueType>;
    const auto& dirichlet_boundary_values = problem.boundary_values();
    const auto boundary_info = XT::Grid::AllDirichletBoundaryInfo<IntersectionType>();
    const BoundaryValueType boundary_values(boundary_info, *basis_functions, dirichlet_boundary_values);
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
    using RhsOperatorType = AdvectionRhsOperator<RhsType>;
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
    using RhsOperatorTimeStepperType =
        typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, TestCaseType::rhs_time_stepper_method>::
            TimeStepperType;
    using TimeStepperType =
        typename Dune::GDT::TimeStepperSplittingFactory<RhsOperatorTimeStepperType,
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
    RhsOperatorType rhs_operator(rhs);

    constexpr double epsilon = 1e-12;
    auto slope = TestCaseType::RealizabilityLimiterChooserType::template make_slope<
        typename ReconstructionOperatorType::MatrixType>(*basis_functions, epsilon);
    ReconstructionOperatorType reconstruction_operator(
        analytical_flux, boundary_values, *slope, Dune::GDT::default_1d_quadrature<double>(1));
    RegularizationOperatorType regularization_operator(analytical_flux);

    RealizabilityLimiterType realizability_limiter(analytical_flux, *basis_functions, epsilon);
    ReconstructionFvOperatorType reconstruction_fv_operator(
        advection_operator, reconstruction_operator, regularization_operator, realizability_limiter);
    NoReconstructionFvOperatorType no_reconstruction_fv_operator(advection_operator, regularization_operator);
    FvOperatorType& fv_operator = FvOperatorChooser<TestCaseType::reconstruction>::choose(no_reconstruction_fv_operator,
                                                                                          reconstruction_fv_operator);

    // ******************************** do the time steps ***********************************************************
    OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    RhsOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
    TimeStepperType timestepper(timestepper_rhs, timestepper_op);
    filename += "_" + ProblemType::static_id();
    filename += "_grid_" + grid_size;
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad" + XT::Common::to_string(num_quad_refinements) + "x" + XT::Common::to_string(quad_order);
    filename += std::is_same<FvOperatorType, ReconstructionFvOperatorType>::value ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->short_id();
    filename += "_m" + Dune::XT::Common::to_string(dimRange);

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

    FieldVector<double, 3> ret(0);
    double& l1norm = ret[0];
    double& l2norm = ret[1];
    double& linfnorm = ret[2];
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
    auto norms = HyperbolicMnDiscretization<TestCaseType>::run();
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