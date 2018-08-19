// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH

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

template <bool reconstruction>
struct FvOperatorChooser
{
  template <class AdvectionOperatorType, class ReconstructionOperatorType>
  static ReconstructionOperatorType& choose(AdvectionOperatorType& /*advection_operator*/,
                                            ReconstructionOperatorType& reconstruction_operator)
  {
    return reconstruction_operator;
  }
};

template <>
struct FvOperatorChooser<false>
{
  template <class AdvectionOperatorType, class ReconstructionOperatorType>
  static AdvectionOperatorType& choose(AdvectionOperatorType& advection_operator,
                                       ReconstructionOperatorType& /*reconstruction_operator*/)
  {
    return advection_operator;
  }
};

template <class TestCaseType>
struct HyperbolicPnTest : public ::testing::Test
{
  void run()
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
    using DomainFieldType = typename EquationType::DomainFieldType;
    using RangeFieldType = typename EquationType::RangeFieldType;
    using RhsType = typename EquationType::RhsType;
    using InitialValueType = typename EquationType::InitialValueType;
    static constexpr size_t dimDomain = BasisfunctionType::dimDomain;
    static constexpr size_t dimRange = BasisfunctionType::dimRange;

    //******************* create grid and FV space ***************************************
    auto grid_config = EquationType::default_grid_cfg();
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GridLayerType grid_layer(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_layer);

    //******************* create EquationType object ***************************************
    const auto quadrature = ProblemType::default_quadrature(grid_config);
    std::shared_ptr<const BasisfunctionType> basis_functions = std::make_shared<const BasisfunctionType>();
    const std::unique_ptr<ProblemType> problem_imp =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_layer, quadrature, grid_config);
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
    project_l2(initial_values, u);

    // ************************* create analytical flux object ***************************************
    using AnalyticalFluxType = typename EquationType::FluxType;
    const AnalyticalFluxType& analytical_flux = problem.flux();

    // ******************** choose flux and rhs operator and timestepper ******************************************
    using RhsOperatorType = AdvectionRhsOperator<RhsType>;

    using AdvectionOperatorType =
        AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType, BasisfunctionType, GridLayerType>;

    using JacobianWrapperType = std::
        conditional_t<std::is_base_of<
                          typename Hyperbolic::Problems::
                              PiecewiseMonomials<DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, dimDomain>,
                          BasisfunctionType>::value,
                      internal::BlockedJacobianWrapper<AnalyticalFluxType>,
                      internal::JacobianWrapper<AnalyticalFluxType>>;
    using ReconstructionOperatorType =
        LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, JacobianWrapperType>;
    using ReconstructionFvOperatorType =
        AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
    using FvOperatorType =
        std::conditional_t<TestCaseType::reconstruction, ReconstructionFvOperatorType, AdvectionOperatorType>;

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

    // *************** choose t_end and initial dt **************************************
    // calculate dx and choose initial dt
    Dune::XT::Grid::Dimensions<GridLayerType> dimensions(grid_layer);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;
    const RangeFieldType t_end = TestCaseType::t_end;
    //    const RangeFieldType t_end = 2.5;

    // *********************** create operators and timesteppers ************************************
    AdvectionOperatorType advection_operator(analytical_flux, boundary_values, *basis_functions);
    RhsOperatorType rhs_operator(rhs);

    MinmodSlope<typename ReconstructionOperatorType::VectorType, typename ReconstructionOperatorType::MatrixType> slope;
    ReconstructionOperatorType reconstruction_operator(
        analytical_flux, boundary_values, slope, Dune::GDT::default_1d_quadrature<double>(1));

    ReconstructionFvOperatorType reconstruction_fv_operator(advection_operator, reconstruction_operator);
    FvOperatorType& fv_operator =
        FvOperatorChooser<TestCaseType::reconstruction>::choose(advection_operator, reconstruction_fv_operator);

    // ******************************** do the time steps ***********************************************************
    OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    RhsOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
    TimeStepperType timestepper(timestepper_rhs, timestepper_op);

    auto visualizer = basis_functions->template visualizer<DiscreteFunctionType>();
    timestepper.solve(t_end,
                      dt,
                      1,
                      0,
                      false,
                      false,
                      true,
                      false,
                      problem_imp->static_id() + "_ord" + (TestCaseType::reconstruction ? "2_" : "1_")
                          + basis_functions->short_id(),
                      visualizer,
                      basis_functions->stringifier());

    RangeFieldType l1norm = 0;
    RangeFieldType l2norm = 0;
    RangeFieldType linfnorm = 0;
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
    EXPECT_DOUBLE_EQ(TestCaseType::ExpectedResultsType::l1norm, l1norm);
    EXPECT_DOUBLE_EQ(TestCaseType::ExpectedResultsType::l2norm, l2norm);
    EXPECT_DOUBLE_EQ(TestCaseType::ExpectedResultsType::linfnorm, linfnorm);
  }
};

#endif // DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH
