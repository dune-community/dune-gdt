// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2019)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH

#include <dune/common/exceptions.hh>

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/la/container.hh>

#include <dune/xt/functions/checkerboard.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/advection-fv-with-reconstruction.hh>
#include <dune/gdt/operators/reconstruction/linear.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/numerical-fluxes/lax-friedrichs.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>
#include "pn-discretization.hh"

#ifndef USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
#  define USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR 0
#endif

template <class Grid, class Problem, bool reconstruct>
struct HyperbolicFvDiscretization
{

  void run(size_t num_save_steps = size_t(-1),
           size_t num_output_steps = size_t(-1),
           std::string grid_size = "",
           size_t overlap_size = 2,
           double t_end = 1.,
           std::string filename = "")
  {
    using namespace Dune;
    using namespace Dune::GDT;


    //******************* get typedefs and constants **********************//
    using G = Grid;
    using ProblemType = Problem;
    static constexpr size_t dimDomain = Grid::dimension;
    static constexpr size_t dimRange = ProblemType::dimRange;
    // using DomainFieldType = typename G::ctype;
    using RangeFieldType = typename ProblemType::RangeFieldType;
    using GV = XT::Grid::PeriodicGridView<typename G::LeafGridView>;
    // using GV = typename G::LeafGridView;
    using I = XT::Grid::extract_intersection_t<GV>;
    using E = XT::Grid::extract_entity_t<GV>;
    using SpaceType = FiniteVolumeSpace<GV, dimRange, 1, RangeFieldType>;
    using AdvectionSourceSpaceType =
        std::conditional_t<reconstruct, DiscontinuousLagrangeSpace<GV, dimRange, RangeFieldType>, SpaceType>;
    using BoundaryValueType = typename ProblemType::BoundaryValueType;
    using MatrixType = typename XT::LA::Container<RangeFieldType>::MatrixType;
    using VectorType = typename XT::LA::Container<RangeFieldType>::VectorType;
    using DiscreteFunctionType = DiscreteFunction<VectorType, GV, dimRange, 1, RangeFieldType>;

    //******************* create grid and FV space ***************************************
    auto grid_config = ProblemType::default_grid_cfg();
    if (!grid_size.empty())
      grid_config["num_elements"] = grid_size;
    grid_config["overlap_size"] = XT::Common::to_string(overlap_size);
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<G>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GV grid_view(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_view);
    const AdvectionSourceSpaceType advection_source_space(grid_view, 1);

    //******************* create EquationType object ***************************************
    const std::unique_ptr<ProblemType> problem_ptr = XT::Common::make_unique<ProblemType>(grid_config);
    const auto& problem = *problem_ptr;
    const auto initial_values = problem.initial_values();
    const auto boundary_values = problem.boundary_values();
    const RangeFieldType CFL = problem.CFL();

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    const size_t vec_size = fv_space.mapper().size();
    // we do very few whole-container operations with this vec, so using that many mutexes improves performance as it
    // avoids locking
    const size_t num_mutexes = XT::Common::threadManager().max_threads() * 100;
    VectorType vec(vec_size, 0., num_mutexes);
    DiscreteFunctionType u(fv_space, vec, "solution");
    // project initial values
    default_interpolation(*initial_values, u, grid_view);

    // ************************* create analytical flux object ***************************************
    using AnalyticalFluxType = typename ProblemType::FluxType;
    const auto analytical_flux = problem.flux();

    // ******************** choose flux and rhs operator and timestepper ******************************************
    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using ReconstructionOperatorType =
#if USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
        LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, GV, MatrixType>;
#else
        PointwiseLinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, GV, VectorType>;
#endif

    using ReconstructionFvOperatorType =
#if USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
        AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
#else
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
#endif
    using FvOperatorType = std::conditional_t<reconstruct, ReconstructionFvOperatorType, AdvectionOperatorType>;
    using OperatorTimeStepperType =
        ExplicitRungeKuttaTimeStepper<FvOperatorType,
                                      DiscreteFunctionType,
                                      TimeStepperMethods::explicit_rungekutta_second_order_ssp>;
    // using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, MomentBasis>;
    using TimeStepperType =
        OperatorTimeStepperType; // StrangSplittingTimeStepper<RhsTimeStepperType, OperatorTimeStepperType>;

    // *************** choose t_end and initial dt **************************************
    // calculate dx and choose initial dt
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;

    // *********************** create operators and timesteppers ************************************
    NumericalLaxFriedrichsFlux<I, dimDomain, dimRange, RangeFieldType> numerical_flux(*analytical_flux, 1.);
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

    MinmodSlope<E, Dune::GDT::internal::EigenvectorWrapper<AnalyticalFluxType>> slope;
    ReconstructionOperatorType reconstruction_operator(*analytical_flux, *boundary_values, fv_space, slope, true);

    ReconstructionFvOperatorType reconstruction_fv_operator(advection_operator, reconstruction_operator);
    FvOperatorType& fv_operator =
        FvOperatorChooser<reconstruct>::choose(advection_operator, reconstruction_fv_operator);

    if (!filename.empty())
      filename += "_";
    filename += ProblemType::static_id();
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += reconstruct ? "_ord2" : "_ord1";

    // ******************************** do the time steps ***********************************************************
    OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
    // RhsTimeStepperType timestepper_rhs(*basis_functions, u, *sigma_a, *sigma_s, *Q);
    // TimeStepperType timestepper(timestepper_rhs, timestepper_op);
    TimeStepperType& timestepper = timestepper_op;

    auto begin_time = std::chrono::steady_clock::now();
    timestepper.solve(t_end, dt, num_save_steps, num_output_steps, false, true, true, false, filename);
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (grid_view.comm().rank() == 0)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;
  }
};


#endif // DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
