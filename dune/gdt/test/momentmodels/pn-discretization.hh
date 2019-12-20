// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH

#include <dune/common/exceptions.hh>

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/xt/functions/checkerboard.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/advection-with-reconstruction.hh>
#include <dune/gdt/operators/reconstruction/linear.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/numerical-fluxes/lax-friedrichs.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

void parse_momentmodel_arguments(int argc,
                                 char** argv,
                                 size_t& num_threads,
                                 size_t& threading_partition_factor,
                                 size_t& num_save_steps,
                                 size_t& num_output_steps,
                                 size_t& quad_order,
                                 size_t& quad_refinements,
                                 std::string& grid_size,
                                 size_t& overlap_size,
                                 double& t_end,
                                 std::string& filename)
{
  using namespace Dune;
  using namespace Dune::GDT;
  MPIHelper::instance(argc, argv);

  // default values
  num_threads = 1;
  threading_partition_factor = 1;
  num_save_steps = 10;
  num_output_steps = num_save_steps;
  quad_order = -1;
  quad_refinements = -1;
  grid_size = "";
  overlap_size = 2;
  t_end = 0.;
  filename = "";

  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--num_threads") {
      if (i + 1 < argc)
        num_threads = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--num_threads option requires one argument.");
    } else if (std::string(argv[i]) == "--threading_partition_factor") {
      if (i + 1 < argc)
        threading_partition_factor = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--threading_partition_factor option requires one argument.");
    } else if (std::string(argv[i]) == "--filename") {
      if (i + 1 < argc)
        filename = std::string(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--filename option requires one argument.");
    } else if (std::string(argv[i]) == "--quad_order") {
      if (i + 1 < argc)
        quad_order = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--quad_order option requires one argument.");
    } else if (std::string(argv[i]) == "--quad_refinements") {
      if (i + 1 < argc)
        quad_refinements = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--quad_refinements option requires one argument.");
    } else if (std::string(argv[i]) == "--num_save_steps") {
      if (i + 1 < argc)
        num_save_steps = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--num_save_steps option requires one argument.");
    } else if (std::string(argv[i]) == "--num_output_steps") {
      if (i + 1 < argc)
        num_output_steps = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--num_output_steps option requires one argument.");
    } else if (std::string(argv[i]) == "--grid_size") {
      if (i + 1 < argc)
        grid_size = argv[++i];
      else
        DUNE_THROW(Dune::IOError, "--grid_size option requires one argument.");
    } else if (std::string(argv[i]) == "--overlap_size") {
      if (i + 1 < argc)
        overlap_size = XT::Common::from_string<size_t>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--overlap_size option requires one argument.");
    } else if (std::string(argv[i]) == "--t_end") {
      if (i + 1 < argc)
        t_end = XT::Common::from_string<double>(argv[++i]);
      else
        DUNE_THROW(Dune::IOError, "--t_end option requires one argument.");
    } else {
      DUNE_THROW(Dune::IOError, "Unknown option " + std::string(argv[i]));
    }
  }

  DXTC_CONFIG.set("threading.partition_factor", threading_partition_factor, true);
  XT::Common::threadManager().set_max_threads(num_threads);
}


template <class MomentBasis, class AnalyticalFluxType>
struct EigenvectorWrapperChooser
{
  using type = Dune::GDT::internal::EigenvectorWrapper<AnalyticalFluxType>;
};

template <class AnalyticalFluxType,
          class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange_or_refinements,
          Dune::GDT::EntropyType entropy>
struct EigenvectorWrapperChooser<Dune::GDT::PartialMomentBasis<DomainFieldType,
                                                               dimDomain,
                                                               RangeFieldType,
                                                               dimRange_or_refinements,
                                                               1,
                                                               dimDomain,
                                                               1,
                                                               entropy>,
                                 AnalyticalFluxType>
{
  using type = Dune::GDT::internal::BlockedEigenvectorWrapper<AnalyticalFluxType>;
};

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

#ifndef USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
#  define USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR 0
#endif

template <class TestCaseType>
struct HyperbolicPnDiscretization
{

  Dune::FieldVector<double, 3> run(size_t num_save_steps = 1,
                                   size_t num_output_steps = 0,
                                   size_t quad_order = size_t(-1),
                                   size_t quad_refinements = size_t(-1),
                                   std::string grid_size = "",
                                   size_t overlap_size = 2,
                                   double t_end = 0.,
                                   std::string filename = "")
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
    using E = XT::Grid::extract_entity_t<GV>;
    using ProblemType = typename TestCaseType::ProblemType;
    using RangeFieldType = typename ProblemType::RangeFieldType;
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
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_config);
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
    typename DiscreteFunctionType::VectorType vec(vec_size, 0., num_mutexes);
    DiscreteFunctionType u(fv_space, vec, "solution");
    // project initial values
    default_interpolation(*initial_values, u, grid_view);

    // ************************* create analytical flux object ***************************************
    using AnalyticalFluxType = typename ProblemType::FluxType;
    const auto analytical_flux = problem.flux();

    // ******************** choose flux and rhs operator and timestepper ******************************************
    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using EigenvectorWrapperType = typename EigenvectorWrapperChooser<MomentBasis, AnalyticalFluxType>::type;
    using ReconstructionOperatorType =
#if USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
        LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, GV, MatrixType, EigenvectorWrapperType>;
#else
        PointwiseLinearReconstructionOperator<AnalyticalFluxType,
                                              BoundaryValueType,
                                              GV,
                                              VectorType,
                                              EigenvectorWrapperType>;
#endif

    using ReconstructionFvOperatorType =
#if USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR
        AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
#else
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
#endif
    using FvOperatorType =
        std::conditional_t<TestCaseType::reconstruction, ReconstructionFvOperatorType, AdvectionOperatorType>;
    using OperatorTimeStepperType =
        ExplicitRungeKuttaTimeStepper<FvOperatorType,
                                      DiscreteFunctionType,
                                      TimeStepperMethods::explicit_rungekutta_second_order_ssp>;
    using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, MomentBasis>;
    using TimeStepperType = StrangSplittingTimeStepper<RhsTimeStepperType, OperatorTimeStepperType>;

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
    NumericalKineticFlux<GV, MomentBasis> numerical_flux(*analytical_flux, *basis_functions);
    //    NumericalLaxFriedrichsFlux<I, dimDomain, dimRange, RangeFieldType> numerical_flux(*analytical_flux, 1.);
    AdvectionOperatorType advection_operator(grid_view, numerical_flux, advection_source_space, fv_space);
    // boundary treatment
    using BoundaryOperator =
        LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;
    using DynamicStateType = typename BoundaryOperator::DynamicStateType;
    LambdaType boundary_lambda =
        [&boundary_values](const I& intersection,
                           const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
                           const AnalyticalFluxType& /*flux*/,
                           const DynamicStateType& /*u*/,
                           DynamicStateType& v,
                           const XT::Common::Parameter& param) {
          boundary_values->evaluate(intersection.geometry().global(xx_in_reference_intersection_coordinates), v, param);
        };
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> filter;
    advection_operator.append(boundary_lambda, {}, filter);

    MinmodSlope<E, EigenvectorWrapperType> slope;
    ReconstructionOperatorType reconstruction_operator(*analytical_flux, *boundary_values, fv_space, slope, true);

    ReconstructionFvOperatorType reconstruction_fv_operator(advection_operator, reconstruction_operator);
    FvOperatorType& fv_operator =
        FvOperatorChooser<TestCaseType::reconstruction>::choose(advection_operator, reconstruction_fv_operator);

    if (!filename.empty())
      filename += "_";
    filename += ProblemType::static_id();
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->pn_name();

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

    FieldVector<double, 3> ret(0);
    double& l1norm = ret[0];
    double& l2norm = ret[1];
    double& linfnorm = ret[2];
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
struct HyperbolicPnTest
  : public HyperbolicPnDiscretization<TestCaseType>
  , public ::testing::Test
{
  void run()
  {
    auto norms =
        HyperbolicPnDiscretization<TestCaseType>::run(1, 0, size_t(-1), size_t(-1), "", 2, TestCaseType::t_end, "test");
    const double l1norm = norms[0];
    const double l2norm = norms[1];
    const double linfnorm = norms[2];
    using ResultsType = typename TestCaseType::ExpectedResultsType;
    EXPECT_NEAR(ResultsType::l1norm, l1norm, ResultsType::l1norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::l2norm, l2norm, ResultsType::l2norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::linfnorm, linfnorm, ResultsType::linfnorm * ResultsType::tol);
  }
};

#endif // DUNE_GDT_TEST_HYPERBOLIC_PN_DISCRETIZATION_HH
