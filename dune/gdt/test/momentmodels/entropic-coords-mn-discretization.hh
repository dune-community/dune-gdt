// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_DISCRETIZATION_HH

#include <chrono>

#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/local/operators/generic.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/advection-fv-entropybased.hh>
#include <dune/gdt/operators/localizable-operator.hh>
#include <dune/gdt/operators/reconstruction/linear_kinetic.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>
#include <dune/gdt/test/momentmodels/hessianinverter.hh>
#include <dune/gdt/test/momentmodels/density_evaluations.hh>
#include <dune/gdt/tools/timestepper/adaptive-rungekutta-kinetic.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

#include "pn-discretization.hh"

template <class TestCaseType>
struct HyperbolicEntropicCoordsMnDiscretization
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
                                                          bool /*disable_thread_cache*/ = false)
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
    //    using E = XT::Grid::extract_entity_t<GV>;
    using I = XT::Grid::extract_intersection_t<GV>;
    using ProblemType = typename TestCaseType::ProblemType;
    using RangeFieldType = typename MomentBasis::RangeFieldType;
    static constexpr size_t dimDomain = MomentBasis::dimDomain;
    static constexpr size_t dimRange = MomentBasis::dimRange;
    static const auto la_backend = TestCaseType::la_backend;
    using MatrixType = typename XT::LA::Container<RangeFieldType, la_backend>::MatrixType;
    using VectorType = typename XT::LA::Container<RangeFieldType, la_backend>::VectorType;
    using GenericFunctionType = XT::Functions::GenericFunction<dimDomain, dimRange, 1, RangeFieldType>;
    using DomainType = FieldVector<RangeFieldType, dimDomain>;
    using RangeType = FieldVector<RangeFieldType, dimRange>;
    using DynamicRangeType = DynamicVector<RangeFieldType>;
    //        static const RangeFieldType scale_factor = 1e2;

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
    const auto initial_values_u = problem.initial_values();
    const auto boundary_values_u = problem.boundary_values();
    const auto boundary_distribution = problem.boundary_distribution();

    using AnalyticalFluxType = typename ProblemType::FluxType;
    constexpr SlopeLimiterType slope =
        TestCaseType::reconstruction ? SlopeLimiterType::minmod : SlopeLimiterType::no_slope;
    using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GV, MomentBasis, slope>;
    using OldEntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
    auto flux = problem.flux();
    auto* entropy_flux = dynamic_cast<OldEntropyFluxType*>(flux.get());
    auto analytical_flux = std::make_unique<EntropyFluxType>(*entropy_flux);
    const RangeFieldType CFL = problem.CFL();

    // calculate boundary values for alpha
    std::map<DomainType, RangeType, XT::Common::FieldVectorFloatLess> alpha_boundary_vals;
    for (const auto& element : Dune::elements(grid_view))
      for (const auto& intersection : Dune::intersections(grid_view, element))
        if (intersection.boundary()) {
          const auto x = intersection.geometry().center();
          const auto u = boundary_values_u->evaluate(x);
          alpha_boundary_vals.insert(std::make_pair(x, analytical_flux->get_alpha(u)));
        }
    GenericFunctionType boundary_values_alpha(
        1, [&](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) {
          ret = alpha_boundary_vals[x];
        });


    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "u_initial");
    DiscreteFunctionType alpha(fv_space, "alpha_initial");
    // project initial values
    default_interpolation(*initial_values_u, u, grid_view);

    // convert initial values to alpha
    // using EntropySolverType = EntropySolver<MomentBasis, SpaceType, MatrixType>;
    // EntropySolverType entropy_solver(*entropy_flux,
    //                                  fv_space,
    //                                  problem.psi_vac() * basis_functions->unit_ball_volume() / 10,
    //                                  filename);
    // entropy_solver.apply(u.dofs().vector(), alpha.dofs().vector(), {});

    const auto u_local_func = u.local_discrete_function();
    const auto alpha_local_func = alpha.local_discrete_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> u_local;
    for (auto&& element : Dune::elements(grid_view)) {
      u_local_func->bind(element);
      alpha_local_func->bind(element);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_local[ii] = u_local_func->dofs().get_entry(ii);
      const auto alpha_local = analytical_flux->get_alpha(u_local);
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_local_func->dofs().set_entry(ii, alpha_local[ii]);
    }

    // ******************** choose flux and rhs operator and timestepper ******************************************

    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using HessianInverterType = EntropicHessianInverter<MomentBasis, SpaceType, slope, MatrixType>;
#if 0
    using ReconstructionOperatorType = PointwiseLinearReconstructionNoCharOperator<
                                                                             GV,
                                                                             BoundaryValueType,
        EntropyFluxType,
                                                                             VectorType,
                                                                             RangeType>;
#else
    using ReconstructionOperatorType =
        PointwiseLinearKineticReconstructionOperator<GV, EntropyFluxType, VectorType, RangeType>;
#endif
    using ReconstructionAdvectionOperatorType =
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
    using NonEntropicAdvectionOperatorType = ReconstructionAdvectionOperatorType;
    //    using FvOperatorType = EntropicCoordinatesOperator<
    //        NonEntropicAdvectionOperatorType,
    //        HessianInverterType>;
    using NonEntropicRhsOperatorType = LocalizableOperator<MatrixType, GV, dimRange>;
    //    using RhsOperatorType = EntropicCoordinatesOperator<NonEntropicRhsOperatorType, HessianInverterType>;
    using DensityOperatorType = DensityEvaluator<MomentBasis, SpaceType, slope, MatrixType>;
    using CombinedOperatorType = EntropicCoordinatesCombinedOperator<DensityOperatorType,
                                                                     NonEntropicAdvectionOperatorType,
                                                                     NonEntropicRhsOperatorType,
                                                                     HessianInverterType>;

    //    using OperatorTimeStepperType =
    //        ExplicitRungeKuttaTimeStepper<FvOperatorType,
    //                                      DiscreteFunctionType,
    //                                      TimeStepperMethods::explicit_euler>;
    //    using RhsTimeStepperType =
    //        ExplicitRungeKuttaTimeStepper<RhsOperatorType,
    //                                      DiscreteFunctionType,
    //                                      TimeStepperMethods::explicit_euler>;

    //    using OperatorTimeStepperType =
    //        KineticAdaptiveRungeKuttaTimeStepper<FvOperatorType,
    //                                      DiscreteFunctionType>;
    //    using RhsTimeStepperType =
    //        KineticAdaptiveRungeKuttaTimeStepper<RhsOperatorType,
    //                                      DiscreteFunctionType>;

    //    using TimeStepperType = FractionalTimeStepper<OperatorTimeStepperType, RhsTimeStepperType>;

    using TimeStepperType = KineticAdaptiveRungeKuttaTimeStepper<CombinedOperatorType,
                                                                 DiscreteFunctionType,
                                                                 EntropyFluxType,
                                                                 TimeStepperMethods::bogacki_shampine>;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);
    RangeFieldType dt = CFL * dx;

    // *********************** create operators and timesteppers ************************************
    NumericalKineticFlux<GV, MomentBasis, EntropyFluxType> numerical_flux(*analytical_flux, *basis_functions);
    AdvectionOperatorType advection_operator(grid_view, numerical_flux, advection_source_space, fv_space);

    // boundary treatment
    using BoundaryOperator =
        LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;

    // calculate boundary kinetic fluxes
    // apply density_operator first to store boundary_evaluations
    const double min_acceptable_density = problem.psi_vac() / 10;
    DensityOperatorType density_operator(*analytical_flux, fv_space, boundary_distribution, min_acceptable_density);
    density_operator.apply(alpha.dofs().vector(), alpha.dofs().vector());

    // store boundary fluxes
    std::map<DomainType, DynamicRangeType, XT::Common::FieldVectorFloatLess> boundary_fluxes;
    for (const auto& element : Dune::elements(grid_view))
      for (const auto& intersection : Dune::intersections(grid_view, element))
        if (intersection.boundary()) {
          const auto x = intersection.geometry().center();
          const auto dd = intersection.indexInInside() / 2;
          const DynamicRangeType boundary_flux =
              problem.kinetic_boundary_flux(x, intersection.centerUnitOuterNormal()[dd], dd);
          boundary_fluxes.insert(std::make_pair(x, boundary_flux));
        }
    GenericFunctionType boundary_kinetic_fluxes(
        1, [&](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) { ret = boundary_fluxes[x]; });

    LambdaType boundary_lambda =
        [&](const I& intersection,
            const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
            const AnalyticalFluxType& /*flux*/,
            const DynamicRangeType& /*u*/,
            DynamicRangeType& v,
            const XT::Common::Parameter& param) {
          boundary_kinetic_fluxes.evaluate(
              intersection.geometry().global(xx_in_reference_intersection_coordinates), v, param);
        };
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> filter;
    advection_operator.append(boundary_lambda, {}, filter);

    ReconstructionOperatorType reconstruction_operator(fv_space, *analytical_flux);
    ReconstructionAdvectionOperatorType reconstruction_advection_operator(advection_operator, reconstruction_operator);

    if (XT::Common::is_zero(t_end))
      t_end = problem.t_end();

    if (!filename.empty())
      filename += "_";
    if (TestCaseType::reconstruction && slope == SlopeLimiterType::minmod)
      filename += "minmod_";
    else if (TestCaseType::reconstruction && slope == SlopeLimiterType::superbee)
      filename += "superbee_";
    filename += ProblemType::static_id();
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad_" + XT::Common::to_string(quad_order);
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->mn_name();

    HessianInverterType hessian_inverter(*analytical_flux, fv_space);
    auto& non_entropic_advection_operator = reconstruction_advection_operator;

    static const RangeType u_iso = basis_functions->u_iso();
    static const RangeType basis_integrated = basis_functions->integrated();
    const auto sigma_a = problem.sigma_a();
    const auto sigma_s = problem.sigma_s();
    const auto Q = problem.Q();
    auto rhs_func = [&](const auto& /*source*/,
                        const auto& /*local_source*/,
                        auto& local_range,
                        const Dune::XT::Common::Parameter& /*param*/) {
      const auto& element = local_range.element();
      const auto center = element.geometry().center();
      const auto u_elem = analytical_flux->get_u(fv_space.grid_view().indexSet().index(element));
      const auto sigma_a_value = sigma_a->evaluate(center)[0];
      const auto sigma_s_value = sigma_s->evaluate(center)[0];
      const auto sigma_t_value = sigma_a_value + sigma_s_value;
      const auto Q_value = Q->evaluate(center)[0];
      auto ret = u_elem;
      ret *= -sigma_t_value;
      ret.axpy(basis_functions->density(u_elem) * sigma_s_value, u_iso);
      ret.axpy(Q_value, basis_integrated);
      auto& range_dofs = local_range.dofs();
      for (size_t ii = 0; ii < dimRange; ++ii)
        range_dofs[ii] += ret[ii];
    };
    NonEntropicRhsOperatorType non_entropic_rhs_operator(grid_view, fv_space, fv_space);
    non_entropic_rhs_operator.append(GenericLocalElementOperator<VectorType, GV, dimRange>(rhs_func));
    //    RhsOperatorType rhs_operator(non_entropic_rhs_operator, hessian_inverter);
    CombinedOperatorType combined_operator(
        density_operator, non_entropic_advection_operator, non_entropic_rhs_operator, hessian_inverter);

    // ******************************** do the time steps ***********************************************************
    //    OperatorTimeStepperType timestepper_op(fv_operator, alpha, -1.0);
    //    RhsTimeStepperType timestepper_rhs(rhs_operator, alpha, 1.0);
    //    TimeStepperType timestepper(timestepper_op, timestepper_rhs);
    TimeStepperType timestepper(combined_operator, *analytical_flux, alpha, 1.);

    auto begin_time = std::chrono::steady_clock::now();
    auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1, [&basis_functions, &analytical_flux](const int /*comp*/, const auto& val) {
          return basis_functions->density(analytical_flux->get_u(val));
        });
    // auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
    //     1, [](const int /*comp*/, const auto& val) {
    //       double ret = 0.;
    //       for (const auto& entry : val)
    //         ret = std::max(std::abs(entry), ret);
    //       return ret;
    //     });

    timestepper.solve(t_end,
                      dt,
                      num_save_steps,
                      num_output_steps,
                      false,
                      true,
                      true,
                      false,
                      filename,
                      *visualizer,
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
struct HyperbolicEntropicCoordsMnTest
  : public HyperbolicEntropicCoordsMnDiscretization<TestCaseType>
  , public ::testing::Test
{
  void run()
  {
    auto norms = HyperbolicEntropicCoordsMnDiscretization<TestCaseType>::run(
                     DXTC_CONFIG.get("num_save_steps", 1),
                     0,
                     TestCaseType::quad_order,
                     TestCaseType::quad_refinements,
                     DXTC_CONFIG.get("grid_size", ""),
                     2,
                     DXTC_CONFIG.get("t_end", TestCaseType::t_end),
                     "test_kinetic_alpha",
                     Dune::GDT::is_full_moment_basis<typename TestCaseType::MomentBasis>::value)
                     .first;
    const double l1norm = norms[0];
    const double l2norm = norms[1];
    const double linfnorm = norms[2];
    using ResultsType = typename TestCaseType::ExpectedResultsType;
    EXPECT_NEAR(ResultsType::l1norm, l1norm, ResultsType::l1norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::l2norm, l2norm, ResultsType::l2norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::linfnorm, linfnorm, ResultsType::linfnorm * ResultsType::tol);
  }
};

#endif // DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_DISCRETIZATION_HH
