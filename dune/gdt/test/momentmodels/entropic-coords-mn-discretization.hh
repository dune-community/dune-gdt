// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_DISCRETIZATION_HH

#include <chrono>

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

template <class GV, class ProblemType, class MapType, class EntropyFluxType>
class BoundaryFluxesFunctor : public Dune::XT::Grid::IntersectionFunctor<GV>
{
  using BaseType = typename Dune::XT::Grid::IntersectionFunctor<GV>;
  using IndexSetType = typename GV::IndexSet;

public:
  using typename BaseType::E;
  using typename BaseType::I;

  BoundaryFluxesFunctor(const ProblemType& problem,
                        MapType& boundary_fluxes_map,
                        EntropyFluxType& analytical_flux,
                        const IndexSetType& index_set)
    : problem_(problem)
    , analytical_flux_(analytical_flux)
    , index_set_(index_set)
    , boundary_distribution_(problem.boundary_distribution())
    , boundary_fluxes_map_(boundary_fluxes_map)
    , mutex_(std::make_shared<std::mutex>())
  {}

  BoundaryFluxesFunctor(const BoundaryFluxesFunctor& other)
    : BaseType(other)
    , problem_(other.problem_)
    , analytical_flux_(other.analytical_flux_)
    , index_set_(other.index_set_)
    , boundary_distribution_(other.boundary_distribution_)
    , boundary_fluxes_map_(other.boundary_fluxes_map_)
    , mutex_(other.mutex_)
  {}

  Dune::XT::Grid::IntersectionFunctor<GV>* copy() override final
  {
    return new BoundaryFluxesFunctor(*this);
  }

  virtual void
  apply_local(const I& intersection, const E& /*inside_element*/, const E& /*outside_element*/) override final
  {
    if (intersection.boundary()) {
      // store boundary fluxes
      const auto x = intersection.geometry().center();
      const auto dd = intersection.indexInInside() / 2;
      const auto n = intersection.centerUnitOuterNormal()[dd];
      auto boundary_flux = problem_.kinetic_boundary_flux(x, n, dd);
      // The boundary_flux calculates <psi b (v[dd]*n)>, we only want to have <psi b v[dd]> because the
      // multiplication with n is done in entropy_flux->evaluate_kinetic_flux(..)
      boundary_flux *= n;
      mutex_->lock();
      boundary_fluxes_map_.insert(std::make_pair(x, boundary_flux));
      mutex_->unlock();
      // store boundary evaluations
      analytical_flux_.store_boundary_evaluations(
          boundary_distribution_(x), index_set_.index(intersection.inside()), intersection.indexInInside());
    }
  }

private:
  const ProblemType& problem_;
  EntropyFluxType& analytical_flux_;
  const IndexSetType& index_set_;
  const typename ProblemType::BoundaryDistributionType boundary_distribution_;
  MapType& boundary_fluxes_map_;
  std::shared_ptr<std::mutex> mutex_;
};


template <class TestCaseType>
struct HyperbolicEntropicCoordsMnDiscretization
{
  // returns: (l1norm, l2norm, linfnorm, MPI rank)
  static std::pair<Dune::FieldVector<double, 3>, int> run(size_t num_save_steps = 1,
                                                          size_t num_output_steps = 0,
                                                          size_t quad_order = TestCaseType::quad_order,
                                                          size_t quad_refinements = TestCaseType::quad_refinements,
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
    // using E = XT::Grid::extract_entity_t<GV>;
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
    std::shared_ptr<const MomentBasis> basis_functions =
        std::make_shared<const MomentBasis>(quad_order, quad_refinements);
    const RangeFieldType psi_vac = DXTC_CONFIG_GET("psi_vac", 1e-6 / basis_functions->unit_ball_volume());
    const std::unique_ptr<ProblemType> problem_ptr =
        std::make_unique<ProblemType>(*basis_functions, grid_view, psi_vac, grid_config);
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
    // entropy_flux->disable_thread_cache();
    auto analytical_flux = std::make_unique<EntropyFluxType>(*entropy_flux);

    // calculate boundary values for alpha
    // std::map<DomainType, RangeType, XT::Common::FieldVectorFloatLess> alpha_boundary_vals;
    // for (const auto& element : Dune::elements(grid_view))
    //   for (const auto& intersection : Dune::intersections(grid_view, element))
    //     if (intersection.boundary()) {
    //       const auto x = intersection.geometry().center();
    //       const auto u = boundary_values_u->evaluate(x);
    //       alpha_boundary_vals.insert(std::make_pair(x, analytical_flux->get_alpha(u)));
    //       std::cout << XT::Common::to_string(alpha_boundary_vals[x]) << std::endl;
    //     }
    // GenericFunctionType boundary_values_alpha(
    //     1, [&](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) {
    //       ret = alpha_boundary_vals[x];
    //     });

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "u_initial");
    // The only operator that needs synchronisation is the advection operator
    VectorType alpha_vec(fv_space.mapper().size(), 0., DXTC_CONFIG_GET("num_alpha_mutexes", 1));
    DiscreteFunctionType alpha(fv_space, alpha_vec, "alpha_initial");
    // project initial values
    default_interpolation(*initial_values_u, u, grid_view);

    const auto u_local_func = u.local_discrete_function();
    const auto alpha_local_func = alpha.local_discrete_function();
    const auto entropy_flux_local_func = entropy_flux->derived_local_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> u_local;
    for (auto&& element : Dune::elements(grid_view)) {
      u_local_func->bind(element);
      alpha_local_func->bind(element);
      entropy_flux_local_func->bind(element);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_local[ii] = u_local_func->dofs().get_entry(ii);
      const auto alpha_local = entropy_flux_local_func->get_alpha(u_local, false)->first;
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_local_func->dofs().set_entry(ii, alpha_local[ii]);
    }

    // enforce min acceptable density for initial values
    const double min_acceptable_density = problem.psi_vac() * basis_functions->unit_ball_volume() / 10;
    // const double min_acceptable_density = problem.psi_vac();
    using DensityOperatorType = DensityEvaluator<MomentBasis, SpaceType, slope, MatrixType>;
    using MinDensitySetterType = MinDensitySetter<MomentBasis, SpaceType, slope, MatrixType>;
    DensityOperatorType density_operator(*analytical_flux, fv_space, boundary_distribution, min_acceptable_density);
    MinDensitySetterType min_density_setter(*analytical_flux, fv_space, min_acceptable_density);
    min_density_setter.apply(alpha.dofs().vector(), alpha.dofs().vector());

    // ******************** choose flux and rhs operator and timestepper ******************************************

    using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
    using HessianInverterType = EntropicHessianInverter<MomentBasis, SpaceType, slope, MatrixType>;
    using ReconstructionOperatorType =
        PointwiseLinearKineticReconstructionOperator<GV, EntropyFluxType, VectorType, RangeType>;
    using ReconstructionAdvectionOperatorType =
        AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
    using FvOperatorType = ReconstructionAdvectionOperatorType;
    using RhsOperatorType = LocalizableOperator<MatrixType, GV, dimRange>;
    using CombinedOperatorType =
        EntropicCoordinatesCombinedOperator<DensityOperatorType, FvOperatorType, RhsOperatorType, HessianInverterType>;

    constexpr TimeStepperMethods time_stepper_type = TimeStepperMethods::bogacki_shampine;
    // constexpr TimeStepperMethods time_stepper_type = TimeStepperMethods::dormand_prince;
    using TimeStepperType = KineticAdaptiveRungeKuttaTimeStepper<CombinedOperatorType,
                                                                 MinDensitySetterType,
                                                                 DiscreteFunctionType,
                                                                 EntropyFluxType,
                                                                 time_stepper_type>;

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);

    // *********************** create operators and timesteppers ************************************
    NumericalKineticFlux<GV, MomentBasis, EntropyFluxType> numerical_flux(*analytical_flux, *basis_functions);
    // do not use parallelisation here, as the advection operator does almost no work (allows to use alpha_vec without
    // mutexes)
    AdvectionOperatorType advection_operator(
        grid_view, numerical_flux, advection_source_space, fv_space, /*use_tbb*/ true);

    // boundary treatment
    using BoundaryOperator =
        LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
    using LambdaType = typename BoundaryOperator::LambdaType;

    // store boundary fluxes
    using BoundaryFluxesMapType = std::map<DomainType, DynamicRangeType, XT::Common::FieldVectorFloatLess>;
    BoundaryFluxesMapType boundary_fluxes;
    BoundaryFluxesFunctor<GV, ProblemType, BoundaryFluxesMapType, EntropyFluxType> boundary_flux_functor(
        problem, boundary_fluxes, *analytical_flux, fv_space.grid_view().indexSet());
    auto walker = XT::Grid::Walker<GV>(fv_space.grid_view());
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> boundary_intersection_filter;
    walker.append(boundary_flux_functor, boundary_intersection_filter);
    walker.walk(true);
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
    advection_operator.append(boundary_lambda, {}, boundary_intersection_filter);

    // create remaining operators
    ReconstructionOperatorType reconstruction_operator(fv_space, *analytical_flux);
    ReconstructionAdvectionOperatorType reconstruction_advection_operator(advection_operator, reconstruction_operator);
    auto& fv_operator = reconstruction_advection_operator;

    const double atol = DXTC_CONFIG.get("timestepper.atol", 1e-3);
    const double rtol = DXTC_CONFIG.get("timestepper.rtol", 1e-3);
    if (XT::Common::is_zero(t_end))
      t_end = problem.t_end();

    if (!filename.empty())
      filename += "_";
    filename += ProblemType::static_id();
    if (TestCaseType::reconstruction && slope == SlopeLimiterType::minmod)
      filename += "_minmod_";
    else if (TestCaseType::reconstruction && slope == SlopeLimiterType::superbee)
      filename += "_superbee_";
    filename += (time_stepper_type == TimeStepperMethods::bogacki_shampine)
                    ? "rk23"
                    : (time_stepper_type == TimeStepperMethods::dormand_prince ? "rk45" : "unknown");
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad_" + XT::Common::to_string(quad_refinements) + "x" + XT::Common::to_string(quad_order);
    filename += "_atol_" + XT::Common::to_string(atol);
    filename += "_rtol_" + XT::Common::to_string(rtol);
    filename += "_threads_" + DXTC_CONFIG.get("threading.max_count", "1") + "x"
                + DXTC_CONFIG.get("threading.partition_factor", "1");
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->mn_name();

    HessianInverterType hessian_inverter(*analytical_flux, fv_space);

    static const RangeType u_iso = basis_functions->u_iso();
    static const RangeType basis_integrated = basis_functions->integrated();
    const auto sigma_a = problem.sigma_a();
    const auto sigma_s = problem.sigma_s();
    const auto Q = problem.Q();
    static std::mutex mutex;
    auto rhs_func = [&](const auto& /*source*/,
                        const auto& /*local_source*/,
                        auto& local_range,
                        const Dune::XT::Common::Parameter& /*param*/) {
      const auto& element = local_range.element();
      const auto center = element.geometry().center();
      const auto& u_elem = analytical_flux->get_precomputed_u(fv_space.grid_view().indexSet().index(element));
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
        range_dofs.add_to_entry(ii, ret[ii]);
    };
    RhsOperatorType rhs_operator(grid_view, fv_space, fv_space, false, true);
    rhs_operator.append(GenericLocalElementOperator<VectorType, GV, dimRange>(rhs_func));
    CombinedOperatorType combined_operator(density_operator, fv_operator, rhs_operator, hessian_inverter);

    // ******************************** do the time steps ***********************************************************
    TimeStepperType timestepper(combined_operator, min_density_setter, *analytical_flux, alpha, true, 1., atol, rtol);

    auto begin_time = std::chrono::steady_clock::now();
    auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1, [&basis_functions, &analytical_flux](const int /*comp*/, const auto& val) {
          return basis_functions->density(analytical_flux->get_u(val));
        });
    const auto u_stringifier = basis_functions->stringifier();
    const auto stringifier = [&u_stringifier, &analytical_flux](const RangeType& val) {
      return u_stringifier(analytical_flux->get_u(val));
    };
    // auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
    //     1, [](const int /*comp*/, const auto& val) {
    //       double ret = 0.;
    //       for (const auto& entry : val)
    //         ret = std::max(std::abs(entry), ret);
    //       return ret;
    //     });

    // The hessian has entries in the order of psi_min, the inverse thus scales with 1/psi_min, and thus the timestep
    // should be psi_min to get an update of order 1
    double initial_dt = dx / 100.; // std::min(dt, min_acceptable_density);
    timestepper.solve(t_end,
                      initial_dt,
                      num_save_steps,
                      num_output_steps,
                      false,
                      true,
                      true,
                      false,
                      true,
                      filename,
                      *visualizer,
                      stringifier);
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
struct HyperbolicEntropicCoordsMnTest
  : public HyperbolicEntropicCoordsMnDiscretization<TestCaseType>
  , public ::testing::Test
{
  void run()
  {
    auto norms = HyperbolicEntropicCoordsMnDiscretization<TestCaseType>::run(
                     DXTC_CONFIG.get("num_save_steps", 1),
                     DXTC_CONFIG.get("num_output_steps", 0),
                     DXTC_CONFIG.get("quad_order", TestCaseType::quad_order),
                     DXTC_CONFIG.get("quad_refinements", TestCaseType::quad_refinements),
                     DXTC_CONFIG.get("grid_size", ""),
                     DXTC_CONFIG.get("overlap_size", 2),
                     DXTC_CONFIG.get("t_end", TestCaseType::t_end),
                     DXTC_CONFIG.get("filename", "timings_kinetic"))
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
