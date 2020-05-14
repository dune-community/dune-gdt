// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_EXAMPLES_MINIMUM_ENTROPY_TRANSFORMED_HH
#define DUNE_GDT_EXAMPLES_MINIMUM_ENTROPY_TRANSFORMED_HH

#ifndef ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING
#  define ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING 0
#endif

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
#include <dune/gdt/test/momentmodels/entropic-coords-mn-discretization.hh>
#include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/checkerboard.hh>
#include <dune/gdt/test/momentmodels/hessianinverter.hh>
#include <dune/gdt/test/momentmodels/density_evaluator.hh>
#include <dune/gdt/test/momentmodels/min_density_setter.hh>
#include <dune/gdt/tools/timestepper/adaptive-rungekutta-kinetic.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

using namespace Dune;
using namespace Dune::GDT;

template <size_t dimDomain>
class CoordinateTransformedBoltzmannSolver
{
  static_assert(dimDomain == 3, "Not yet implemented for other dimensions!");

public:
  // set dimensions
  static constexpr size_t momentOrder_or_refinements = (dimDomain == 1) ? 20 : 2;
  // choose basis
  using MomentBasis = HatFunctionMomentBasis<double, dimDomain, double, momentOrder_or_refinements, 1, dimDomain>;
  //   using MomentBasis = PartialMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3, 1>;
  // choose timestepper
  static constexpr TimeStepperMethods time_stepper_type = TimeStepperMethods::bogacki_shampine;

  static constexpr bool reconstruct = false;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using DomainFieldType = RangeFieldType;
  using GridType = YaspGrid<dimDomain, EquidistantOffsetCoordinates<DomainFieldType, dimDomain>>;
  using GV = typename GridType::LeafGridView;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using SpaceType = FiniteVolumeSpace<GV, dimRange, 1, RangeFieldType>;
  using AdvectionSourceSpaceType =
      std::conditional_t<reconstruct, DiscontinuousLagrangeSpace<GV, dimRange, RangeFieldType>, SpaceType>;
  static constexpr auto la_backend = Dune::XT::LA::Backends::common_dense;
  using MatrixType = typename XT::LA::Container<RangeFieldType, la_backend>::MatrixType;
  using VectorType = typename Dune::XT::LA::Container<RangeFieldType, la_backend>::VectorType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GV, dimRange, 1, RangeFieldType>;
  using GenericFunctionType = XT::Functions::GenericFunction<dimDomain, dimRange, 1, RangeFieldType>;
  using DomainType = FieldVector<RangeFieldType, dimDomain>;
  using RangeType = FieldVector<RangeFieldType, dimRange>;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using ProblemType = CheckerboardMn<GV, MomentBasis>;
  using AnalyticalFluxType = typename ProblemType::FluxType;
  static constexpr SlopeLimiterType slope = reconstruct ? SlopeLimiterType::minmod : SlopeLimiterType::no_slope;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GV, MomentBasis, slope>;
  using OldEntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  using DensityOperatorType = DensityEvaluator<MomentBasis, SpaceType, slope, MatrixType>;
  using MinDensitySetterType = MinDensitySetter<MomentBasis, SpaceType, slope, MatrixType>;
  using AdvectionOperatorType = AdvectionFvOperator<MatrixType, GV, dimRange>;
  using ReconstructionOperatorType =
      PointwiseLinearKineticReconstructionOperator<GV, EntropyFluxType, VectorType, RangeType>;
  using ReconstructionAdvectionOperatorType =
      AdvectionWithPointwiseReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
  using FvOperatorType = ReconstructionAdvectionOperatorType;
#if ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING
  using MasslumpedOperatorType = EntropicCoordinatesMasslumpedOperator<FvOperatorType, ProblemType>;
#else
  using HessianInverterType = EntropicHessianInverter<MomentBasis, SpaceType, slope, MatrixType>;
  using RhsOperatorType = LocalizableOperator<MatrixType, GV, dimRange>;
  using CombinedOperatorType =
      EntropicCoordinatesCombinedOperator<DensityOperatorType, FvOperatorType, RhsOperatorType, HessianInverterType>;
#endif
#if ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING
  using TimeStepperType = KineticAdaptiveRungeKuttaTimeStepper<MasslumpedOperatorType,
                                                               MinDensitySetterType,
                                                               DiscreteFunctionType,
                                                               EntropyFluxType,
                                                               time_stepper_type>;
#else
  using TimeStepperType = KineticAdaptiveRungeKuttaTimeStepper<CombinedOperatorType,
                                                               MinDensitySetterType,
                                                               DiscreteFunctionType,
                                                               EntropyFluxType,
                                                               time_stepper_type>;
#endif
  using BoundaryOperator = LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, VectorType, GV, dimRange>;
  using LambdaType = typename BoundaryOperator::LambdaType;
  using BoundaryFluxesMapType = std::map<DomainType, DynamicRangeType, XT::Common::FieldVectorFloatLess>;
  using KineticNumericalFluxType = NumericalKineticFlux<GV, MomentBasis, EntropyFluxType>;
  using BoundaryDistributionType = typename ProblemType::BoundaryDistributionType;
  using SolutionType = typename TimeStepperType::DiscreteSolutionType;
  using SolutionVectorsVectorType = std::vector<VectorType>;
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;

  CoordinateTransformedBoltzmannSolver(const std::string output_dir = "kinetic_transformed",
                                       const size_t num_save_steps = 10,
                                       const size_t grid_size = 20,
                                       const bool visualize_solution = true,
                                       const bool silent = false,
                                       const RangeFieldType sigma_s_scattering = 1,
                                       const RangeFieldType sigma_s_absorbing = 0,
                                       const RangeFieldType sigma_a_scattering = 0,
                                       const RangeFieldType sigma_a_absorbing = 10)
  {
    auto num_save_steps_copy = num_save_steps;
    if (num_save_steps > 1e6) // hack to allow for size_t(-1) when called from the python bindings
      num_save_steps_copy = size_t(-1);
    init(output_dir,
         num_save_steps_copy,
         grid_size,
         visualize_solution,
         silent,
         sigma_s_scattering,
         sigma_s_absorbing,
         sigma_a_scattering,
         sigma_a_absorbing);
  }

  double current_time() const
  {
    return timestepper_->current_time();
  }

  double t_end() const
  {
    return problem_->t_end();
  }

  void set_current_time(const double time)
  {
    timestepper_->current_time() = time;
  }

  void set_current_solution(const VectorType& vec)
  {
    timestepper_->current_solution().dofs().vector() = vec;
  }

  bool linear() const
  {
    return false;
  }

  SolutionVectorsVectorType solve()
  {
    if (!silent_)
      std::cout << "Solving... " << std::endl;
    auto begin_time = std::chrono::steady_clock::now();
    const auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1,
        [&](const int /*comp*/, const auto& val) { return basis_functions_->density(analytical_flux_->get_u(val)); });
    const auto u_stringifier = basis_functions_->stringifier();
    const auto stringifier = [&](const RangeType& val) { return u_stringifier(analytical_flux_->get_u(val)); };
    timestepper_->solve(t_end_,
                        initial_dt_,
                        num_save_steps_,
                        silent_ ? size_t(0) : size_t(-1),
                        true,
                        visualize_solution_,
                        DXTC_CONFIG_GET("write_txt", false),
                        false,
                        true,
                        filename_,
                        *visualizer,
                        stringifier);
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (!silent_)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;
    timestepper_->write_timings(filename_);

    std::vector<VectorType> ret;
    for (const auto& pair : timestepper_->solution())
      ret.push_back(pair.second.dofs().vector());
    return ret;
  }

  //   SolutionVectorsVectorType next_n_timesteps(const size_t n) const
  //   {
  //     if (!silent_)
  //       std::cout << "Calculating next " << XT::Common::to_string(n) << " time steps... " << std::endl;
  //     DXTC_TIMINGS.start("fv.solve");
  //     SolutionType solution;
  //     if (linear_)
  //       pn_timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
  //     else
  //       timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
  //     DXTC_TIMINGS.stop("fv.solve");
  //     if (!silent_)
  //       std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
  //     std::vector<VectorType> ret;
  //     for (const auto& pair : solution)
  //       ret.push_back(pair.second.dofs().vector());
  //     return ret;
  //   }

  VectorType apply_operator(const VectorType& source, const double time, const double dt) const
  {
    VectorType ret(source);
    combined_operator_->apply(source, ret, {{"t", {time}}, {"dt", {dt}}});
    return ret;
  }

  VectorType apply_restricted_operator(const VectorType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "");
    // VectorType ret(restricted_op_output_dofs_->size(), 0.);
    // RangeType u_entity;
    // RangeType u_neighbor;
    // const auto* mn_flux = dynamic_cast<const EntropyFluxType*>(flux_.get());

    // double min_acceptable_density = 1e-9;
    // size_t jj = 0;
    // for (size_t ii = 0; ii < restricted_op_entities_->size(); ++ii) {
    //   const auto& entity = (*restricted_op_entities_)[ii];
    //   RangeType ret_entity(0.), local_ret(0.);
    //   for (size_t kk = 0; kk < dimRange; ++kk)
    //     u_entity[kk] = source[jj * dimRange + kk];
    //   ++jj;
    //   basis_functions_->ensure_min_density(u_entity, min_acceptable_density);
    //   for (auto&& intersection : intersections(*grid_view_, entity)) {
    //     const auto intersection_center = intersection.geometry().center();
    //     if (intersection.neighbor()) {
    //       for (size_t kk = 0; kk < dimRange; ++kk)
    //         u_neighbor[kk] = source[jj * dimRange + kk];
    //       basis_functions_->ensure_min_density(u_neighbor, min_acceptable_density);
    //       ++jj;
    //     } else if (intersection.boundary()) {
    //       u_neighbor = boundary_values_->evaluate(intersection_center);
    //     } else {
    //       DUNE_THROW(MathError, "This should not happen!");
    //     }
    //     assert(intersection.indexInInside() >= 0);
    //     size_t direction = static_cast<size_t>(intersection.indexInInside()) / 2;
    //     const auto local_intersection_center = intersection.geometry().local(intersection_center);
    //     auto n_ij = intersection.unitOuterNormal(local_intersection_center);
    //     const auto neighbor = intersection.neighbor() ? intersection.outside() : entity;
    //     mn_flux->evaluate_kinetic_flux(entity, neighbor, u_entity, u_neighbor, n_ij, direction, local_ret);
    //     ret_entity += local_ret * intersection.geometry().integrationElement(local_intersection_center);
    //   } // intersections
    //   ret_entity /= entity.geometry().volume();
    //   for (const auto& pair : (*restricted_op_entity_dofs_to_output_dofs_)[ii]) {
    //     assert(ret[pair.second] == 0. && "More than one output dofs maps to same vector location!");
    //     ret[pair.second] = ret_entity[pair.first];
    //   }
    // } // entities
    // return ret;
  }

  //  bool is_realizable(const VectorType& vector) const
  //  {
  //    const DiscreteFunctionType discrete_function(*fv_space_, vector);
  //    for (const auto& entity : elements(fv_space_->grid_view())) {
  //      const auto local_vec =
  //      discrete_function.local_function(entity)->evaluate(entity.geometry.local(entity.geometry().center())); if
  //      (!basis_functions_.is_realizable(local_vec))
  //        return false;
  //    }
  //    return true;
  //  }

  void set_parameters(const RangeFieldType sigma_s_scattering = 1,
                      const RangeFieldType sigma_s_absorbing = 0,
                      const RangeFieldType sigma_a_scattering = 0,
                      const RangeFieldType sigma_a_absorbing = 10)
  {
    static const RangeType u_iso = basis_functions_->u_iso();
    static const RangeType basis_integrated = basis_functions_->integrated();
    std::shared_ptr<ParameterFunctionType> sigma_s(
        problem_->create_parameter_function(sigma_s_absorbing, sigma_s_scattering, sigma_s_scattering));
    std::shared_ptr<ParameterFunctionType> sigma_a(
        problem_->create_parameter_function(sigma_a_absorbing, sigma_a_scattering, sigma_a_scattering));
    std::shared_ptr<ParameterFunctionType> Q(problem_->Q());
    auto rhs_func = [&, sigma_a, sigma_s, Q](const auto& /*source*/,
                                             const auto& /*local_source*/,
                                             auto& local_range,
                                             const XT::Common::Parameter& /*param*/) {
      const auto& element = local_range.element();
      const auto center = element.geometry().center();
      const auto& u_elem = analytical_flux_->get_precomputed_u(grid_view_->indexSet().index(element));
      const auto sigma_a_value = sigma_a->evaluate(center)[0];
      const auto sigma_s_value = sigma_s->evaluate(center)[0];
      const auto sigma_t_value = sigma_a_value + sigma_s_value;
      const auto Q_value = Q->evaluate(center)[0];
      auto ret = u_elem;
      ret *= -sigma_t_value;
      ret.axpy(basis_functions_->density(u_elem) * sigma_s_value, u_iso);
      ret.axpy(Q_value, basis_integrated);
      auto& range_dofs = local_range.dofs();
      for (size_t ii = 0; ii < dimRange; ++ii)
        range_dofs.add_to_entry(ii, ret[ii]);
    };
    rhs_operator_ = std::make_shared<RhsOperatorType>(*grid_view_, *fv_space_, *fv_space_, false, false);
    rhs_operator_->append(GenericLocalElementOperator<VectorType, GV, dimRange>(rhs_func));
  }

  void create_rhs_operator(const RangeFieldType sigma_s_scattering = 1,
                           const RangeFieldType sigma_s_absorbing = 0,
                           const RangeFieldType sigma_a_scattering = 0,
                           const RangeFieldType sigma_a_absorbing = 10)
  {
    return set_parameters(sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
  }

  VectorType get_initial_values() const
  {
    return *initial_values_alpha_;
  }

  bool finished() const
  {
    return XT::Common::FloatCmp::ge(timestepper_->current_time(), t_end_);
  }

  void init(const std::string output_dir = "boltzmann",
            const size_t num_save_steps = 10,
            const size_t grid_size = 50,
            const bool visualize_solution = true,
            const bool silent = false,
            const double sigma_s_scattering = 1.,
            const double sigma_s_absorbing = 0.,
            const double sigma_a_scattering = 0.,
            const double sigma_a_absorbing = 10.)
  {
    silent_ = silent;
    if (!silent_)
      std::cout << "Setting problem parameters ...";
    visualize_solution_ = visualize_solution;
    filename_ = output_dir;
    num_save_steps_ = num_save_steps;
    // create grid
    auto grid_config = ProblemType::default_grid_cfg();
    grid_config["num_elements"] = XT::Common::to_string(grid_size);
    grid_ =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_->comm().size() == 1 || grid_->overlapSize(0) > 0);
    grid_view_ = std::make_shared<GV>(grid_->leafGridView());

    // create space
    fv_space_ = std::make_shared<const SpaceType>(*grid_view_);

    basis_functions_ = std::make_shared<const MomentBasis>();
    static const RangeFieldType psi_vac = DXTC_CONFIG_GET("psi_vac", 1e-6 / basis_functions_->unit_ball_volume());
    problem_ = std::make_shared<const ProblemType>(*basis_functions_, *grid_view_, psi_vac, grid_config, true, 1e-09);

    // fluxes
    flux_ = std::shared_ptr<AnalyticalFluxType>(problem_->flux());
    auto* entropy_flux = dynamic_cast<OldEntropyFluxType*>(flux_.get());
    analytical_flux_ = std::make_shared<EntropyFluxType>(*entropy_flux);

    // create a discrete function for the solution
    u_ = std::make_shared<DiscreteFunctionType>(*fv_space_, "u_initial");
    // The only operator that needs synchronisation is the advection operator
    alpha_vec_ = std::make_shared<VectorType>(fv_space_->mapper().size(), 0., DXTC_CONFIG_GET("num_alpha_mutexes", 1));
    alpha_ = std::make_shared<DiscreteFunctionType>(*fv_space_, *alpha_vec_, "alpha");
    // project initial values
    const auto initial_values_u = problem_->initial_values();
    if (!silent_)
      std::cout << "Projecting initial u values...";
    default_interpolation(*initial_values_u, *u_, *grid_view_);
    if (!silent_) {
      std::cout << " done " << std::endl;
      std::cout << "Transforming initial u values to alpha...";
    }
    initial_values_alpha_ = std::make_shared<VectorType>(*alpha_vec_);
    const auto u_local_func = u_->local_discrete_function();
    const auto alpha_local_func = alpha_->local_discrete_function();
    const auto entropy_flux_local_func = entropy_flux->derived_local_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> u_local;
    for (auto&& element : Dune::elements(*grid_view_)) {
      u_local_func->bind(element);
      alpha_local_func->bind(element);
      entropy_flux_local_func->bind(element);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_local[ii] = u_local_func->dofs().get_entry(ii);
      const auto alpha_local = entropy_flux_local_func->get_alpha(u_local, false)->first;
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_local_func->dofs().set_entry(ii, alpha_local[ii]);
    }
    if (!silent_)
      std::cout << " done " << std::endl;

    const auto first_entity = *grid_view_->template begin<0>();
    const auto first_intersection = *grid_view_->ibegin(first_entity);
    dx_ = first_entity.geometry().volume() / first_intersection.geometry().volume();
    initial_dt_ = dx_ / 100.;
    t_end_ = problem_->t_end();

    // create operators
    numerical_flux_ = std::make_shared<KineticNumericalFluxType>(*analytical_flux_, *basis_functions_);
    advection_source_space_ = std::make_shared<AdvectionSourceSpaceType>(*grid_view_);
    advection_operator_ = std::make_shared<AdvectionOperatorType>(
        *grid_view_, *numerical_flux_, *advection_source_space_, *fv_space_, /*use_tbb*/ false);
    reconstruction_operator_ = std::make_shared<ReconstructionOperatorType>(*fv_space_, *analytical_flux_);
    fv_operator_ =
        std::make_shared<ReconstructionAdvectionOperatorType>(*advection_operator_, *reconstruction_operator_);
    const double min_acceptable_density =
        DXTC_CONFIG_GET("rho_min", problem_->psi_vac() * basis_functions_->unit_ball_volume() / 10);
    boundary_distribution_ = std::make_shared<BoundaryDistributionType>(problem_->boundary_distribution());
    density_operator_ = std::make_shared<DensityOperatorType>(
        *analytical_flux_, *fv_space_, *boundary_distribution_, min_acceptable_density);
    min_density_setter_ = std::make_shared<MinDensitySetterType>(*analytical_flux_, *fv_space_, min_acceptable_density);

#if ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING
    MasslumpedOperatorType masslumped_operator(fv_operator, problem, dx, boundary_fluxes);
#else
    hessian_inverter_ = std::make_shared<HessianInverterType>(*analytical_flux_, *fv_space_);
    create_rhs_operator(sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
    combined_operator_ =
        std::make_shared<CombinedOperatorType>(*density_operator_, *fv_operator_, *rhs_operator_, *hessian_inverter_);
#endif

    // boundary treatment
    // store boundary fluxes
    if (!silent_)
      std::cout << "Calculating boundary fluxes... ";
    boundary_fluxes_ = std::make_shared<BoundaryFluxesMapType>();
    BoundaryFluxesFunctor<GV, ProblemType, BoundaryFluxesMapType, EntropyFluxType> boundary_flux_functor(
        *problem_, *boundary_fluxes_, *analytical_flux_, grid_view_->indexSet());
    auto walker = XT::Grid::Walker<GV>(*grid_view_);
    XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GV> boundary_intersection_filter;
    walker.append(boundary_flux_functor, boundary_intersection_filter);
    walker.walk(false);
    boundary_kinetic_fluxes_ = std::make_shared<GenericFunctionType>(
        1, [this](const DomainType& x, DynamicRangeType& ret, const XT::Common::Parameter&) {
          ret = (*boundary_fluxes_)[x];
        });
    boundary_lambda_ = std::make_shared<LambdaType>(
        [this](const I& intersection,
               const FieldVector<RangeFieldType, dimDomain - 1>& xx_in_reference_intersection_coordinates,
               const AnalyticalFluxType& /*flux*/,
               const DynamicRangeType& /*u*/,
               DynamicRangeType& v,
               const XT::Common::Parameter& param) {
          boundary_kinetic_fluxes_->evaluate(
              intersection.geometry().global(xx_in_reference_intersection_coordinates), v, param);
        });
    advection_operator_->append(*boundary_lambda_, {}, boundary_intersection_filter);
    if (!silent_)
      std::cout << " done " << std::endl;

    // create timestepper
    constexpr double default_tol = dimDomain == 1 ? 1e-3 : 1e-2;
    atol_ = DXTC_CONFIG.get("timestepper.atol", default_tol);
    rtol_ = DXTC_CONFIG.get("timestepper.rtol", default_tol);
#if ENTROPY_FLUX_HATFUNCTIONS_USE_MASSLUMPING
    TimeStepperType timestepper(masslumped_operator, min_density_setter, *analytical_flux, alpha, true, 1., atol, rtol);
#else
    timestepper_ = std::make_shared<TimeStepperType>(
        *combined_operator_, *min_density_setter_, *analytical_flux_, *alpha_, true, 1., atol_, rtol_);
#endif
  } // void init()

  void reset()
  {
    *alpha_vec_ = get_initial_values();
    alpha_ = std::make_shared<DiscreteFunctionType>(*fv_space_, *alpha_vec_, "alpha");
    timestepper_ = std::make_shared<TimeStepperType>(
        *combined_operator_, *min_density_setter_, *analytical_flux_, *alpha_, true, 1., atol_, rtol_);
    TimeStepperType::reset_static_variables();
  }

  void prepare_restricted_operator(const std::vector<size_t>& /*output_dofs*/)
  {
    DUNE_THROW(NotImplemented, "");
    // if (!restricted_op_output_dofs_ || *restricted_op_output_dofs_ != output_dofs) {
    //   restricted_op_output_dofs_ = std::make_shared<std::vector<size_t>>(output_dofs);
    //   restricted_op_input_dofs_ = std::make_shared<std::vector<size_t>>();
    //   restricted_op_entity_dofs_to_output_dofs_ = std::make_shared<std::vector<std::map<size_t, size_t>>>();
    //   const auto& mapper = fv_space_->mapper();
    //   DynamicVector<size_t> global_dofs_entity;
    //   DynamicVector<size_t> global_dofs_neighbor;
    //   // calculate entities corresponding to dofs in restricted operator
    //   restricted_op_entities_ = std::make_shared<std::vector<E>>();
    //   for (auto&& entity : elements(*grid_view_)) {
    //     // we only want to add each entity once, but there may be several output_dofs per entity
    //     bool entity_added = false;
    //     mapper.global_indices(entity, global_dofs_entity);
    //     // check if any output dof matches a dof on this entity
    //     for (size_t ll = 0; ll < output_dofs.size(); ++ll) {
    //       for (size_t kk = 0; kk < global_dofs_entity.size(); ++kk) {
    //         if (global_dofs_entity[kk] == output_dofs[ll]) {
    //           if (!entity_added) {
    //             restricted_op_entities_->push_back(entity);
    //             for (auto&& global_dof : global_dofs_entity)
    //               restricted_op_input_dofs_->push_back(global_dof);
    //             for (auto&& intersection : intersections(*grid_view_, entity)) {
    //               if (intersection.neighbor()) {
    //                 mapper.global_indices(intersection.outside(), global_dofs_neighbor);
    //                 for (auto&& global_dof : global_dofs_neighbor)
    //                   restricted_op_input_dofs_->push_back(global_dof);
    //               }
    //             } // intersections
    //             restricted_op_entity_dofs_to_output_dofs_->emplace_back();
    //             entity_added = true;
    //           } // if (!entity_added)
    //           restricted_op_entity_dofs_to_output_dofs_->back().insert(std::make_pair(kk, ll));
    //         } // if (output dof found)
    //       } // kk
    //     } // ll
    //   } // entities
    // }
  }

  // std::vector<size_t> restricted_op_input_dofs() const
  // {
  //   return *restricted_op_input_dofs_;
  // }

  // size_t restricted_op_input_dofs_size() const
  // {
  //   return restricted_op_input_dofs_->size();
  // }

private:
  std::shared_ptr<const GridType> grid_;
  std::shared_ptr<const GV> grid_view_;
  std::shared_ptr<const MomentBasis> basis_functions_;
  std::shared_ptr<const ProblemType> problem_;
  std::shared_ptr<const SpaceType> fv_space_;
  std::shared_ptr<const AdvectionSourceSpaceType> advection_source_space_;
  std::shared_ptr<DiscreteFunctionType> u_;
  std::shared_ptr<VectorType> alpha_vec_;
  std::shared_ptr<VectorType> initial_values_alpha_;
  std::shared_ptr<DiscreteFunctionType> alpha_;
  std::shared_ptr<KineticNumericalFluxType> numerical_flux_;
  std::shared_ptr<AdvectionOperatorType> advection_operator_;
  std::shared_ptr<ReconstructionOperatorType> reconstruction_operator_;
  std::shared_ptr<ReconstructionAdvectionOperatorType> fv_operator_;
  std::shared_ptr<RhsOperatorType> rhs_operator_;
  std::shared_ptr<HessianInverterType> hessian_inverter_;
  std::shared_ptr<BoundaryDistributionType> boundary_distribution_;
  std::shared_ptr<DensityOperatorType> density_operator_;
  std::shared_ptr<MinDensitySetterType> min_density_setter_;
  std::shared_ptr<CombinedOperatorType> combined_operator_;
  // std::shared_ptr<std::vector<E>> restricted_op_entities_;
  // std::shared_ptr<std::vector<size_t>> restricted_op_input_dofs_;
  // std::shared_ptr<std::vector<size_t>> restricted_op_output_dofs_;
  // std::shared_ptr<std::vector<std::map<size_t, size_t>>> restricted_op_entity_dofs_to_output_dofs_;
  std::shared_ptr<TimeStepperType> timestepper_;
  std::shared_ptr<AnalyticalFluxType> flux_;
  std::shared_ptr<EntropyFluxType> analytical_flux_;
  std::shared_ptr<BoundaryFluxesMapType> boundary_fluxes_;
  std::shared_ptr<LambdaType> boundary_lambda_;
  std::shared_ptr<GenericFunctionType> boundary_kinetic_fluxes_;
  double t_end_;
  double initial_dt_;
  double current_dt_;
  double dx_;
  bool silent_;
  bool visualize_solution_;
  std::string filename_;
  size_t num_save_steps_;
  double atol_;
  double rtol_;
};

#endif // DUNE_GDT_EXAMPLES_MINIMUM_ENTROPY_TRANSFORMED_HH
