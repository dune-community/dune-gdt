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
#include <dune/xt/common/exceptions.hh>
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
#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>
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

//! struct to be used as comparison function e.g. in a std::map<Entity, EntityLess>
template <class GV>
struct EntityLess
{
  using IndexSet = typename GV::IndexSet;
  using E = typename GV::Grid::template Codim<0>::Entity;

  EntityLess(const IndexSet& index_set)
    : index_set_(index_set)
  {}

  bool operator()(const E& a, const E& b) const
  {
    return index_set_.index(a) < index_set_.index(b);
  }

  const IndexSet& index_set_;
};


template <class ProblemType,
          class GridType =
              YaspGrid<ProblemType::dimDomain, EquidistantOffsetCoordinates<double, ProblemType::dimDomain>>>
class CoordinateTransformedMnSolver
{
public:
  static constexpr size_t dimDomain = ProblemType::dimDomain;

  // choose timestepper
  static constexpr TimeStepperMethods time_stepper_type = TimeStepperMethods::bogacki_shampine;
  static constexpr bool reconstruct = false;
  using MomentBasis = typename ProblemType::MomentBasis;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using DomainFieldType = RangeFieldType;
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
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GV, dimRange, 1, RangeFieldType>;
  using GenericFunctionType = XT::Functions::GenericFunction<dimDomain, dimRange, 1, RangeFieldType>;
  using DomainType = FieldVector<RangeFieldType, dimDomain>;
  using RangeType = FieldVector<RangeFieldType, dimRange>;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
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
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;

  static std::vector<double> default_parameters()
  {
    return (dimDomain == 3) ? std::vector<double>{1, 0, 0, 10} : std::vector<double>{1, 0, 0, 2, 10};
  }

  // The parameters vector should contain:
  // - Sourcebeam test: (sigma_a_left, sigma_a_right, sigma_s_left, sigma_s_middle, sigma_s_right)
  // - Checkerboard test: (sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing)
  CoordinateTransformedMnSolver(const std::string output_dir = "kinetic_transformed",
                                const size_t num_save_steps = 10,
                                const size_t grid_size = DXTC_CONFIG_GET("grid_size", 21),
                                const bool visualize_solution = true,
                                const bool silent = false,
                                const std::vector<double>& parameters = default_parameters())
  {
    auto num_save_steps_copy = num_save_steps;
    if (num_save_steps > 1e6) // hack to allow for size_t(-1) when called from the python bindings
      num_save_steps_copy = size_t(-1);
    init(output_dir, num_save_steps_copy, grid_size, visualize_solution, silent, parameters);
  }

  double current_time() const
  {
    return timestepper_->current_time();
  }

  double t_end() const
  {
    return t_end_;
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

  // First vector in returned pair contains the time points, the second vector contains the solution vectors
  std::pair<std::vector<double>, std::vector<VectorType>> solve()
  {
    if (!silent_)
      std::cout << "Solving... " << std::endl;
    auto begin_time = std::chrono::steady_clock::now();
    const auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1,
        [&](const int /*comp*/, const auto& val) { return basis_functions_->density(analytical_flux_->get_u(val)); });
    const auto u_stringifier = basis_functions_->stringifier();
    const auto stringifier = [&](const RangeType& val) { return u_stringifier(analytical_flux_->get_u(val)); };
    typename TimeStepperType::DiscreteSolutionType solution;
    timestepper_->solve(t_end_,
                        initial_dt_,
                        num_save_steps_,
                        silent_ ? size_t(0) : size_t(-1),
                        true,
                        visualize_solution_,
                        false,
                        DXTC_CONFIG_GET("write_txt", false),
                        false,
                        true,
                        filename_,
                        solution,
                        *visualizer,
                        stringifier,
                        timestepper_->dummy_solution());
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (!silent_)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;
    if (visualize_solution_)
      timestepper_->write_timings(filename_);

    // const auto interval_length = t_end_ / num_intervals;
    // for (size_t ii = 0; ii < num_intervals; ++ii) {
    //   const auto left_boundary = ii * interval_length;
    //   const auto right_boundary = (ii+1) * interval_length;
    //   auto begin_it = solution.lower_bound(left_boundary);
    //   const auto end_it = solution.upper_bound(right_boundary);
    //   const auto num_elements_in_interval = std::distance(begin_it, end_it);
    //   if (num_elements_in_interval > max_vectors_per_interval) {
    std::pair<std::vector<double>, std::vector<VectorType>> ret;
    for (auto it = solution.begin(); it != solution.end();) {
      ret.first.push_back(it->first);
      ret.second.push_back(it->second.dofs().vector());
      it = solution.erase(it);
    }
    return ret;
  }

  //   std::vector<VectorType> next_n_timesteps(const size_t n) const
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

  VectorType apply_operator(const VectorType& source) const
  {
    try {
      VectorType ret(source);
      static const XT::Common::Parameter unused_parameter;
      combined_operator_->apply(source, ret, unused_parameter);
      return ret;
    } catch (const Dune::MathError&) {
      return VectorType{};
    }
  }

  VectorType apply_restricted_operator(const VectorType& source) const
  {
    static const XT::Common::Parameter unused_parameter;
    try {
      // copy source values to full-dimensional vector
      const auto& input_dofs = *restricted_op_input_dofs_;
      for (size_t kk = 0; kk < input_dofs.size(); ++kk)
        alpha_vec_->set_entry(input_dofs[kk], source.get_entry(kk));
      // apply restricted operator
      combined_operator_->apply_range(*alpha_vec_,
                                      *tmp_fulldimensional_vec_,
                                      unused_parameter,
                                      *restricted_op_output_entities_,
                                      *restricted_op_input_entities_);
      // copy output_dofs to return vector
      const auto& output_dofs = *restricted_op_output_dofs_;
      VectorType ret(output_dofs.size(), 0.);
      for (size_t jj = 0; jj < output_dofs.size(); ++jj)
        ret[jj] = (*tmp_fulldimensional_vec_)[output_dofs[jj]];
      return ret;
    } catch (const Dune::MathError&) {
      return VectorType{};
    }
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

  void set_parameters(const std::vector<double>& parameters)
  {
    static const RangeType u_iso = basis_functions_->u_iso();
    static const RangeType basis_integrated = basis_functions_->integrated();
    std::shared_ptr<ParameterFunctionType> sigma_s;
    std::shared_ptr<ParameterFunctionType> sigma_a;
    if constexpr (dimDomain == 1) {
      DUNE_THROW_IF(parameters.size() != 5, Dune::InvalidStateException, "Wrong parameter size!");
      const double sigma_a_left = parameters[0];
      const double sigma_a_right = parameters[1];
      const double sigma_s_left = parameters[2];
      const double sigma_s_middle = parameters[3];
      const double sigma_s_right = parameters[4];
      sigma_s = std::shared_ptr(problem_->create_sigma_s_function(sigma_s_left, sigma_s_middle, sigma_s_right));
      sigma_a = std::shared_ptr(problem_->create_sigma_a_function(sigma_a_left, sigma_a_right));
    } else {
      DUNE_THROW_IF(parameters.size() != 4, Dune::InvalidStateException, "Wrong parameter size!");
      const double sigma_s_scattering = parameters[0];
      const double sigma_s_absorbing = parameters[1];
      const double sigma_a_scattering = parameters[2];
      const double sigma_a_absorbing = parameters[3];
      sigma_s = std::shared_ptr(
          problem_->create_parameter_function(sigma_s_absorbing, sigma_s_scattering, sigma_s_scattering));
      sigma_a = std::shared_ptr(
          problem_->create_parameter_function(sigma_a_absorbing, sigma_a_scattering, sigma_a_scattering));
    }
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
    combined_operator_ =
        std::make_shared<CombinedOperatorType>(*density_operator_, *fv_operator_, *rhs_operator_, *hessian_inverter_);
    timestepper_ = std::make_shared<TimeStepperType>(
        *combined_operator_, *min_density_setter_, *analytical_flux_, *alpha_, true, 1., atol_, rtol_);
  }

  void create_combined_operator(const std::vector<double>& parameters)
  {
    return set_parameters(parameters);
  }

  VectorType get_initial_values() const
  {
    return *initial_values_alpha_;
  }

  bool finished() const
  {
    return XT::Common::FloatCmp::ge(timestepper_->current_time(), t_end_);
  }

  void visualize(const VectorType& alpha_vec, const std::string& prefix)
  {
    const ConstDiscreteFunctionType alpha(*fv_space_, alpha_vec, "alpha");
    const auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1,
        [&](const int /*comp*/, const auto& val) { return basis_functions_->density(analytical_flux_->get_u(val)); });
    alpha.visualize(*grid_view_, prefix, false, VTK::appendedraw, {}, *visualizer);
  }

  void init(const std::string output_dir,
            const size_t num_save_steps,
            const size_t grid_size,
            const bool visualize_solution,
            const bool silent,
            const std::vector<double>& parameters)
  {
    silent_ = silent;
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

    basis_functions_ =
        std::make_shared<const MomentBasis>(QuadratureChooser<MomentBasis, dimDomain == 1>::quad_order,
                                            QuadratureChooser<MomentBasis, dimDomain == 1>::quad_refinements);
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
    tmp_fulldimensional_vec_ =
        std::make_shared<VectorType>(fv_space_->mapper().size(), 0., DXTC_CONFIG_GET("num_alpha_mutexes", 1));
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
    initial_values_alpha_ = std::make_shared<VectorType>(*alpha_vec_);
    if (!silent_)
      std::cout << " done " << std::endl;

    const auto first_entity = *grid_view_->template begin<0>();
    const auto first_intersection = *grid_view_->ibegin(first_entity);
    dx_ = first_entity.geometry().volume() / first_intersection.geometry().volume();
    initial_dt_ = dx_ / 100.;
    t_end_ = DXTC_CONFIG_GET("t_end", problem_->t_end());

    // create operators
    numerical_flux_ = std::make_shared<KineticNumericalFluxType>(*analytical_flux_, *basis_functions_);
    advection_source_space_ = std::make_shared<AdvectionSourceSpaceType>(*grid_view_);
    advection_operator_ = std::make_shared<AdvectionOperatorType>(
        *grid_view_, *numerical_flux_, *advection_source_space_, *fv_space_, /*use_tbb*/ false);
    outside_indices_to_ignore_ = std::make_shared<std::map<size_t, std::set<size_t>>>();
    restricted_advection_operator_ = std::make_shared<AdvectionOperatorType>(*grid_view_,
                                                                             *numerical_flux_,
                                                                             *advection_source_space_,
                                                                             *fv_space_,
                                                                             /*use_tbb*/ false,
                                                                             XT::Grid::ApplyOn::NoIntersections<GV>(),
                                                                             true,
                                                                             *outside_indices_to_ignore_);
    reconstruction_operator_ = std::make_shared<ReconstructionOperatorType>(*fv_space_, *analytical_flux_);
    fv_operator_ = std::make_shared<ReconstructionAdvectionOperatorType>(
        *advection_operator_, *reconstruction_operator_, restricted_advection_operator_);
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
    create_combined_operator(parameters);
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
    restricted_advection_operator_->append(*boundary_lambda_, {}, boundary_intersection_filter);
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

  VectorType u_from_alpha(const VectorType& alpha_vec) const
  {
    DiscreteFunctionType u(*fv_space_, "u_initial");
    const ConstDiscreteFunctionType alpha(*fv_space_, alpha_vec, "alpha");
    // project initial values
    auto u_local_func = u.local_discrete_function();
    const auto alpha_local_func = alpha.local_discrete_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> alpha_local, u_local;
    for (auto&& element : Dune::elements(*grid_view_)) {
      u_local_func->bind(element);
      alpha_local_func->bind(element);
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_local[ii] = alpha_local_func->dofs().get_entry(ii);
      u_local = analytical_flux_->get_u(alpha_local);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_local_func->dofs().set_entry(ii, u_local[ii]);
    }
    return u.dofs().vector();
  }

  void
  maybe_insert_entity(const E& entity, const std::vector<size_t>& output_dofs, DynamicVector<size_t> global_dofs_entity)
  {
    // check if any output dof matches a dof on this entity
    fv_space_->mapper().global_indices(entity, global_dofs_entity);
    for (size_t ll = 0; ll < output_dofs.size(); ++ll) {
      for (size_t kk = 0; kk < global_dofs_entity.size(); ++kk) {
        if (global_dofs_entity[kk] == output_dofs[ll]) {
          restricted_op_output_entities_->insert(entity);
          restricted_op_input_entities_->insert(entity);
          for (auto&& intersection : intersections(*grid_view_, entity))
            if (intersection.neighbor())
              restricted_op_input_entities_->insert(intersection.outside());
          return;
        } // if (output dof found)
      } // kk
    } // ll
  }

  // Mostly copied from boltzmann.hh withou
  void prepare_restricted_operator(const std::vector<size_t>& output_dofs)
  {
    if (!restricted_op_output_dofs_ || *restricted_op_output_dofs_ != output_dofs) {
      restricted_op_output_dofs_ = std::make_shared<std::vector<size_t>>(output_dofs);
      restricted_op_input_dofs_ = std::make_shared<std::vector<size_t>>();
      const auto& mapper = fv_space_->mapper();
      DynamicVector<size_t> global_dofs;
      // calculate entities corresponding to dofs in restricted operator
      restricted_op_input_entities_ = std::make_shared<std::set<E, EntityLess<GV>>>(grid_view_->indexSet());
      restricted_op_output_entities_ = std::make_shared<std::set<E, EntityLess<GV>>>(grid_view_->indexSet());
      for (auto&& entity : elements(*grid_view_))
        maybe_insert_entity(entity, output_dofs, global_dofs);
      // find dofs corresponding to input entities
      for (const auto& entity : *restricted_op_input_entities_) {
        mapper.global_indices(entity, global_dofs);
        for (const auto& global_dof : global_dofs)
          restricted_op_input_dofs_->push_back(global_dof);
      }
      // if there are adjacent entities in the output_entities set, we have to make sure we do not apply the advection
      // operator to the intersection between these entities twice
      for (auto&& entity : *restricted_op_output_entities_) {
        const auto inside_index = grid_view_->indexSet().index(entity);
        (*outside_indices_to_ignore_)[inside_index] = std::set<size_t>();
        for (auto&& intersection : intersections(*grid_view_, entity)) {
          if (intersection.neighbor() && restricted_op_output_entities_->count(intersection.outside())) {
            const auto outside_index = grid_view_->indexSet().index(intersection.outside());
            if (inside_index > outside_index)
              (*outside_indices_to_ignore_)[inside_index].insert(outside_index);
          }
        } // intersections
      } // entities
    }
  } // ... prepare_restricted_operator(...)

  std::vector<size_t> restricted_op_input_dofs() const
  {
    return *restricted_op_input_dofs_;
  }

  size_t restricted_op_input_dofs_size() const
  {
    return restricted_op_input_dofs_->size();
  }

  double dx() const
  {
    return dx_;
  }

private:
  std::shared_ptr<const GridType> grid_;
  std::shared_ptr<const GV> grid_view_;
  std::shared_ptr<const MomentBasis> basis_functions_;
  std::shared_ptr<const ProblemType> problem_;
  std::shared_ptr<const SpaceType> fv_space_;
  std::shared_ptr<const AdvectionSourceSpaceType> advection_source_space_;
  std::shared_ptr<DiscreteFunctionType> u_;
  std::shared_ptr<VectorType> alpha_vec_;
  std::shared_ptr<VectorType> tmp_fulldimensional_vec_;
  std::shared_ptr<VectorType> initial_values_alpha_;
  std::shared_ptr<DiscreteFunctionType> alpha_;
  std::shared_ptr<KineticNumericalFluxType> numerical_flux_;
  std::shared_ptr<AdvectionOperatorType> advection_operator_;
  std::shared_ptr<AdvectionOperatorType> restricted_advection_operator_;
  std::shared_ptr<ReconstructionOperatorType> reconstruction_operator_;
  std::shared_ptr<ReconstructionAdvectionOperatorType> fv_operator_;
  std::shared_ptr<RhsOperatorType> rhs_operator_;
  std::shared_ptr<HessianInverterType> hessian_inverter_;
  std::shared_ptr<BoundaryDistributionType> boundary_distribution_;
  std::shared_ptr<DensityOperatorType> density_operator_;
  std::shared_ptr<MinDensitySetterType> min_density_setter_;
  std::shared_ptr<CombinedOperatorType> combined_operator_;
  std::shared_ptr<std::set<E, EntityLess<GV>>> restricted_op_input_entities_;
  std::shared_ptr<std::set<E, EntityLess<GV>>> restricted_op_output_entities_;
  std::shared_ptr<std::vector<size_t>> restricted_op_input_dofs_;
  std::shared_ptr<std::vector<size_t>> restricted_op_output_dofs_;
  std::shared_ptr<std::map<size_t, std::set<size_t>>> outside_indices_to_ignore_;
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