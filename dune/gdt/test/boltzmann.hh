// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_BOLTZMANN_HH
#define DUNE_GDT_TEST_BOLTZMANN_HH

#include <sys/resource.h>

#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <random>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/generic.hh>
#include <dune/gdt/operators/advection-fv-entropybased.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>
#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>

#include <dune/gdt/test/momentmodels/kinetictransport/checkerboard.hh>
#include <dune/gdt/test/momentmodels/pn-discretization.hh>

using namespace Dune::GDT;
using namespace Dune::XT;

typedef std::mt19937 RandomNumberGeneratorType;

std::vector<double>
create_random_sigma_s_or_a(const double lower_bound, const double upper_bound, const size_t dimDomain)
{
  static RandomNumberGeneratorType rng{std::random_device()()};
  std::uniform_real_distribution<double> distribution(lower_bound, upper_bound);
  std::vector<double> ret;
  if (dimDomain == 2) {
    for (size_t row = 0; row < 7; ++row)
      for (size_t col = 0; col < 7; ++col)
        ret[7 * row + col] = distribution(rng);
  } else if (dimDomain == 3) {
    for (size_t plane = 0; plane < 7; ++plane)
      for (size_t row = 0; row < 7; ++row)
        for (size_t col = 0; col < 7; ++col) {
          ret[49 * plane + 7 * row + col] = distribution(rng);
        } // cols
  } // if (dimDomain == 2)
  return ret;
}

template <size_t dimDomain>
class BoltzmannSolver
{
public:
  // set dimensions
  static const size_t momentOrder_or_refinements = (dimDomain == 2) ? 15 : 2;

  using MomentBasis = typename std::conditional_t<
      dimDomain == 2,
      Dune::GDT::SphericalHarmonicsMomentBasis<double, double, momentOrder_or_refinements, 2, true>,
      //      Dune::GDT::PartialMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3, 1>>;
      Dune::GDT::HatFunctionMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3>>;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  static constexpr auto time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto rhs_time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step;
  using DomainFieldType = typename MomentBasis::DomainFieldType;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using RangeType = typename MomentBasis::RangeType;
  using GridType = Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<DomainFieldType, dimDomain>>;
  using GridViewType = typename GridType::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GridViewType>;
  using SpaceType = FiniteVolumeSpace<GridViewType, dimRange, 1, RangeFieldType>;
  using VectorType = typename Dune::XT::LA::Container<RangeFieldType, LA::Backends::common_dense>::VectorType;
  using MatrixType = typename Dune::XT::LA::Container<RangeFieldType, LA::Backends::common_dense>::MatrixType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using PnProblemType = Dune::GDT::CheckerboardPn<E, MomentBasis>;
  using ProblemType = Dune::GDT::CheckerboardMn<GridViewType, MomentBasis>;
  using ConfigType = Common::Configuration;
  using AnalyticalFluxType = typename ProblemType::FluxType;
  using EntropyFluxType = typename ProblemType::ActualFluxType;
  using InitialValueType = typename ProblemType::InitialValueType;
  using BoundaryValueType = typename ProblemType::BoundaryValueType;
  using GridProviderFactoryType = Dune::XT::Grid::CubeGridProviderFactory<GridType>;
  using GridProviderType = Dune::XT::Grid::GridProvider<GridType>;
  using AdvectionOperatorType = Dune::GDT::AdvectionFvOperator<MatrixType, GridViewType, dimRange>;
  using EigenvectorWrapperType = typename EigenvectorWrapperChooser<MomentBasis, AnalyticalFluxType>::type;
  using EntropySolverType = Dune::GDT::EntropySolver<MomentBasis, SpaceType, MatrixType>;
  using FvOperatorType = Dune::GDT::EntropyBasedMomentFvOperator<AdvectionOperatorType, EntropySolverType>;
  using PnFvOperatorType = AdvectionOperatorType;
  using RhsOperatorType = GenericOperator<MatrixType, GridViewType, dimRange>;
  using FluxTimeStepperType =
      ExplicitRungeKuttaTimeStepper<FvOperatorType, DiscreteFunctionType, TimeStepperMethods::explicit_euler>;
  using PnFluxTimeStepperType =
      ExplicitRungeKuttaTimeStepper<AdvectionOperatorType, DiscreteFunctionType, TimeStepperMethods::explicit_euler>;
  using KineticNumericalFluxType = NumericalKineticFlux<GridViewType, MomentBasis>;
  using RhsTimeStepperType =
      ExplicitRungeKuttaTimeStepper<RhsOperatorType, DiscreteFunctionType, TimeStepperMethods::explicit_euler>;
  using TimeStepperType = FractionalTimeStepper<FluxTimeStepperType, RhsTimeStepperType>;
  using PnTimeStepperType = FractionalTimeStepper<PnFluxTimeStepperType, RhsTimeStepperType>;
  using SolutionType = typename TimeStepperType::DiscreteSolutionType;
  using SolutionVectorsVectorType = std::vector<VectorType>;
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;

  BoltzmannSolver(const std::string output_dir = "boltzmann",
                  const size_t num_save_steps = 10,
                  const size_t grid_size = 50,
                  const bool visualize_solution = true,
                  const bool silent = false,
                  const RangeFieldType sigma_s_scattering = 1,
                  const RangeFieldType sigma_s_absorbing = 0,
                  const RangeFieldType sigma_a_scattering = 0,
                  const RangeFieldType sigma_a_absorbing = 10,
                  const bool linear = true)
  {
    auto num_save_steps_copy = num_save_steps;
    if (num_save_steps > 1e6) // hack to allow for size_t(-1) when called from the python bindings
      num_save_steps_copy = size_t(-1);
    init(output_dir,
         num_save_steps_copy,
         grid_size,
         visualize_solution,
         true,
         sigma_s_scattering,
         sigma_s_absorbing,
         sigma_a_scattering,
         sigma_a_absorbing,
         linear);
    silent_ = silent;
  }

  double current_time() const
  {
    return linear_ ? pn_timestepper_->current_time() : timestepper_->current_time();
  }

  double t_end() const
  {
    return linear_ ? pn_problem_->t_end() : problem_->t_end();
  }

  void set_current_time(const double time)
  {
    if (linear_) {
      pn_timestepper_->current_time() = time;
    } else {
      timestepper_->current_time() = time;
    }
  }

  void set_current_solution(const VectorType& vec)
  {
    if (linear_) {
      pn_timestepper_->current_solution().dofs().vector() = vec;
      pn_flux_timestepper_->current_solution().dofs().vector() = vec;
    } else {
      timestepper_->current_solution().dofs().vector() = vec;
      flux_timestepper_->current_solution().dofs().vector() = vec;
    }
    rhs_timestepper_->current_solution().dofs().vector() = vec;
  }

  double time_step_length() const
  {
    return dt_;
  }

  bool linear() const
  {
    return linear_;
  }

  SolutionVectorsVectorType solve(const bool with_half_steps = false)
  {
    if (!silent_)
      std::cout << "Solving... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    const auto visualizer = basis_functions_->visualizer();
    if (linear_) {
      pn_timestepper_->solve(t_end_,
                             dt_,
                             num_save_steps_,
                             silent_ ? size_t(0) : size_t(-1),
                             true,
                             false,
                             with_half_steps,
                             false,
                             false,
                             file_path_,
                             *visualizer);
    } else {
      timestepper_->solve(t_end_,
                          dt_,
                          num_save_steps_,
                          silent_ ? size_t(0) : size_t(-1),
                          true,
                          false,
                          with_half_steps,
                          false,
                          false,
                          file_path_,
                          *visualizer);
    }
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    if (visualize_solution_) {
      if (!silent_)
        std::cout << "Visualizing... ";
      linear_ ? pn_timestepper_->visualize_solution(file_path_, *visualizer)
              : timestepper_->visualize_solution(file_path_, *visualizer);
      if (!silent_)
        std::cout << " done" << std::endl;
    }
    std::vector<VectorType> ret;
    for (const auto& pair : (linear_ ? pn_timestepper_->solution() : timestepper_->solution()))
      ret.push_back(pair.second.dofs().vector());
    return ret;
  }

  SolutionVectorsVectorType next_n_time_steps(const size_t n, const bool with_half_steps = false) const
  {
    if (!silent_)
      std::cout << "Calculating next " << Common::to_string(n) << " time steps... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    SolutionType solution;
    if (linear_)
      pn_timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
    else
      timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    std::vector<VectorType> ret;
    for (const auto& pair : solution)
      ret.push_back(pair.second.dofs().vector());
    return ret;
  }

  VectorType apply_kinetic_operator(VectorType source, const double time, const double dt) const
  {
    VectorType ret(source);
    if (linear_) {
      pn_kinetic_operator_->apply(source, ret, {{"t", {time}}, {"dt", {dt}}});
    } else {
      kinetic_operator_->apply(source, ret, {{"t", {time}}, {"dt", {dt}}});
    }
    return ret;
  }

  VectorType apply_restricted_kinetic_operator(VectorType source) const
  {
    if (linear_)
      DUNE_THROW(Dune::NotImplemented, "This needs a Mn operator!");
    VectorType ret(restricted_op_output_dofs_->size(), 0.);
    RangeType u_entity;
    RangeType u_neighbor;
    const auto* mn_flux = dynamic_cast<const EntropyFluxType*>(flux_.get());

    double min_acceptable_density = 1e-9;
    size_t jj = 0;
    for (size_t ii = 0; ii < restricted_op_entities_->size(); ++ii) {
      const auto& entity = (*restricted_op_entities_)[ii];
      RangeType ret_entity(0.);
      for (size_t kk = 0; kk < dimRange; ++kk)
        u_entity[kk] = source[jj * dimRange + kk];
      ++jj;
      basis_functions_->ensure_min_density(u_entity, min_acceptable_density);
      for (auto&& intersection : Dune::intersections(*grid_view_, entity)) {
        const auto intersection_center = intersection.geometry().center();
        if (intersection.neighbor()) {
          for (size_t kk = 0; kk < dimRange; ++kk)
            u_neighbor[kk] = source[jj * dimRange + kk];
          basis_functions_->ensure_min_density(u_neighbor, min_acceptable_density);
          ++jj;
        } else if (intersection.boundary()) {
          u_neighbor = boundary_values_->evaluate(intersection_center);
        } else {
          DUNE_THROW(Dune::MathError, "This should not happen!");
        }
        assert(intersection.indexInInside() >= 0);
        size_t direction = static_cast<size_t>(intersection.indexInInside()) / 2;
        const auto local_intersection_center = intersection.geometry().local(intersection_center);
        auto n_ij = intersection.unitOuterNormal(local_intersection_center);
        const auto neighbor = intersection.neighbor() ? intersection.outside() : entity;
        ret_entity += mn_flux->evaluate_kinetic_flux(entity, neighbor, u_entity, u_neighbor, n_ij, direction)
                      * intersection.geometry().integrationElement(local_intersection_center);
      } // intersections
      ret_entity /= entity.geometry().volume();
      for (const auto& pair : (*restricted_op_entity_dofs_to_output_dofs_)[ii]) {
        assert(ret[pair.second] == 0. && "More than one output dofs maps to same vector location!");
        ret[pair.second] = ret_entity[pair.first];
      }
    } // entities
    return ret;
  }

  //  bool is_realizable(const VectorType& vector) const
  //  {
  //    const DiscreteFunctionType discrete_function(*fv_space_, vector);
  //    for (const auto& entity : Dune::elements(fv_space_->grid_view())) {
  //      const auto local_vec =
  //      discrete_function.local_function(entity)->evaluate(entity.geometry.local(entity.geometry().center())); if
  //      (!basis_functions_.is_realizable(local_vec))
  //        return false;
  //    }
  //    return true;
  //  }

  // for explicit euler timestepping
  VectorType apply_rhs_operator(VectorType source, const double time) const
  {
    const DiscreteFunctionType source_function(*fv_space_, source);
    VectorType ret(source);
    DiscreteFunctionType range_function(*fv_space_, ret);
    rhs_operator_->apply(source_function, range_function, time);
    return ret;
  }

  VectorType apply_rhs_operator(VectorType source,
                                const double time,
                                const RangeFieldType sigma_s_scattering,
                                const RangeFieldType sigma_s_absorbing = 0,
                                const RangeFieldType sigma_a_scattering = 1,
                                const RangeFieldType sigma_a_absorbing = 10)
  {
    set_rhs_operator_parameters(sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
    return apply_rhs_operator(source, time);
  }

  void set_rhs_operator_parameters(const RangeFieldType sigma_s_scattering = 1,
                                   const RangeFieldType sigma_s_absorbing = 0,
                                   const RangeFieldType sigma_a_scattering = 0,
                                   const RangeFieldType sigma_a_absorbing = 10)
  {
    const auto u_iso = basis_functions_->u_iso();
    const auto basis_integrated = basis_functions_->integrated();
    std::shared_ptr<ParameterFunctionType> sigma_s(
        problem_->create_parameter_function(sigma_s_absorbing, sigma_s_scattering, sigma_s_scattering));
    std::shared_ptr<ParameterFunctionType> sigma_a(
        problem_->create_parameter_function(sigma_a_absorbing, sigma_a_scattering, sigma_a_scattering));
    std::shared_ptr<ParameterFunctionType> Q(problem_->Q());
    auto rhs_func = [=](const auto& /*source*/,
                        const auto& local_source,
                        auto& local_range,
                        const Dune::XT::Common::Parameter& /*param*/) {
      const auto& element = local_range.element();
      local_source->bind(element);
      const auto center = element.geometry().center();
      const auto u = local_source->evaluate(center);
      const auto sigma_a_value = sigma_a->evaluate(center)[0];
      const auto sigma_s_value = sigma_s->evaluate(center)[0];
      const auto sigma_t_value = sigma_a_value + sigma_s_value;
      const auto Q_value = Q->evaluate(center)[0];
      auto ret = u * (-sigma_t_value);
      ret += u_iso * basis_functions_->density(u) * sigma_s_value;
      ret += basis_integrated * Q_value;
      for (size_t ii = 0; ii < local_range.dofs().size(); ++ii)
        local_range.dofs()[ii] = ret[ii];
    };
    rhs_operator_ = std::make_shared<const RhsOperatorType>(
        *fv_space_, *fv_space_, std::vector<typename RhsOperatorType::GenericElementFunctionType>(1, rhs_func));
    rhs_timestepper_->set_operator(*rhs_operator_);
    sigma_t_max_ = calculate_max_sigma_t(*sigma_s, *sigma_a);
    dt_ = std::min(problem_->CFL() * dx_,
                   Dune::XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
  }

  VectorType get_initial_values() const
  {
    DiscreteFunctionType ret(*fv_space_, "discrete_initial_values");
    const auto& problem = linear_ ? pn_problem_ : problem_;
    default_interpolation(*problem->initial_values(), ret, *grid_view_);
    return ret.dofs().vector();
  }

  bool finished() const
  {
    const auto current_time = linear_ ? pn_timestepper_->current_time() : timestepper_->current_time();
    return Common::FloatCmp::ge(current_time, t_end_);
  }

  void init(const std::string output_dir = "boltzmann",
            const size_t num_save_steps = 10,
            const size_t grid_size = 50,
            const bool visualize_solution = true,
            const bool silent = false,
            const double sigma_s_scattering = 1.,
            const double sigma_s_absorbing = 0.,
            const double sigma_a_scattering = 0.,
            const double sigma_a_absorbing = 10.,
            const bool linear = true)
  {
#if HAVE_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
      MPI_Init(NULL, NULL);
#endif
    silent_ = silent;
    if (!silent_)
      std::cout << "Setting problem parameters ...";
    visualize_solution_ = visualize_solution;
    file_path_ = output_dir;
    num_save_steps_ = num_save_steps;
    linear_ = linear;
    // setup threadmanager
    DXTC_CONFIG.set("global.datadir", output_dir, true);
    DXTC_TIMINGS.set_outputdir(output_dir);

    auto grid_config = ProblemType::default_grid_cfg();
    grid_config["num_elements"] = "[" + Common::to_string(grid_size);
    for (size_t ii = 1; ii < dimDomain; ++ii)
      grid_config["num_elements"] += " " + Common::to_string(grid_size);
    grid_config["num_elements"] += "]";

    basis_functions_ = std::make_shared<const MomentBasis>();

    // create grid
    auto grid_provider = GridProviderFactoryType::create(grid_config, Dune::MPIHelper::getCommunicator());
    grid_ = grid_provider.grid_ptr();

    // make a product finite volume space on the leaf grid
    grid_view_ = std::make_shared<const GridViewType>(grid_->leafGridView());
    fv_space_ = std::make_shared<const SpaceType>(*grid_view_);

    problem_ = std::make_shared<const ProblemType>(
        *basis_functions_, *grid_view_, grid_config, ProblemType::default_boundary_cfg());

    pn_problem_ =
        std::make_shared<const PnProblemType>(*basis_functions_, grid_config, ProblemType::default_boundary_cfg());

    // allocate a discrete function for the concentration
    u_ = std::make_shared<DiscreteFunctionType>(*fv_space_, "solution");

    // project initial values
    if (!silent_) {
      std::cout << " done " << std::endl;
      std::cout << "Projecting initial values...";
    }
    const auto initial_values = problem_->initial_values();
    default_interpolation(*initial_values, *u_, *grid_view_);
    if (!silent_)
      std::cout << " done " << std::endl;

    // calculate dx and choose t_end and initial dt
    Grid::Dimensions<GridViewType> dimensions(*grid_view_);
    if (dimDomain != 2 && dimDomain != 3)
      DUNE_THROW(Dune::NotImplemented, "Need to adjust dx calculation!");
    dx_ = dimensions.entity_width.max() / (dimDomain == 2 ? std::sqrt(2) : std::sqrt(3));
    dt_ = std::min(problem_->CFL() * dx_,
                   Dune::XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
    t_end_ = problem_->t_end();

    // create Operators
    flux_ = std::shared_ptr<const AnalyticalFluxType>(problem_->flux());
    pn_flux_ = std::shared_ptr<const AnalyticalFluxType>(pn_problem_->flux());
    boundary_values_ = std::shared_ptr<const BoundaryValueType>(problem_->boundary_values());
    numerical_flux_ = std::make_shared<KineticNumericalFluxType>(*flux_, *basis_functions_);
    pn_numerical_flux_ = std::make_shared<KineticNumericalFluxType>(*pn_flux_, *basis_functions_);
    advection_operator_ =
        std::make_shared<AdvectionOperatorType>(*grid_view_, *numerical_flux_, *fv_space_, *fv_space_);
    pn_kinetic_operator_ =
        std::make_shared<AdvectionOperatorType>(*grid_view_, *pn_numerical_flux_, *fv_space_, *fv_space_);
    entropy_solver_ =
        std::make_shared<EntropySolverType>(*(dynamic_cast<const EntropyFluxType*>(flux_.get())),
                                            *fv_space_,
                                            problem_->psi_vac() * basis_functions_->unit_ball_volume() / 10,
                                            file_path_);
    kinetic_operator_ = std::make_shared<const FvOperatorType>(*advection_operator_, *entropy_solver_);

    // create timestepper
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
    pn_flux_timestepper_ = std::make_shared<PnFluxTimeStepperType>(*pn_kinetic_operator_, *u_, -1.0);
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(*rhs_operator_, *u_);
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
    pn_timestepper_ = std::make_shared<PnTimeStepperType>(*pn_flux_timestepper_, *rhs_timestepper_);

    set_rhs_operator_parameters(sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
  } // void init()

  void reset()
  {
    u_ = std::make_shared<DiscreteFunctionType>(u_->space(), "solution");
    default_interpolation(*problem_->initial_values(), *u_, *grid_view_);
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
    pn_flux_timestepper_ = std::make_shared<PnFluxTimeStepperType>(*pn_kinetic_operator_, *u_, -1.0);
    FluxTimeStepperType::reset_static_variables();
    PnFluxTimeStepperType::reset_static_variables();
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(*rhs_operator_, *u_);
    RhsTimeStepperType::reset_static_variables();
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
    pn_timestepper_ = std::make_shared<PnTimeStepperType>(*pn_flux_timestepper_, *rhs_timestepper_);
    TimeStepperType::reset_static_variables();
    PnTimeStepperType::reset_static_variables();
  }

  double calculate_max_sigma_t(const ParameterFunctionType& /*sigma_s*/, const ParameterFunctionType& /*sigma_a*/) const
  {
    // temporarily disable as python side can't handle different time step length yet.
#if 0
    std::vector<double> sigma_t(sigma_s.size());
    std::transform(sigma_s.begin(), sigma_s.end(), sigma_a.begin(), sigma_t.begin(), std::plus<double>());
    return *std::max_element(sigma_t.begin(), sigma_t.end());
#endif
    return 16.; // for sigma_a, sigma_s in [0,8].
  }

  void prepare_restricted_operator(const std::vector<size_t>& output_dofs)
  {
    if (linear_)
      DUNE_THROW(Dune::NotImplemented, "This needs a Mn operator!");
    if (!restricted_op_output_dofs_ || *restricted_op_output_dofs_ != output_dofs) {
      restricted_op_output_dofs_ = std::make_shared<std::vector<size_t>>(output_dofs);
      restricted_op_input_dofs_ = std::make_shared<std::vector<size_t>>();
      restricted_op_entity_dofs_to_output_dofs_ = std::make_shared<std::vector<std::map<size_t, size_t>>>();
      const auto& mapper = fv_space_->mapper();
      Dune::DynamicVector<size_t> global_dofs_entity;
      Dune::DynamicVector<size_t> global_dofs_neighbor;
      // calculate entities corresponding to dofs in restricted operator
      restricted_op_entities_ = std::make_shared<std::vector<E>>();
      for (auto&& entity : Dune::elements(*grid_view_)) {
        // we only want to add each entity once, but there may be several output_dofs per entity
        bool entity_added = false;
        mapper.global_indices(entity, global_dofs_entity);
        // check if any output dof matches a dof on this entity
        for (size_t ll = 0; ll < output_dofs.size(); ++ll) {
          for (size_t kk = 0; kk < global_dofs_entity.size(); ++kk) {
            if (global_dofs_entity[kk] == output_dofs[ll]) {
              if (!entity_added) {
                restricted_op_entities_->push_back(entity);
                for (auto&& global_dof : global_dofs_entity)
                  restricted_op_input_dofs_->push_back(global_dof);
                for (auto&& intersection : Dune::intersections(*grid_view_, entity)) {
                  if (intersection.neighbor()) {
                    mapper.global_indices(intersection.outside(), global_dofs_neighbor);
                    for (auto&& global_dof : global_dofs_neighbor)
                      restricted_op_input_dofs_->push_back(global_dof);
                  }
                } // intersections
                restricted_op_entity_dofs_to_output_dofs_->emplace_back();
                entity_added = true;
              } // if (!entity_added)
              restricted_op_entity_dofs_to_output_dofs_->back().insert(std::make_pair(kk, ll));
            } // if (output dof found)
          } // kk
        } // ll
      } // entities
    }
  }

  std::vector<size_t> restricted_op_input_dofs() const
  {
    return *restricted_op_input_dofs_;
  }

  size_t restricted_op_input_dofs_size() const
  {
    return restricted_op_input_dofs_->size();
  }

private:
  std::shared_ptr<const GridType> grid_;
  std::shared_ptr<const GridViewType> grid_view_;
  std::shared_ptr<const MomentBasis> basis_functions_;
  std::shared_ptr<const ProblemType> problem_;
  std::shared_ptr<const PnProblemType> pn_problem_;
  std::shared_ptr<const SpaceType> fv_space_;
  std::shared_ptr<DiscreteFunctionType> u_;
  std::shared_ptr<KineticNumericalFluxType> numerical_flux_;
  std::shared_ptr<KineticNumericalFluxType> pn_numerical_flux_;
  std::shared_ptr<const AdvectionOperatorType> advection_operator_;
  std::shared_ptr<const EntropySolverType> entropy_solver_;
  std::shared_ptr<const FvOperatorType> kinetic_operator_;
  std::shared_ptr<const AdvectionOperatorType> pn_kinetic_operator_;
  std::shared_ptr<std::vector<E>> restricted_op_entities_;
  std::shared_ptr<std::vector<size_t>> restricted_op_input_dofs_;
  std::shared_ptr<std::vector<size_t>> restricted_op_output_dofs_;
  std::shared_ptr<std::vector<std::map<size_t, size_t>>> restricted_op_entity_dofs_to_output_dofs_;
  std::shared_ptr<const RhsOperatorType> rhs_operator_;
  std::shared_ptr<FluxTimeStepperType> flux_timestepper_;
  std::shared_ptr<PnFluxTimeStepperType> pn_flux_timestepper_;
  std::shared_ptr<RhsTimeStepperType> rhs_timestepper_;
  std::shared_ptr<TimeStepperType> timestepper_;
  std::shared_ptr<PnTimeStepperType> pn_timestepper_;
  std::shared_ptr<const AnalyticalFluxType> flux_;
  std::shared_ptr<const AnalyticalFluxType> pn_flux_;
  std::shared_ptr<const BoundaryValueType> boundary_values_;
  double t_end_;
  double dt_;
  double dx_;
  bool silent_;
  bool visualize_solution_;
  std::string file_path_;
  size_t num_save_steps_;
  size_t linear_;
  double sigma_t_max_;
};

#endif // DUNE_GDT_TEST_BOLTZMANN_HH