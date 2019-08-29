// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

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

// using BoltzmannSolver2d = BoltzmannSolver<2>;
using BoltzmannSolver3d = BoltzmannSolver<3>;

// Python bindings
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

namespace {

// Converts a std::pair instance to a Python tuple.
template <typename T1, typename T2>
struct std_pair_to_tuple
{
  static PyObject* convert(std::pair<T1, T2> const& p)
  {
    return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr());
  }
  static PyTypeObject const* get_pytype()
  {
    return &PyTuple_Type;
  }
};

// Helper for convenience.
template <typename T1, typename T2>
struct std_pair_to_python_converter
{
  std_pair_to_python_converter()
  {
    boost::python::to_python_converter<std::pair<T1, T2>,
                                       std_pair_to_tuple<T1, T2>,
                                       true // std_pair_to_tuple has get_pytype
                                       >();
  }
};

// Converts a std::vector instance to a Python list.
template <typename T1>
struct std_vector_to_list
{
  static PyObject* convert(std::vector<T1> const& vec)
  {
    boost::python::list list;
    for (auto&& entry : vec) {
      list.append(entry);
    }
    return boost::python::incref(list.ptr());
  }

  static PyTypeObject const* get_pytype()
  {
    return &PyList_Type;
  }
};

template <typename T1>
struct std_vector_to_python_converter
{
  std_vector_to_python_converter()
  {
    boost::python::to_python_converter<std::vector<T1>, std_vector_to_list<T1>, true>();
  }
};


// The iterable_converter is copied from https://stackoverflow.com/a/15940413
// Converts any iterable type from python to C++
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter
{
  /// @note Registers converter from a python iterable type to the
  ///       provided type.
  template <typename Container>
  iterable_converter& from_python()
  {
    boost::python::converter::registry::push_back(&iterable_converter::convertible,
                                                  &iterable_converter::construct<Container>,
                                                  boost::python::type_id<Container>());

    // Support chaining.
    return *this;
  }

  /// @brief Check if PyObject is iterable.
  static void* convertible(PyObject* object)
  {
    return PyObject_GetIter(object) ? object : nullptr;
  }

  /// @brief Convert iterable PyObject to C++ container type.
  ///
  /// Container Concept requirements:
  ///
  ///   * Container::value_type is CopyConstructable.
  ///   * Container can be constructed and populated with two iterators.
  ///     I.e. Container(begin, end)
  template <typename Container>
  static void construct(PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    namespace python = boost::python;
    // Object is a borrowed reference, so create a handle indicting it is
    // borrowed for proper reference counting.
    python::handle<> handle(python::borrowed(object));

    // Obtain a handle to the memory block that the converter has allocated
    // for the C++ type.
    typedef python::converter::rvalue_from_python_storage<Container> storage_type;
    void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

    typedef python::stl_input_iterator<typename Container::value_type> iterator;

    // Allocate the C++ type into the converter's memory block, and assign
    // its handle to the converter's convertible variable.  The C++
    // container is populated by passing the begin and end iterators of
    // the python object to the container's constructor.
    new (storage) Container(iterator(python::object(handle)), // begin
                            iterator()); // end
    data->convertible = storage;
  }
};


} // namespace


template <class Vec>
struct VectorExporter
{
  typedef typename Vec::ScalarType ScalarType;
  typedef typename LA::VectorInterface<typename Vec::Traits, ScalarType> VectorInterfaceType;
  typedef typename VectorInterfaceType::derived_type derived_type;
  typedef typename Vec::RealType RealType;

  static object buffer(Vec& slf)
  {
    PyObject* py_buf =
        PyMemoryView_FromMemory((char*)&slf[0], slf.size() * sizeof(ScalarType), PyBUF_WRITE); // for python3
    //    PyObject* py_buf = PyBuffer_FromReadWriteMemory(&slf[0], slf.size() * sizeof(ScalarType)); // for python2
    object retval = object(handle<>(py_buf));
    return retval;
  }

  static std::shared_ptr<Vec> create_from_buffer(PyObject* memory_view, const size_t buffer_pos, const size_t vec_size)
  {
    Py_buffer* buffer = PyMemoryView_GET_BUFFER(memory_view);
    ScalarType* cxx_buf = (ScalarType*)buffer->buf;
    return std::make_shared<Vec>(vec_size, cxx_buf + buffer_pos, 0);
  }

  static void export_(const std::string& classname)
  {
    boost::python::type_info info = boost::python::type_id<std::pair<size_t, double>>();
    const boost::python::converter::registration* reg = boost::python::converter::registry::query(info);
    if (reg == nullptr) {
      std_pair_to_python_converter<size_t, double>();
    } else if ((*reg).m_to_python == nullptr) {
      std_pair_to_python_converter<size_t, double>();
    }

    void (Vec::*sub_void)(const derived_type&, derived_type&) const = &Vec::sub;
    derived_type (Vec::*sub_vec)(const derived_type&) const = &Vec::sub;

    void (Vec::*add_void)(const derived_type&, derived_type&) const = &Vec::add;
    derived_type (Vec::*add_vec)(const derived_type&) const = &Vec::add;

    class_<Vec, std::shared_ptr<Vec>>(classname.c_str())
        .def(init<const size_t, const ScalarType, optional<const size_t>>())
        .def("create_from_buffer", &create_from_buffer)
        .staticmethod("create_from_buffer")
        .def("size", &Vec::size)
        .def("add_to_entry", &Vec::add_to_entry)
        .def("__setitem__", &Vec::set_entry)
        .def("__getitem__", &Vec::get_entry)
        .def("l1_norm", &Vec::l1_norm)
        .def("l2_norm", &Vec::l2_norm)
        .def("sup_norm", &Vec::sup_norm)
        .def("standard_deviation", &Vec::standard_deviation)
        .def("set_all", &Vec::set_all)
        .def("valid", &Vec::valid)
        .def("dim", &Vec::size)
        .def("mean", &Vec::mean)
        .def("amax", &Vec::amax)
        .def("sub", sub_void)
        .def("sub", sub_vec)
        .def("add", add_void)
        .def("add", add_vec)
        .def("__add__", add_vec)
        .def("__sub__", sub_vec)
        .def("__iadd__", &Vec::iadd)
        .def("__isub__", &Vec::isub)
        .def("dot", &Vec::dot)
        .def("__mul__", &Vec::dot)
        .def("buffer", &buffer)
        .def("scal", &Vec::scal)
        .def("axpy", &Vec::axpy)
        .def("copy", &Vec::copy);
  }
};


#include <dune/xt/common/disable_warnings.hh>
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads2d, BoltzmannSolver2d::init, 0, 9)
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_time_steps_overloads2d, BoltzmannSolver2d::next_n_time_steps, 1, 2)
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads2d, BoltzmannSolver2d::apply_rhs_operator, 3, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads3d, BoltzmannSolver3d::init, 0, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_time_steps_overloads3d, BoltzmannSolver3d::next_n_time_steps, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads3d, BoltzmannSolver3d::apply_rhs_operator, 3, 6)
#include <dune/xt/common/reenable_warnings.hh>


BOOST_PYTHON_MODULE(libboltzmann)
{
  typedef typename BoltzmannSolver3d::VectorType VectorType;
  typedef typename BoltzmannSolver3d::RangeFieldType RangeFieldType;
#if 0
  // 2d
  VectorType (BoltzmannSolver2d::*apply_rhs_without_params2d)(VectorType, const double) const =
      &BoltzmannSolver2d::apply_rhs_operator;
  VectorType (BoltzmannSolver2d::*apply_rhs_with_params2d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver2d::apply_rhs_operator;

  class_<BoltzmannSolver2d>("BoltzmannSolver2d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double>>())
      .def("init", &BoltzmannSolver2d::init, init_overloads2d())
      .def("solve", &BoltzmannSolver2d::solve)
      .def("next_n_time_steps", &BoltzmannSolver2d::next_n_time_steps, next_n_time_steps_overloads2d())
      .def("reset", &BoltzmannSolver2d::reset)
      .def("finished", &BoltzmannSolver2d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver2d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver2d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver2d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params2d)
      .def("apply_rhs_operator", apply_rhs_with_params2d, apply_rhs_overloads2d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver2d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver2d::get_initial_values)
      .def("current_time", &BoltzmannSolver2d::current_time)
      .def("set_current_time", &BoltzmannSolver2d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver2d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver2d::time_step_length)
      .def("t_end", &BoltzmannSolver2d::t_end);

  class_<typename BoltzmannSolver2d::SolutionVectorsVectorType>("SolutionVectorsVectorType")
      .def(vector_indexing_suite<typename BoltzmannSolver2d::SolutionVectorsVectorType>())
      .def("size", &BoltzmannSolver2d::SolutionVectorsVectorType::size);
#endif

  // 3d
  VectorType (BoltzmannSolver3d::*apply_rhs_without_params3d)(VectorType, const double) const =
      &BoltzmannSolver3d::apply_rhs_operator;
  VectorType (BoltzmannSolver3d::*apply_rhs_with_params3d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver3d::apply_rhs_operator;

  class_<BoltzmannSolver3d>("BoltzmannSolver3d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double,
                                          const bool>>())
      .def("init", &BoltzmannSolver3d::init, init_overloads3d())
      .def("linear", &BoltzmannSolver3d::linear)
      .def("solve", &BoltzmannSolver3d::solve)
      .def("next_n_time_steps", &BoltzmannSolver3d::next_n_time_steps, next_n_time_steps_overloads3d())
      .def("reset", &BoltzmannSolver3d::reset)
      .def("finished", &BoltzmannSolver3d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver3d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver3d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver3d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params3d)
      .def("apply_rhs_operator", apply_rhs_with_params3d, apply_rhs_overloads3d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver3d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver3d::get_initial_values)
      .def("current_time", &BoltzmannSolver3d::current_time)
      .def("set_current_time", &BoltzmannSolver3d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver3d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver3d::time_step_length)
      .def("t_end", &BoltzmannSolver3d::t_end);

  VectorExporter<typename LA::CommonDenseVector<double>>::export_("CommonDenseVector");

  iterable_converter().from_python<std::vector<double>>();
  iterable_converter().from_python<std::vector<LA::CommonDenseVector<double>>>();
  iterable_converter().from_python<std::vector<size_t>>();

  std_vector_to_python_converter<double>();
  std_vector_to_python_converter<LA::CommonDenseVector<double>>();
  std_vector_to_python_converter<size_t>();
}


int main(int argc, char* argv[])
{
  try {
    // parse options
    if (argc == 1)
      std::cout << "The following options are available: " << argv[0]
                << " [-output_dir DIR -num_save_steps INT -gridsize INT "
                << "  -sigma_s_1 FLOAT -sigma_s_2 FLOAT -sigma_a_1 FLOAT -sigma_a_2 FLOAT"
                << " --no_visualization --silent --random_parameters]" << std::endl;

    size_t num_save_steps = 10;
    size_t grid_size = 20;
    bool visualize = true;
    bool silent = false;
    bool random_parameters = false;
    bool parameters_given = false;
    std::string output_dir;
    double sigma_s_lower = 0, sigma_s_upper = 8, sigma_a_lower = 0, sigma_a_upper = 8;
    double sigma_s_1 = 1, sigma_s_2 = 0, sigma_a_1 = 0, sigma_a_2 = 10;
    for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-output_dir") {
        if (i + 1 < argc) {
          output_dir = argv[++i];
        } else {
          std::cerr << "-output_dir option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-num_save_steps") {
        if (i + 1 < argc) {
          num_save_steps = Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-num_save_steps option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "--no_visualization") {
        visualize = false;
      } else if (std::string(argv[i]) == "--silent") {
        silent = true;
      } else if (std::string(argv[i]) == "--random_parameters") {
        if (parameters_given) {
          std::cerr << "You specified a value for at least one parameter so you can't use --random_parameters!"
                    << std::endl;
          return 1;
        }
        random_parameters = true;
        RandomNumberGeneratorType rng{std::random_device()()};
        std::uniform_real_distribution<double> sigma_s_dist(sigma_s_lower, sigma_s_upper);
        std::uniform_real_distribution<double> sigma_a_dist(sigma_a_lower, sigma_a_upper);
        sigma_s_1 = sigma_s_dist(rng);
        sigma_s_2 = sigma_s_dist(rng);
        sigma_a_1 = sigma_a_dist(rng);
        sigma_a_2 = sigma_a_dist(rng);
      } else if (std::string(argv[i]) == "-gridsize") {
        if (i + 1 < argc) {
          grid_size = Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-gridsize option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_1 = Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_2 = Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_2 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_1 = Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_2 = Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_2 option requires one argument." << std::endl;
          return 1;
        }
      } else {
        std::cerr << "Unknown option " << std::string(argv[i]) << std::endl;
        return 1;
      }
    }

    std::ofstream parameterfile;
    parameterfile.open(output_dir + "_parameters.txt");
    parameterfile << "Gridsize: " << Common::to_string(grid_size) + " x " + Common::to_string(grid_size) << std::endl;

    // run solver
    std::shared_ptr<BoltzmannSolver3d> solver;
    parameterfile << "Domain was composed of two materials, parameters were: " << std::endl
                  << "First material: sigma_s = " + Common::to_string(sigma_s_1)
                         + ", sigma_a = " + Common::to_string(sigma_a_1)
                  << std::endl
                  << "Second material: sigma_s = " + Common::to_string(sigma_s_2)
                         + ", sigma_a = " + Common::to_string(sigma_a_2)
                  << std::endl;

    solver = std::make_shared<BoltzmannSolver3d>(
        output_dir, num_save_steps, grid_size, visualize, silent, sigma_s_1, sigma_s_2, sigma_a_1, sigma_a_2);

    DXTC_TIMINGS.start("solve_all");
    std::vector<size_t> output_dofs{2728, 3868, 4468, 929};
    solver->prepare_restricted_operator(output_dofs);
    using VectorType = typename LA::Container<double, LA::Backends::common_dense>::VectorType;
    auto initial_vals = solver->get_initial_values();
    RandomNumberGeneratorType rng{std::random_device()()};
    std::uniform_real_distribution<double> distribution(1, 1000);
    for (auto&& val : initial_vals)
      val *= 1e4 * distribution(rng);
    auto source_dofs = solver->restricted_op_input_dofs();
    std::cout << source_dofs << std::endl;
#if 0
    VectorType initial_vals_restr{
        3.00887845e-05, 7.40090567e-05, 7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 1.00000000e-08,
        3.79301179e-05, 1.42264780e-05, 1.51590332e-05, 1.00000000e-08, 6.04617301e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 1.00000000e-08, 3.39443780e-04, 1.51590332e-05, 1.42264780e-05, 3.79301179e-05, 1.00000000e-08,
        1.00000000e-08, 6.04617301e-05, 1.47623908e-05, 9.55564780e-06, 9.55564780e-06, 1.47623908e-05, 1.00000000e-08,
        3.30163508e-05, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.20284294e-04, 1.00000000e-08, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08,
        1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.00000000e-08, 1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02,
        3.62665787e-02, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.50692562e-04, 2.63568632e-05, 6.74968960e-05, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 2.63568632e-05, 6.74968960e-05, 1.00000000e-08,
        1.50692562e-04, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05, 7.40090567e-05, 3.39443780e-04, 1.00000000e-08};
#endif
    VectorType initial_vals_restr(source_dofs.size());
    for (size_t kk = 0; kk < source_dofs.size(); ++kk)
      initial_vals_restr[kk] = initial_vals[source_dofs[kk]];
    auto output1 = solver->apply_restricted_kinetic_operator(initial_vals_restr);
    // for (size_t ii = 0; ii < 1000; ++ii) {
    const size_t ii = 0;
    const auto actual_initial_vals = initial_vals_restr * ii / 1000.;
    output1 = solver->apply_restricted_kinetic_operator(actual_initial_vals);
    // }
    auto output2 = solver->apply_kinetic_operator(initial_vals, 0, solver->time_step_length());
    for (size_t kk = 0; kk < output_dofs.size(); ++kk)
      EXPECT_NEAR(output1[kk], output2[output_dofs[kk]], 1e-10);
    const auto result = solver->solve();
    std::cout << " Result = " << std::accumulate(result.back().begin(), result.back().end(), 0.) << std::endl;
    DXTC_TIMINGS.stop("solve_all");
    parameterfile << "Elapsed time: " << DXTC_TIMINGS.walltime("solve_all") / 1000.0 << " s" << std::endl;
    parameterfile.close();

    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported: " << e.what() << std::endl;
    std::abort();
  }
} // ... main(...)
