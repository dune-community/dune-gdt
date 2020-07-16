// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_EXAMPLES_BOLTZMANN_HH
#define DUNE_GDT_EXAMPLES_BOLTZMANN_HH

#include <sys/resource.h>

#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <random>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/timings.hh>

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/local/operators/generic.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/advection-fv-entropybased.hh>
#include <dune/gdt/operators/localizable-operator.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/planesource.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/sourcebeam.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/shadow.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/checkerboard.hh>
#include <dune/gdt/test/momentmodels/pn-discretization.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

using namespace Dune;
using namespace Dune::GDT;

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

template <class ProblemType,
          class GridType =
              YaspGrid<ProblemType::dimDomain, EquidistantOffsetCoordinates<double, ProblemType::dimDomain>>>
class BoltzmannSolver
{
public:
  // set dimensions
  // GDT::SphericalHarmonicsMomentBasis<double, double, momentOrder_or_refinements, 2, true>,
  // GDT::HatFunctionMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3>>;

  using MomentBasis = typename ProblemType::MomentBasis;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  static constexpr auto time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto rhs_time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step;
  using DomainFieldType = typename MomentBasis::DomainFieldType;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using RangeType = typename MomentBasis::RangeType;
  using GV = typename GridType::LeafGridView;
  using E = XT::Grid::extract_entity_t<GV>;
  static constexpr bool is_pn =
      std::is_same_v<
          ProblemType,
          GDT::CheckerboardPn<
              E,
              MomentBasis>> || std::is_same_v<ProblemType, GDT::PlaneSourcePn<E, MomentBasis>> || std::is_same_v<ProblemType, GDT::SourceBeamPn<E, MomentBasis>> || std::is_same_v<ProblemType, GDT::ShadowPn<E, MomentBasis>>;
  using SpaceType = FiniteVolumeSpace<GV, dimRange, 1, RangeFieldType>;
  using VectorType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::common_dense>::VectorType;
  using MatrixType = typename XT::LA::Container<RangeFieldType, XT::LA::Backends::common_dense>::MatrixType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GV, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GV, dimRange, 1, RangeFieldType>;
  using ConfigType = XT::Common::Configuration;
  using AnalyticalFluxType = typename ProblemType::FluxType;
  using EntropyFluxType = typename ProblemType::ActualFluxType;
  using InitialValueType = typename ProblemType::InitialValueType;
  using BoundaryValueType = typename ProblemType::BoundaryValueType;
  using GridProviderFactoryType = XT::Grid::CubeGridProviderFactory<GridType>;
  using GridProviderType = XT::Grid::GridProvider<GridType>;
  using AdvectionOperatorType = GDT::AdvectionFvOperator<MatrixType, GV, dimRange>;
  using EigenvectorWrapperType = typename EigenvectorWrapperChooser<MomentBasis, AnalyticalFluxType>::type;
  using EntropySolverType = GDT::EntropySolver<MomentBasis, SpaceType, MatrixType>;
  using FvOperatorType =
      std::conditional_t<is_pn,
                         AdvectionOperatorType,
                         GDT::EntropyBasedMomentFvOperator<AdvectionOperatorType, EntropySolverType>>;
  using RhsOperatorType = LocalizableOperator<MatrixType, GV, dimRange>;
  using FluxTimeStepperType =
      ExplicitRungeKuttaTimeStepper<FvOperatorType, DiscreteFunctionType, TimeStepperMethods::explicit_euler>;
  using KineticNumericalFluxType = NumericalKineticFlux<GV, MomentBasis>;
  using RhsTimeStepperType =
      ExplicitRungeKuttaTimeStepper<RhsOperatorType, DiscreteFunctionType, TimeStepperMethods::explicit_euler>;
  using TimeStepperType = FractionalTimeStepper<FluxTimeStepperType, RhsTimeStepperType>;
  using SolutionType = typename TimeStepperType::DiscreteSolutionType;
  using SolutionVectorsVectorType = std::vector<VectorType>;
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;

  static std::vector<double> default_parameters()
  {
    if constexpr (std::is_same_v<ProblemType, SourceBeamMn<GV, MomentBasis>>)
      return std::vector<double>{1, 0, 0, 2, 10};
    else
      return std::vector<double>{1, 0, 0, 10};
  }

  BoltzmannSolver(const std::string output_dir = "boltzmann",
                  const size_t num_save_steps = 10,
                  const size_t grid_size = 50,
                  const bool visualize_solution = true,
                  const bool silent = false,
                  const std::vector<double>& parameters = default_parameters(),
                  const double dt = -1.)
  {
    std::cout << "dt in is: " << dt << std::endl;
    auto num_save_steps_copy = num_save_steps;
    if (num_save_steps > 1e6) // hack to allow for size_t(-1) when called from the python bindings
      num_save_steps_copy = size_t(-1);
    init(output_dir, num_save_steps_copy, grid_size, visualize_solution, silent, parameters, dt);
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
    flux_timestepper_->current_solution().dofs().vector() = vec;
    rhs_timestepper_->current_solution().dofs().vector() = vec;
  }

  double time_step_length() const
  {
    return dt_;
  }

  bool linear() const
  {
    return is_pn;
  }

  SolutionVectorsVectorType solve(const bool with_half_steps = false)
  {
    if (!silent_)
      std::cout << "Solving... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    const auto visualizer = basis_functions_->visualizer();
    timestepper_->solve(t_end_,
                        dt_,
                        num_save_steps_,
                        silent_ ? size_t(0) : size_t(-1),
                        true,
                        false,
                        with_half_steps,
                        false,
                        false,
                        false,
                        file_path_,
                        *visualizer);
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    if (visualize_solution_) {
      if (!silent_)
        std::cout << "Visualizing... ";
      timestepper_->visualize_solution(file_path_, *visualizer);
      if (!silent_)
        std::cout << " done" << std::endl;
    }
    std::vector<VectorType> ret;
    for (const auto& pair : timestepper_->solution())
      ret.push_back(pair.second.dofs().vector());
    return ret;
  }

  SolutionVectorsVectorType next_n_timesteps(const size_t n, const bool with_half_steps = false) const
  {
    if (!silent_)
      std::cout << "Calculating next " << XT::Common::to_string(n) << " time steps... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    SolutionType solution;
    timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    std::vector<VectorType> ret;
    for (const auto& pair : solution)
      ret.push_back(pair.second.dofs().vector());
    return ret;
  }

  VectorType apply_kinetic_operator(const VectorType& source, const double time, const double dt) const
  {
    VectorType ret(source);
    kinetic_operator_->apply(source, ret, {{"t", {time}}, {"dt", {dt}}});
    return ret;
  }

  VectorType apply_restricted_kinetic_operator(const VectorType& source) const
  {
    if (is_pn)
      DUNE_THROW(NotImplemented, "This needs a Mn operator!");
    VectorType ret(restricted_op_output_dofs_->size(), 0.);
    RangeType u_entity;
    RangeType u_neighbor;
    const auto* mn_flux = dynamic_cast<const EntropyFluxType*>(flux_.get());

    double min_acceptable_density = 1e-9;
    size_t jj = 0;
    for (size_t ii = 0; ii < restricted_op_entities_->size(); ++ii) {
      const auto& entity = (*restricted_op_entities_)[ii];
      RangeType ret_entity(0.), local_ret(0.);
      for (size_t kk = 0; kk < dimRange; ++kk)
        u_entity[kk] = source[jj * dimRange + kk];
      ++jj;
      basis_functions_->ensure_min_density(u_entity, min_acceptable_density);
      for (auto&& intersection : intersections(*grid_view_, entity)) {
        const auto intersection_center = intersection.geometry().center();
        if (intersection.neighbor()) {
          for (size_t kk = 0; kk < dimRange; ++kk)
            u_neighbor[kk] = source[jj * dimRange + kk];
          basis_functions_->ensure_min_density(u_neighbor, min_acceptable_density);
          ++jj;
        } else if (intersection.boundary()) {
          u_neighbor = boundary_values_->evaluate(intersection_center);
        } else {
          DUNE_THROW(MathError, "This should not happen!");
        }
        assert(intersection.indexInInside() >= 0);
        size_t direction = static_cast<size_t>(intersection.indexInInside()) / 2;
        const auto local_intersection_center = intersection.geometry().local(intersection_center);
        auto n_ij = intersection.unitOuterNormal(local_intersection_center);
        const auto neighbor = intersection.neighbor() ? intersection.outside() : entity;
        mn_flux->evaluate_kinetic_flux(entity, neighbor, u_entity, u_neighbor, n_ij, direction, local_ret);
        ret_entity += local_ret * intersection.geometry().integrationElement(local_intersection_center);
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
  //    for (const auto& entity : elements(fv_space_->grid_view())) {
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

  VectorType apply_rhs_operator(VectorType source, const double time, const std::vector<double>& parameters)
  {
    set_parameters(parameters);
    return apply_rhs_operator(source, time);
  }

  void set_parameters(const std::vector<double>& parameters)
  {
    static const RangeType u_iso = basis_functions_->u_iso();
    static const RangeType basis_integrated = basis_functions_->integrated();
    std::shared_ptr<ParameterFunctionType> sigma_s;
    std::shared_ptr<ParameterFunctionType> sigma_a;
    if constexpr (std::is_same<ProblemType, SourceBeamMn<GV, MomentBasis>>::value) {
      DUNE_THROW_IF(parameters.size() != 5, Dune::InvalidStateException, "Wrong parameter size!");
      const double sigma_a_left = parameters[0];
      const double sigma_a_right = parameters[1];
      const double sigma_s_left = parameters[2];
      const double sigma_s_middle = parameters[3];
      const double sigma_s_right = parameters[4];
      sigma_s = std::shared_ptr(problem_->create_sigma_s_function(sigma_s_left, sigma_s_middle, sigma_s_right));
      sigma_a = std::shared_ptr(problem_->create_sigma_a_function(sigma_a_left, sigma_a_right));
    } else if constexpr (std::is_same<ProblemType, PlaneSourceMn<GV, MomentBasis>>::value) {
      DUNE_THROW_IF(parameters.size() != 4, Dune::InvalidStateException, "Wrong parameter size!");
      const double sigma_s_left = parameters[0];
      const double sigma_s_right = parameters[1];
      const double sigma_a_left = parameters[2];
      const double sigma_a_right = parameters[3];
      sigma_s = std::shared_ptr(problem_->create_parameter_function(sigma_s_left, sigma_s_right));
      sigma_a = std::shared_ptr(problem_->create_parameter_function(sigma_a_left, sigma_a_right));
    } else if constexpr (std::is_same<ProblemType, CheckerboardMn<GV, MomentBasis>>::value) {
      DUNE_THROW_IF(parameters.size() != 4, Dune::InvalidStateException, "Wrong parameter size!");
      const double sigma_s_scattering = parameters[0];
      const double sigma_s_absorbing = parameters[1];
      const double sigma_a_scattering = parameters[2];
      const double sigma_a_absorbing = parameters[3];
      sigma_s = std::shared_ptr(
          problem_->create_parameter_function(sigma_s_absorbing, sigma_s_scattering, sigma_s_scattering));
      sigma_a = std::shared_ptr(
          problem_->create_parameter_function(sigma_a_absorbing, sigma_a_scattering, sigma_a_scattering));
    } else {
      DUNE_THROW(Dune::NotImplemented, "");
    }

    std::shared_ptr<ParameterFunctionType> Q(problem_->Q());
    auto rhs_func = [this, sigma_a, sigma_s, Q](const auto& /*source*/,
                                                const auto& local_source,
                                                auto& local_range,
                                                const XT::Common::Parameter& /*param*/) {
      const auto& element = local_range.element();
      local_source[0]->bind(element);
      const auto center = element.geometry().center();
      const auto u = local_source[0]->evaluate(element.geometry().local(center));
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
    rhs_operator_ = std::make_shared<RhsOperatorType>(fv_space_->grid_view(), *fv_space_, *fv_space_);
    rhs_operator_->append(GenericLocalElementOperator<VectorType, GV, dimRange>(rhs_func, 1));
    rhs_timestepper_->set_operator(*rhs_operator_);
    // sigma_t_max_ = calculate_max_sigma_t(*sigma_s, *sigma_a);
    // dt_ =
    // std::min(problem_->CFL() * dx_, XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
  }

  VectorType get_initial_values() const
  {
    DiscreteFunctionType ret(*fv_space_, "discrete_initial_values");
    default_interpolation(*problem_->initial_values(), ret, *grid_view_);
    return ret.dofs().vector();
  }

  bool finished() const
  {
    const auto current_time = timestepper_->current_time();
    return XT::Common::FloatCmp::ge(current_time, t_end_);
  }

  void visualize(const VectorType& u_vec, const std::string& prefix)
  {
    const ConstDiscreteFunctionType u(*fv_space_, u_vec, "u");
    const auto rho_visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1, [&](const int /*comp*/, const auto& val) { return basis_functions_->density(val); });
    auto vtk_writer = u.create_vtkwriter(fv_space_->grid_view(), false);
    // visualize rho
    const auto rho_adapter =
        std::make_shared<XT::Functions::VisualizationAdapter<GV, dimRange, 1, double>>(u, *rho_visualizer, "rho");
    vtk_writer->addVertexData(rho_adapter);
    // visualize components of u
    std::vector<std::shared_ptr<XT::Functions::VisualizerInterface<dimRange, 1, double>>> comp_visualizers(dimRange);
    for (int ii = 0; ii < static_cast<int>(dimRange); ++ii) {
      comp_visualizers[ii] = std::make_shared<XT::Functions::ComponentVisualizer<dimRange, 1, double>>(ii);
      const auto adapter = std::make_shared<XT::Functions::VisualizationAdapter<GV, dimRange, 1, double>>(
          u, *comp_visualizers[ii], "u" + XT::Common::to_string(ii));
      vtk_writer->addVertexData(adapter);
    }
    u.write_visualization(*vtk_writer, prefix);
  }

  void init(const std::string output_dir,
            const size_t num_save_steps,
            const size_t grid_size,
            const bool visualize_solution,
            const bool silent,
            const std::vector<double>& parameters,
            const double dt)
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
    // setup threadmanager
    DXTC_CONFIG.set("global.datadir", output_dir, true);
    DXTC_TIMINGS.set_outputdir(output_dir);

    auto grid_config = ProblemType::default_grid_cfg();
    grid_config["num_elements"] = XT::Common::to_string(grid_size);

    basis_functions_ = std::make_shared<const MomentBasis>();

    // create grid
    auto grid_provider = GridProviderFactoryType::create(grid_config, MPIHelper::getCommunicator());
    grid_ = grid_provider.grid_ptr();

    // make a product finite volume space on the leaf grid
    grid_view_ = std::make_shared<const GV>(grid_->leafGridView());
    fv_space_ = std::make_shared<const SpaceType>(*grid_view_);

    static const RangeFieldType psi_vac = DXTC_CONFIG_GET("psi_vac", 1e-6 / basis_functions_->unit_ball_volume());
    if constexpr (is_pn) {
      problem_ = std::make_shared<const ProblemType>(*basis_functions_, psi_vac, grid_config);
    } else {
      problem_ = std::make_shared<const ProblemType>(*basis_functions_, *grid_view_, psi_vac, grid_config);
    }

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

    // calculate dx and choose t_end and dt
    const auto first_entity = *grid_view_->template begin<0>();
    const auto first_intersection = *grid_view_->ibegin(first_entity);
    dx_ = first_entity.geometry().volume() / first_intersection.geometry().volume();
    if (dt < 0.) {
      dt_ = std::min(problem_->CFL() * dx_,
                     XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
    } else {
      dt_ = dt;
    }
    t_end_ = problem_->t_end();

    // create Operators
    flux_ = std::shared_ptr<const AnalyticalFluxType>(problem_->flux());
    boundary_values_ = std::shared_ptr<const BoundaryValueType>(problem_->boundary_values());
    numerical_flux_ = std::make_shared<KineticNumericalFluxType>(*flux_, *basis_functions_);
    advection_operator_ =
        std::make_shared<AdvectionOperatorType>(*grid_view_, *numerical_flux_, *fv_space_, *fv_space_);
    entropy_solver_ =
        std::make_shared<EntropySolverType>(*(dynamic_cast<const EntropyFluxType*>(flux_.get())),
                                            *fv_space_,
                                            problem_->psi_vac() * basis_functions_->unit_ball_volume() / 10,
                                            file_path_);
    if constexpr (is_pn) {
      kinetic_operator_ = advection_operator_;
    } else {
      kinetic_operator_ = std::make_shared<const FvOperatorType>(*advection_operator_, *entropy_solver_);
    }

    // create timestepper
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
    // here, rhs_operator is still undefined, will be set in call to set_parameters
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(nullptr, *u_);
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
    set_parameters(parameters);
  } // void init()

  void reset()
  {
    u_ = std::make_shared<DiscreteFunctionType>(u_->space(), "solution");
    default_interpolation(*problem_->initial_values(), *u_, *grid_view_);
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
    FluxTimeStepperType::reset_static_variables();
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(*rhs_operator_, *u_);
    RhsTimeStepperType::reset_static_variables();
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
    TimeStepperType::reset_static_variables();
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
    if constexpr (is_pn)
      DUNE_THROW(NotImplemented, "This needs a Mn operator!");
    if (!restricted_op_output_dofs_ || *restricted_op_output_dofs_ != output_dofs) {
      restricted_op_output_dofs_ = std::make_shared<std::vector<size_t>>(output_dofs);
      restricted_op_input_dofs_ = std::make_shared<std::vector<size_t>>();
      restricted_op_entity_dofs_to_output_dofs_ = std::make_shared<std::vector<std::map<size_t, size_t>>>();
      const auto& mapper = fv_space_->mapper();
      DynamicVector<size_t> global_dofs_entity;
      DynamicVector<size_t> global_dofs_neighbor;
      // calculate entities corresponding to dofs in restricted operator
      restricted_op_entities_ = std::make_shared<std::vector<E>>();
      for (auto&& entity : elements(*grid_view_)) {
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
                for (auto&& intersection : intersections(*grid_view_, entity)) {
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
  std::shared_ptr<const GV> grid_view_;
  std::shared_ptr<const MomentBasis> basis_functions_;
  std::shared_ptr<const ProblemType> problem_;
  std::shared_ptr<const SpaceType> fv_space_;
  std::shared_ptr<DiscreteFunctionType> u_;
  std::shared_ptr<KineticNumericalFluxType> numerical_flux_;
  std::shared_ptr<const AdvectionOperatorType> advection_operator_;
  std::shared_ptr<const EntropySolverType> entropy_solver_;
  std::shared_ptr<const FvOperatorType> kinetic_operator_;
  std::shared_ptr<std::vector<E>> restricted_op_entities_;
  std::shared_ptr<std::vector<size_t>> restricted_op_input_dofs_;
  std::shared_ptr<std::vector<size_t>> restricted_op_output_dofs_;
  std::shared_ptr<std::vector<std::map<size_t, size_t>>> restricted_op_entity_dofs_to_output_dofs_;
  std::shared_ptr<RhsOperatorType> rhs_operator_;
  std::shared_ptr<FluxTimeStepperType> flux_timestepper_;
  std::shared_ptr<RhsTimeStepperType> rhs_timestepper_;
  std::shared_ptr<TimeStepperType> timestepper_;
  std::shared_ptr<const AnalyticalFluxType> flux_;
  std::shared_ptr<const BoundaryValueType> boundary_values_;
  double t_end_;
  double dt_;
  double dx_;
  bool silent_;
  bool visualize_solution_;
  std::string file_path_;
  size_t num_save_steps_;
  double sigma_t_max_;
};

#endif // DUNE_GDT_EXAMPLES_BOLTZMANN_HH
