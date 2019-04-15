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
#include <dune/xt/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/checkerboard.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/spherical_harmonics.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/partial_moments.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/hatfunctions.hh>

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

#define MATEXP_ISO 0
#define EXPLICIT_EULER 1

template <size_t dimDomain>
class BoltzmannSolver
{
public:
  // set dimensions
  static const size_t momentOrder_or_refinements = (dimDomain == 2) ? 15 : 0;

  using BasisfunctionType = typename std::conditional_t<
      dimDomain == 2,
      Dune::GDT::SphericalHarmonicsMomentBasis<double, double, momentOrder_or_refinements, 2, true>,
      //      Dune::GDT::PartialMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3, 1>>;
      Dune::GDT::HatFunctionMomentBasis<double, 3, double, momentOrder_or_refinements, 1, 3>>;
  static constexpr size_t dimRange = BasisfunctionType::dimRange;
  static constexpr auto time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto rhs_time_stepper_method = TimeStepperMethods::explicit_euler;
  static constexpr auto time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step;
  using DomainFieldType = typename BasisfunctionType::DomainFieldType;
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  using RangeType = typename BasisfunctionType::RangeType;
  using GridType = Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<DomainFieldType, dimDomain>>;
  using GridLayerType = typename GridType::LeafGridView;
  using SpaceType = FvSpace<GridLayerType, RangeFieldType, dimRange, 1>;
  using VectorType = typename LA::Container<RangeFieldType, LA::Backends::common_dense>::VectorType;
  using DiscreteFunctionType = DiscreteFunction<SpaceType, VectorType>;
  //    using ProblemType = Dune::GDT::Hyperbolic::Problems::KineticTransport::
  //        CheckerboardPn<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
  using ProblemType = Dune::GDT::Hyperbolic::Problems::KineticTransport::
      CheckerboardMn<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
  using EquationType = Hyperbolic::Problems::KineticEquation<ProblemType>;
  using EntityType = typename GridType::template Codim<0>::Entity;

  using ConfigType = Common::Configuration;
  using AnalyticalFluxType = typename ProblemType::FluxType;
  using ActualAnalyticalFluxType = typename ProblemType::ActualFluxType;
  using RhsType = typename ProblemType::RhsType;
  using InitialValueType = typename ProblemType::InitialValueType;
  using BoundaryValueType = typename ProblemType::BoundaryValueType;
  using ActualBoundaryValueType = typename ProblemType::ActualBoundaryValueType;
  using GridProviderFactoryType = Grid::CubeGridProviderFactory<GridType>;
  using GridProviderType = Grid::GridProvider<GridType>;
  using MnFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, DiscreteFunctionType>;

  using ConstantFunctionType =
      Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using LfOperatorType =
      typename Dune::GDT::AdvectionLaxFriedrichsOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType>;
  using KineticOperatorType = typename Dune::GDT::
      AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType, BasisfunctionType, GridLayerType>;
  //  using FluxTimeStepperType =
  //      typename TimeStepperFactory<LfOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType;
  using FluxTimeStepperType =
      typename TimeStepperFactory<KineticOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType;
  using RhsOperatorType = typename Dune::GDT::AdvectionRhsOperator<RhsType>;
#if EXPLICIT_EULER
  using RhsTimeStepperType =
      typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, rhs_time_stepper_method>::TimeStepperType;
#elif MATEXP_ISO
  using RhsTimeStepperType = KineticIsotropicTimeStepper<DiscreteFunctionType, BasisfunctionType>;
#endif
  using TimeStepperType =
      typename Dune::GDT::TimeStepperSplittingFactory<FluxTimeStepperType,
                                                      RhsTimeStepperType,
                                                      time_stepper_splitting_method>::TimeStepperType;


  using SolutionType = typename TimeStepperType::DiscreteSolutionType;
  using SolutionVectorsVectorType = std::vector<VectorType>;

  BoltzmannSolver(const std::string output_dir,
                  const size_t num_save_steps,
                  const size_t grid_size,
                  const bool visualize_solution,
                  const bool silent,
                  const std::string sigma_s_in,
                  const std::string sigma_a_in)
  {
    auto num_save_steps_copy = num_save_steps;
    if (num_save_steps > 1e6) // hack to allow for size_t(-1) when called from the python bindings
      num_save_steps_copy = size_t(-1);
    auto random_sigma_s = create_random_sigma_s_or_a(0.0, 10.0, dimDomain);
    auto random_sigma_a = create_random_sigma_s_or_a(0.0, 10.0, dimDomain);
    auto sigma_s = sigma_s_in.empty() ? random_sigma_s : Common::from_string<std::vector<double>>(sigma_s_in);
    auto sigma_a = sigma_a_in.empty() ? random_sigma_a : Common::from_string<std::vector<double>>(sigma_a_in);
    init(output_dir, num_save_steps_copy, grid_size, visualize_solution, true, sigma_s, sigma_a);
    silent_ = silent;
  }

  BoltzmannSolver(const std::string output_dir = "boltzmann",
                  const size_t num_save_steps = 10,
                  const size_t grid_size = 50,
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

    std::vector<double> sigma_s(static_cast<size_t>(std::pow(7, dimDomain)));
    std::vector<double> sigma_a(sigma_s);
    create_sigma_s_and_a(
        sigma_s, sigma_a, sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
    init(output_dir, num_save_steps_copy, grid_size, visualize_solution, true, sigma_s, sigma_a);
    silent_ = silent;
  }

  static void create_sigma_s_and_a(std::vector<double>& sigma_s,
                                   std::vector<double>& sigma_a,
                                   const RangeFieldType sigma_s_scattering,
                                   const RangeFieldType sigma_s_absorbing,
                                   const RangeFieldType sigma_a_scattering,
                                   const RangeFieldType sigma_a_absorbing)
  {
    if (dimDomain == 2) {
      for (size_t row = 0; row < 7; ++row)
        for (size_t col = 0; col < 7; ++col) {
          sigma_s[7 * row + col] = ProblemType::is_absorbing(row, col) ? sigma_s_absorbing : sigma_s_scattering;
          sigma_a[7 * row + col] = ProblemType::is_absorbing(row, col) ? sigma_a_absorbing : sigma_a_scattering;
        } // cols
    } else if (dimDomain == 3) {
      for (size_t plane = 0; plane < 7; ++plane)
        for (size_t row = 0; row < 7; ++row)
          for (size_t col = 0; col < 7; ++col) {
            sigma_s[49 * plane + 7 * row + col] =
                ProblemType::is_absorbing(plane, row, col) ? sigma_s_absorbing : sigma_s_scattering;
            sigma_a[49 * plane + 7 * row + col] =
                ProblemType::is_absorbing(plane, row, col) ? sigma_a_absorbing : sigma_a_scattering;
          } // cols
    } // if (dimDomain == 2)
  }

  double current_time() const
  {
    return timestepper_->current_time();
  }

  double t_end() const
  {
    return 3.2;
  }

  void set_current_time(const double time)
  {
    timestepper_->current_time() = time;
  }

  void set_current_solution(const VectorType& vec)
  {
    timestepper_->current_solution().vector() = vec;
    rhs_timestepper_->current_solution().vector() = vec;
    flux_timestepper_->current_solution().vector() = vec;
  }

  double time_step_length() const
  {
    return dt_;
  }

  SolutionVectorsVectorType solve(const bool with_half_steps = false)
  {
    if (!silent_)
      std::cout << "Solving... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    const auto visualizer = basis_functions_->template visualizer<DiscreteFunctionType>();
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
                        visualizer);
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    if (visualize_solution_) {
      if (!silent_)
        std::cout << "Visualizing... ";
      timestepper_->visualize_solution(file_path_, visualizer);
      if (!silent_)
        std::cout << " done" << std::endl;
    }
    std::vector<VectorType> ret;
    for (const auto& pair : timestepper_->solution())
      ret.push_back(pair.second.vector());
    return ret;
  }

  SolutionVectorsVectorType next_n_time_steps(const size_t n, const bool with_half_steps = false) const
  {
    if (!silent_)
      std::cout << "Calculating next " << Common::to_string(n) << " time steps... " << std::endl;
    DXTC_TIMINGS.start("fv.solve");
    SolutionType solution;
    timestepper_->next_n_steps(n, t_end_, dt_, !silent_, with_half_steps, solution);
    DXTC_TIMINGS.stop("fv.solve");
    if (!silent_)
      std::cout << "Solving took: " << DXTC_TIMINGS.walltime("fv.solve") / 1000.0 << "s" << std::endl;
    std::vector<VectorType> ret;
    for (const auto& pair : solution)
      ret.push_back(pair.second.vector());
    return ret;
  }

  VectorType apply_LF_operator(VectorType source, const double time, const double dt) const
  {
    const DiscreteFunctionType source_function(*fv_space_, source);
    VectorType ret(source);
    DiscreteFunctionType range_function(*fv_space_, ret);
    lf_operator_->apply(source_function, range_function, {{"t", {time}}, {"dt", {dt}}});
    return ret;
  }

  VectorType apply_kinetic_operator(VectorType source, const double time, const double dt) const
  {
    const DiscreteFunctionType source_function(*fv_space_, source);
    VectorType ret(source);
    DiscreteFunctionType range_function(*fv_space_, ret);
    kinetic_operator_->apply(source_function, range_function, {{"t", {time}}, {"dt", {dt}}});
    return ret;
  }

  VectorType apply_restricted_kinetic_operator(VectorType source) const
  {
    VectorType ret(restricted_op_output_dofs_->size(), 0.);
    RangeType u_entity;
    RangeType u_neighbor;
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
        const auto x_local_entity = entity.geometry().local(intersection.geometry().center());
        if (intersection.neighbor()) {
          for (size_t kk = 0; kk < dimRange; ++kk)
            u_neighbor[kk] = source[jj * dimRange + kk];
          basis_functions_->ensure_min_density(u_neighbor, min_acceptable_density);
          ++jj;
        } else if (intersection.boundary()) {
          const auto local_boundary_vals = boundary_values_->local_function(entity);
          u_neighbor = local_boundary_vals->evaluate(intersection, x_local_entity, u_entity);
        } else {
          DUNE_THROW(Dune::MathError, "This should not happen!");
        }
        assert(intersection.indexInInside() >= 0);
        size_t direction = static_cast<size_t>(intersection.indexInInside()) / 2;
        const auto local_intersection_center = intersection.geometry().local(intersection.geometry().center());
        auto n_ij = intersection.unitOuterNormal(local_intersection_center);
        const auto neighbor = intersection.neighbor() ? intersection.outside() : entity;
        ret_entity += mn_flux_->evaluate_kinetic_flux(entity,
                                                      x_local_entity,
                                                      u_entity,
                                                      neighbor,
                                                      neighbor.geometry().local(intersection.geometry().center()),
                                                      u_neighbor,
                                                      n_ij,
                                                      direction,
                                                      {},
                                                      {})
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
//    for (const auto& entity : Dune::elements(fv_space_->grid_layer())) {
//      const auto local_vec =
//      discrete_function.local_function(entity)->evaluate(entity.geometry.local(entity.geometry().center())); if
//      (!basis_functions_.is_realizable(local_vec))
//        return false;
//    }
//    return true;
//  }

// for kinetic isotropic matrix exponential timestepper
#if MATEXP_ISO
  void set_rhs_timestepper_parameters(const RangeFieldType sigma_s_scattering = 1,
                                      const RangeFieldType sigma_s_absorbing = 0,
                                      const RangeFieldType sigma_a_scattering = 0,
                                      const RangeFieldType sigma_a_absorbing = 10)
  {
    std::vector<double> sigma_s(static_cast<size_t>(std::pow(7, dimDomain)));
    std::vector<double> sigma_a(sigma_s);
    create_sigma_s_and_a(
        sigma_s, sigma_a, sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
    auto param = ProblemType::default_param();
    param.set("sigma_s", sigma_s, true);
    param.set("sigma_a", sigma_a, true);
    problem_ = std::make_shared<const ProblemType>(
        *basis_functions_, *grid_view_, problem_->grid_config(), problem_->boundary_config(), param);
    equation_ = std::make_shared<const EquationType>(*problem_);
    rhs_timestepper_->set_params(problem_->get_sigma_a(), problem_->get_sigma_s(), problem_->get_Q());
  }

  VectorType apply_rhs_timestepper(VectorType source, const double time, const double dt) const
  {
    const DiscreteFunctionType source_function(*fv_space_, source);
    VectorType ret(source);
    DiscreteFunctionType range_function(*fv_space_, ret);
    rhs_timestepper_->apply(source_function, range_function, time, dt);
    return ret;
  }
#endif

  // for explicit euler timestepping
#if EXPLICIT_EULER
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
    std::vector<double> sigma_s(static_cast<size_t>(std::pow(7, dimDomain)));
    std::vector<double> sigma_a(sigma_s);
    create_sigma_s_and_a(
        sigma_s, sigma_a, sigma_s_scattering, sigma_s_absorbing, sigma_a_scattering, sigma_a_absorbing);
    auto param = ProblemType::default_param();
    param.set("sigma_s", sigma_s, true);
    param.set("sigma_a", sigma_a, true);
    problem_ = std::make_shared<const ProblemType>(
        *basis_functions_, *grid_view_, problem_->grid_config(), problem_->boundary_config(), param);
    equation_ = std::make_shared<const EquationType>(*problem_);
    rhs_operator_ = std::make_shared<const RhsOperatorType>(equation_->rhs());
    rhs_timestepper_->set_operator(*rhs_operator_);
    sigma_t_max_ = calculate_max_sigma_t(sigma_s, sigma_a);
    dt_ = std::min(equation_->CFL() * dx_,
                   Dune::XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
  }
#endif

  VectorType get_initial_values() const
  {
    DiscreteFunctionType ret(*fv_space_, "discrete_initial_values");
    project(equation_->initial_values(), ret);
    return ret.vector();
  }

  bool finished() const
  {
    return Common::FloatCmp::ge(timestepper_->current_time(), t_end_);
  }

  void init(const std::string output_dir = "boltzmann",
            const size_t num_save_steps = 10,
            const size_t grid_size = 50,
            const bool visualize_solution = true,
            const bool silent = false,
            const std::vector<double>& sigma_s = ProblemType::create_sigma_s(),
            const std::vector<double>& sigma_a = ProblemType::create_sigma_a())
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
    grid_config["num_elements"] = "[" + Common::to_string(grid_size);
    for (size_t ii = 1; ii < dimDomain; ++ii)
      grid_config["num_elements"] += " " + Common::to_string(grid_size);
    grid_config["num_elements"] += "]";

    auto param = ProblemType::default_param();
    param.set("sigma_s", sigma_s, true);
    param.set("sigma_a", sigma_a, true);
    basis_functions_ = std::make_shared<const BasisfunctionType>();

    sigma_t_max_ = calculate_max_sigma_t(sigma_s, sigma_a);

    // create grid
    auto grid_provider = GridProviderFactoryType::create(grid_config, Dune::MPIHelper::getCommunicator());
    grid_ = grid_provider.grid_ptr();

    // make a product finite volume space on the leaf grid
    grid_view_ = std::make_shared<const GridLayerType>(grid_->leafGridView());
    fv_space_ = std::make_shared<const SpaceType>(*grid_view_);

    problem_ = std::make_shared<const ProblemType>(
        *basis_functions_, *grid_view_, grid_config, ProblemType::default_boundary_cfg(), param);
    equation_ = std::make_shared<const EquationType>(*problem_);

    // allocate a discrete function for the concentration
    u_ = std::make_shared<DiscreteFunctionType>(*fv_space_, "solution");

    // project initial values
    if (!silent_) {
      std::cout << " done " << std::endl;
      std::cout << "Projecting initial values...";
    }
    project(equation_->initial_values(), *u_);
    if (!silent_)
      std::cout << " done " << std::endl;

    // calculate dx and choose t_end and initial dt
    Grid::Dimensions<GridLayerType> dimensions(*grid_view_);
    if (dimDomain != 2 && dimDomain != 3)
      DUNE_THROW(Dune::NotImplemented, "Need to adjust dx calculation!");
    dx_ = dimensions.entity_width.max() / (dimDomain == 2 ? std::sqrt(2) : std::sqrt(3));
    dt_ = std::min(equation_->CFL() * dx_,
                   Dune::XT::Common::FloatCmp::ne(sigma_t_max_, 0.) ? 0.99 * 1. / sigma_t_max_ : 1e100);
    t_end_ = equation_->t_end();

    // create Operators
    dx_function_ = std::make_shared<const ConstantFunctionType>(dx_);
    flux_ = std::shared_ptr<const AnalyticalFluxType>(problem_->create_flux());
    boundary_values_ = std::shared_ptr<const BoundaryValueType>(problem_->create_boundary_values());
    lf_operator_ = std::make_shared<const LfOperatorType>(*flux_, *boundary_values_, *dx_function_);
    kinetic_operator_ = std::make_shared<const KineticOperatorType>(*flux_, *boundary_values_, *basis_functions_);
    rhs_operator_ = std::make_shared<const RhsOperatorType>(equation_->rhs());

    // create timestepper
    //    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*lf_operator_, *u_, -1.0);
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
#if EXPLICIT_EULER
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(*rhs_operator_, *u_);
#elif MATEXP_ISO
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(
        *basis_functions_, *u_, problem_->get_sigma_a(), problem_->get_sigma_s(), problem_->get_Q());
#endif
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
  } // void init()

  void reset()
  {
    u_ = std::make_shared<DiscreteFunctionType>(u_->space(), "solution");
    project(equation_->initial_values(), *u_);
    //    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*lf_operator_, *u_, -1.0);
    flux_timestepper_ = std::make_shared<FluxTimeStepperType>(*kinetic_operator_, *u_, -1.0);
#if EXPLICIT_EULER
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(*rhs_operator_, *u_);
#elif MATEXP_ISO
    rhs_timestepper_ = std::make_shared<RhsTimeStepperType>(
        *basis_functions_, *u_, problem_->get_sigma_a(), problem_->get_sigma_s(), problem_->get_Q());
#endif
    timestepper_ = std::make_shared<TimeStepperType>(*flux_timestepper_, *rhs_timestepper_);
  }

  double calculate_max_sigma_t(const std::vector<double>& /*sigma_s*/, const std::vector<double>& /*sigma_a*/) const
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
    if (!restricted_op_output_dofs_ || *restricted_op_output_dofs_ != output_dofs) {
      mn_flux_ = std::make_shared<MnFluxType>(*basis_functions_, *grid_view_);
      restricted_op_output_dofs_ = std::make_shared<std::vector<size_t>>(output_dofs);
      restricted_op_input_dofs_ = std::make_shared<std::vector<size_t>>();
      restricted_op_entity_dofs_to_output_dofs_ = std::make_shared<std::vector<std::map<size_t, size_t>>>();
      const auto& mapper = fv_space_->mapper();
      Dune::DynamicVector<size_t> global_dofs_entity;
      Dune::DynamicVector<size_t> global_dofs_neighbor;
      // calculate entities corresponding to dofs in restricted operator
      restricted_op_entities_ = std::make_shared<std::vector<EntityType>>();
      for (auto&& entity : Dune::elements(*grid_view_)) {
        // we only want to add each entity once, but there may be several output_dofs per entity
        bool entity_added = false;
        mapper.globalIndices(entity, global_dofs_entity);
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
                    mapper.globalIndices(intersection.outside(), global_dofs_neighbor);
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
  std::shared_ptr<const GridLayerType> grid_view_;
  std::shared_ptr<const BasisfunctionType> basis_functions_;
  std::shared_ptr<const ProblemType> problem_;
  std::shared_ptr<const EquationType> equation_;
  std::shared_ptr<const SpaceType> fv_space_;
  std::shared_ptr<DiscreteFunctionType> u_;
  std::shared_ptr<const LfOperatorType> lf_operator_;
  std::shared_ptr<const KineticOperatorType> kinetic_operator_;
  std::shared_ptr<MnFluxType> mn_flux_;
  std::shared_ptr<std::vector<EntityType>> restricted_op_entities_;
  std::shared_ptr<std::vector<size_t>> restricted_op_input_dofs_;
  std::shared_ptr<std::vector<size_t>> restricted_op_output_dofs_;
  std::shared_ptr<std::vector<std::map<size_t, size_t>>> restricted_op_entity_dofs_to_output_dofs_;
  std::shared_ptr<const RhsOperatorType> rhs_operator_;
  std::shared_ptr<FluxTimeStepperType> flux_timestepper_;
  std::shared_ptr<RhsTimeStepperType> rhs_timestepper_;
  std::shared_ptr<TimeStepperType> timestepper_;
  std::shared_ptr<const AnalyticalFluxType> flux_;
  std::shared_ptr<const BoundaryValueType> boundary_values_;
  std::shared_ptr<const ConstantFunctionType> dx_function_;
  double t_end_;
  double dt_;
  double dx_;
  bool silent_;
  bool visualize_solution_;
  std::string file_path_;
  size_t num_save_steps_;
  double sigma_t_max_;
};

using BoltzmannSolver2d = BoltzmannSolver<2>;
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

    class_<Vec>(classname.c_str())
        .def(init<const size_t, const ScalarType, optional<const size_t>>())
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
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads2d, BoltzmannSolver2d::init, 0, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_time_steps_overloads2d, BoltzmannSolver2d::next_n_time_steps, 1, 2)
#if EXPLICIT_EULER
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads2d, BoltzmannSolver2d::apply_rhs_operator, 3, 6)
#endif

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads3d, BoltzmannSolver3d::init, 0, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_time_steps_overloads3d, BoltzmannSolver3d::next_n_time_steps, 1, 2)
#if EXPLICIT_EULER
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads3d, BoltzmannSolver3d::apply_rhs_operator, 3, 6)
#endif
#include <dune/xt/common/reenable_warnings.hh>


BOOST_PYTHON_MODULE(libboltzmann)
{
  typedef typename BoltzmannSolver2d::VectorType VectorType;
  typedef typename BoltzmannSolver2d::RangeFieldType RangeFieldType;
  // 2d
#if EXPLICIT_EULER
  VectorType (BoltzmannSolver2d::*apply_rhs_without_params2d)(VectorType, const double) const =
      &BoltzmannSolver2d::apply_rhs_operator;
  VectorType (BoltzmannSolver2d::*apply_rhs_with_params2d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver2d::apply_rhs_operator;
#endif

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
      .def(init<const std::string,
                const size_t,
                const size_t,
                const bool,
                const bool,
                const std::string,
                const std::string>())
      .def("init", &BoltzmannSolver2d::init, init_overloads2d())
      .def("solve", &BoltzmannSolver2d::solve)
      .def("next_n_time_steps", &BoltzmannSolver2d::next_n_time_steps, next_n_time_steps_overloads2d())
      .def("reset", &BoltzmannSolver2d::reset)
      .def("finished", &BoltzmannSolver2d::finished)
      .def("apply_LF_operator", &BoltzmannSolver2d::apply_LF_operator)
      .def("apply_kinetic_operator", &BoltzmannSolver2d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver2d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver2d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs_size)
#if EXPLICIT_EULER
      .def("apply_rhs_operator", apply_rhs_without_params2d)
      .def("apply_rhs_operator", apply_rhs_with_params2d, apply_rhs_overloads2d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver2d::set_rhs_operator_parameters)
#elif MATEXP_ISO
      .def("apply_rhs_timestepper", &BoltzmannSolver2d::apply_rhs_timestepper)
      .def("set_rhs_timestepper_parameters", &BoltzmannSolver2d::set_rhs_timestepper_parameters)
#endif
      .def("get_initial_values", &BoltzmannSolver2d::get_initial_values)
      .def("current_time", &BoltzmannSolver2d::current_time)
      .def("set_current_time", &BoltzmannSolver2d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver2d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver2d::time_step_length)
      .def("t_end", &BoltzmannSolver2d::t_end);

  class_<typename BoltzmannSolver2d::SolutionVectorsVectorType>("SolutionVectorsVectorType")
      .def(vector_indexing_suite<typename BoltzmannSolver2d::SolutionVectorsVectorType>())
      .def("size", &BoltzmannSolver2d::SolutionVectorsVectorType::size);

  // 3d
#if EXPLICIT_EULER
  VectorType (BoltzmannSolver3d::*apply_rhs_without_params3d)(VectorType, const double) const =
      &BoltzmannSolver3d::apply_rhs_operator;
  VectorType (BoltzmannSolver3d::*apply_rhs_with_params3d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver3d::apply_rhs_operator;
#endif

  class_<BoltzmannSolver3d>("BoltzmannSolver3d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double>>())
      .def(init<const std::string,
                const size_t,
                const size_t,
                const bool,
                const bool,
                const std::string,
                const std::string>())
      .def("init", &BoltzmannSolver3d::init, init_overloads3d())
      .def("solve", &BoltzmannSolver3d::solve)
      .def("next_n_time_steps", &BoltzmannSolver3d::next_n_time_steps, next_n_time_steps_overloads3d())
      .def("reset", &BoltzmannSolver3d::reset)
      .def("finished", &BoltzmannSolver3d::finished)
      .def("apply_LF_operator", &BoltzmannSolver3d::apply_LF_operator)
      .def("apply_kinetic_operator", &BoltzmannSolver3d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver3d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver3d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs_size)
#if EXPLICIT_EULER
      .def("apply_rhs_operator", apply_rhs_without_params3d)
      .def("apply_rhs_operator", apply_rhs_with_params3d, apply_rhs_overloads3d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver3d::set_rhs_operator_parameters)
#elif MATEXP_ISO
      .def("apply_rhs_timestepper", &BoltzmannSolver3d::apply_rhs_timestepper)
      .def("set_rhs_timestepper_parameters", &BoltzmannSolver3d::set_rhs_timestepper_parameters)
#endif
      .def("get_initial_values", &BoltzmannSolver3d::get_initial_values)
      .def("current_time", &BoltzmannSolver3d::current_time)
      .def("set_current_time", &BoltzmannSolver3d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver3d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver3d::time_step_length)
      .def("t_end", &BoltzmannSolver3d::t_end);


  VectorExporter<typename LA::CommonDenseVector<double>>::export_("CommonDenseVector");
  VectorExporter<typename LA::IstlDenseVector<double>>::export_("IstlDenseVector");

  iterable_converter().from_python<std::vector<double>>();
  iterable_converter().from_python<std::vector<size_t>>();

  std_vector_to_python_converter<double>();
  std_vector_to_python_converter<size_t>();
}


int main(int argc, char* argv[])
{
  static const size_t dimDomain = 2;

  try {
    // parse options
    if (argc == 1)
      std::cout << "The following options are available: " << argv[0]
                << " [-output_dir DIR -num_save_steps INT -gridsize INT "
                << "  -sigma_s_1 FLOAT -sigma_s_2 FLOAT -sigma_a_1 FLOAT -sigma_a_2 FLOAT"
                << " --no_visualization --silent --random_parameters --totally_random_parameters]" << std::endl;

    size_t num_save_steps = 10;
    size_t grid_size = 20;
    bool visualize = true;
    bool silent = false;
    bool random_parameters = false;
    bool totally_random_parameters = false;
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
        if (totally_random_parameters) {
          std::cerr << "Options --random_parameters and --totally-random_parameters are not compatible!" << std::endl;
          return 1;
        }
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
      } else if (std::string(argv[i]) == "--totally_random_parameters") {
        if (random_parameters) {
          std::cerr << "Options --random_parameters and --totally-random_parameters are not compatible!" << std::endl;
          return 1;
        }
        if (parameters_given) {
          std::cerr << "You specified a value for at least one parameter so you can't use --totally_random_parameters!"
                    << std::endl;
          return 1;
        }
        totally_random_parameters = true;
      } else if (std::string(argv[i]) == "-gridsize") {
        if (i + 1 < argc) {
          grid_size = Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-gridsize option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_1") {
        if (random_parameters || totally_random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters or --totally_random_parameters!" << std::endl;
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
        if (random_parameters || totally_random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters or --totally_random_parameters!" << std::endl;
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
        if (random_parameters || totally_random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters or --totally_random_parameters!" << std::endl;
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
        if (random_parameters || totally_random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters or --totally_random_parameters!" << std::endl;
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
    if (totally_random_parameters) {
      const auto sigma_s_matrix = create_random_sigma_s_or_a(sigma_s_lower, sigma_s_upper, dimDomain);
      const auto sigma_a_matrix = create_random_sigma_s_or_a(sigma_a_lower, sigma_a_upper, dimDomain);
      parameterfile << "Random parameters chosen on each square of the 7x7 checkerboard domain were: " << std::endl
                    << "sigma_s: " << Common::to_string(sigma_s_matrix) << std::endl
                    << "sigma_a: " << Common::to_string(sigma_a_matrix) << std::endl;
      solver = std::make_shared<BoltzmannSolver3d>(output_dir,
                                                   num_save_steps,
                                                   grid_size,
                                                   visualize,
                                                   silent,
                                                   Common::to_string(sigma_s_matrix),
                                                   Common::to_string(sigma_a_matrix));
    } else {
      parameterfile << "Domain was composed of two materials, parameters were: " << std::endl
                    << "First material: sigma_s = " + Common::to_string(sigma_s_1)
                           + ", sigma_a = " + Common::to_string(sigma_a_1)
                    << std::endl
                    << "Second material: sigma_s = " + Common::to_string(sigma_s_2)
                           + ", sigma_a = " + Common::to_string(sigma_a_2)
                    << std::endl;

      solver = std::make_shared<BoltzmannSolver3d>(
          output_dir, num_save_steps, grid_size, visualize, silent, sigma_s_1, sigma_s_2, sigma_a_1, sigma_a_2);
    }

    DXTC_TIMINGS.start("solve_all");
    std::vector<size_t> output_dofs{15, 3, 27, 4, 200, 533};
    solver->prepare_restricted_operator(output_dofs);
    using VectorType = typename LA::Container<double, LA::Backends::common_dense>::VectorType;
    auto initial_vals = solver->get_initial_values();
    RandomNumberGeneratorType rng{std::random_device()()};
    std::uniform_real_distribution<double> distribution(1, 100);
    for (auto&& val : initial_vals)
      val *= 1e2 * distribution(rng);
    auto source_dofs = solver->restricted_op_input_dofs();
    VectorType initial_vals_restr(source_dofs.size());
    for (size_t kk = 0; kk < source_dofs.size(); ++kk)
      initial_vals_restr[kk] = initial_vals[source_dofs[kk]];
    auto output1 = solver->apply_restricted_kinetic_operator(initial_vals_restr);
    auto output2 = solver->apply_kinetic_operator(initial_vals, 0, solver->time_step_length());
    for (size_t kk = 0; kk < output_dofs.size(); ++kk)
      std::cout << output1[kk] - output2[output_dofs[kk]] << std::endl;
    solver->solve();
    DXTC_TIMINGS.stop("solve_all");
    parameterfile << "Elapsed time: " << DXTC_TIMINGS.walltime("solve_all") / 1000.0 << " s" << std::endl;
    parameterfile.close();

    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported: " << e.what() << std::endl;
    std::abort();
  }
} // ... main(...)
