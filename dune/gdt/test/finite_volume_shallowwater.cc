// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <functional>
#include <vector>
#include <iostream>

#include "config.h"

#if HAVE_TBB
#include <tbb/task_scheduler_init.h>
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/string.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>

#include <dune/gdt/test/hyperbolic/problems/shallowwater.hh>

void trim(std::vector<std::string>& v)
{
  for (auto& s : v) {
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
    s = (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
  }
}


int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  MPIHelper::instance(argc, argv);

  // ***************** parse arguments and set up MPI and TBB
  size_t num_threads = 1;
  size_t threading_partition_factor = 1;
  size_t num_save_steps = 100;
  std::string grid_size("100"), overlap_size("2");
  double t_end = 0;
  double epsilon = 1e-10;
  std::string filename;
  double CFL = 0;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--num_threads") {
      if (i + 1 < argc) {
        num_threads = Dune::XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "--num_threads option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--CFL") {
      if (i + 1 < argc) {
        CFL = XT::Common::from_string<double>(argv[++i]);
      } else {
        std::cerr << "--CFL option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--threading_partition_factor") {
      if (i + 1 < argc) {
        threading_partition_factor = XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "--threading_partition_factor option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--filename") {
      if (i + 1 < argc) {
        filename = std::string(argv[++i]);
      } else {
        std::cerr << "--filename option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--num_save_steps") {
      if (i + 1 < argc) {
        num_save_steps = Dune::XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "--num_save_steps option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--grid_size") {
      if (i + 1 < argc) {
        grid_size = argv[++i];
      } else {
        std::cerr << "--grid_size option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--overlap_size") {
      if (i + 1 < argc) {
        overlap_size = argv[++i];
      } else {
        std::cerr << "--overlap_size option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--t_end") {
      if (i + 1 < argc) {
        t_end = XT::Common::from_string<double>(argv[++i]);
      } else {
        std::cerr << "--t_end option requires one argument." << std::endl;
        return 1;
      }
    } else {
      std::cerr << "Warning: Unrecognized option '" << argv[i] << "' ignored." << std::endl;
    }
  }
  (void)epsilon;

  DXTC_CONFIG.set("threading.partition_factor", threading_partition_factor, true);
  Dune::XT::Common::threadManager().set_max_threads(num_threads);
#if HAVE_TBB
  tbb::task_scheduler_init tbb_init(boost::numeric_cast<int>(num_threads));
#endif

  // ********************* choose dimensions ************************
  //  static const size_t dimDomain = 1;
  static const size_t dimDomain = 2;
  static const size_t dimRange = dimDomain == 1 ? 2 : 3;
  static const bool periodic_boundaries = false; // use periodic boundaries ?
  using DomainFieldType = double;
  using RangeFieldType = double;

  // ********************* choose timesteppers ************************
  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  const auto rhs_time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  const auto time_stepper_splitting_method = TimeStepperSplittingMethods::strang;

  // ********************* choose GridType ************************
  using GridType = Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<DomainFieldType, dimDomain>>;
  using GridLayerType = std::conditional_t<periodic_boundaries,
                                           typename XT::Grid::PeriodicGridView<GridType::LeafGridView, true>,
                                           typename GridType::LeafGridView>;
  using EntityType = typename GridLayerType::template Codim<0>::Entity;

  //******************** choose SpaceType and VectorType *****************************************
  using SpaceType = FvProductSpace<GridLayerType, RangeFieldType, dimRange, 1>;
  using VectorType = typename Dune::XT::LA::Container<double, Dune::XT::LA::default_backend>::VectorType;
  using DiscreteFunctionType = DiscreteFunction<SpaceType, VectorType>;

  //******************** choose ProblemType ***********************************************
  using ProblemType = Dune::GDT::Hyperbolic::Problems::ShallowWater<EntityType, DomainFieldType, DiscreteFunctionType>;

  //******************* get typedefs and constants from ProblemType **********************//
  using RhsType = typename ProblemType::RhsType;
  using InitialValueType = typename ProblemType::InitialValueType;
  using IntersectionType = typename GridLayerType::Intersection;

  //******************* create grid and FV space ***************************************
  auto grid_config = ProblemType::default_grid_cfg();
  grid_config["num_elements"] = grid_size;
  grid_config["overlap_size"] = overlap_size;
  const auto grid_ptr =
      Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
  assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
  const GridLayerType grid_layer(grid_ptr->leafGridView());
  const SpaceType fv_space(grid_layer);

  //******************* create ProblemType object ***************************************
  const ProblemType problem(grid_config);
  const InitialValueType& initial_values = problem.initial_values();
  const auto& dirichlet_boundary_values = problem.boundary_values();
  const auto boundary_info = XT::Grid::AllDirichletBoundaryInfo<IntersectionType>();
  using BoundaryValueType =
      LocalizableFunctionBasedLocalizableDirichletBoundaryValue<GridLayerType, typename ProblemType::BoundaryValueType>;
  const BoundaryValueType boundary_values(boundary_info, dirichlet_boundary_values);
  const RhsType& rhs = problem.rhs();

  // ***************** project initial values to discrete function *********************
  // create a discrete function for the solution
  DiscreteFunctionType u(fv_space, "solution");
  // project initial values
  project_l2(initial_values, u);

  // ************************* create analytical flux object ***************************************
  using AnalyticalFluxType = typename ProblemType::FluxType;
  const AnalyticalFluxType& analytical_flux = problem.flux();

  // ******************** choose flux and rhs operator and timestepper ******************************************
  typedef AdvectionRhsOperator<RhsType> RhsOperatorType;

  using AdvectionOperatorType = AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType>;

  using ReconstructionOperatorType =
      LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, SlopeLimiters::minmod>;
  //  using FvOperatorType = AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
  using FvOperatorType = AdvectionOperatorType;

  using OperatorTimeStepperType =
      typename TimeStepperFactory<FvOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType;
  using RhsOperatorTimeStepperType =
      typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, rhs_time_stepper_method>::TimeStepperType;
  using TimeStepperType = Dune::GDT::TimeStepperSplittingFactory<RhsOperatorTimeStepperType,
                                                                 OperatorTimeStepperType,
                                                                 time_stepper_splitting_method>::TimeStepperType;

  // *************** choose initial dt **************************************
  // calculate dx and choose initial dt
  Dune::XT::Grid::Dimensions<typename SpaceType::GridLayerType> dimensions(grid_layer);
  RangeFieldType dx = dimensions.entity_width.max();
  if (dimDomain == 2)
    dx /= std::sqrt(2);
  if (dimDomain == 3)
    dx /= std::sqrt(3);
  CFL = (CFL == 0.) ? problem.CFL() : CFL;
  RangeFieldType dt = CFL * dx;
  t_end = XT::Common::FloatCmp::eq(t_end, 0.) ? problem.t_end() : t_end;

  // *********************** create operators and timesteppers ************************************
  AdvectionOperatorType advection_operator(analytical_flux, boundary_values);
  RhsOperatorType rhs_operator(rhs);
  ReconstructionOperatorType reconstruction_operator(analytical_flux, boundary_values);
  //  FvOperatorType fv_operator(advection_operator, reconstruction_operator);
  FvOperatorType& fv_operator = advection_operator;

  // ******************************** do the time steps ***********************************************************
  OperatorTimeStepperType timestepper_op(fv_operator, u, -1.0);
  RhsOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
  TimeStepperType timestepper(timestepper_rhs, timestepper_op);
  //    TimeStepperType timestepper(timestepper_op, timestepper_rhs);
  //  TimeStepperType timestepper(fv_operator, u, -1.0);
  filename += "_" + ProblemType::static_id();
  filename += Dune::XT::Common::to_string(dimRange);

  timestepper.solve(t_end,
                    dt,
                    num_save_steps,
                    /* save_solution = */ false, // Save vector of calculated timesteps?
                    /* output_progress = */ true, // Progress written to std::cout?
                    /* visualize */ true, // vtp Output?
                    /* write_to_file */ true, // txt Output?
                    filename);

  return 0;
}
