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

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/istl.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshwriter.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
//#include <dune/gdt/test/hyperbolic/problems/momentmodels/fokkerplanck/checkerboard3d.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/fokkerplanck/twobeams.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/fokkerplanck/sourcebeam.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/twobeams.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/sourcebeam.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/planesource.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/pointsource.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/linesource.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/lebedevquadrature.hh>
#include <dune/gdt/test/hyperbolic/problems/transport.hh>

//! struct to be used as comparison function e.g. in a std::map<FieldVector<...>, ..., FieldVectorLess>
struct CmpStruct
{
  template <class FieldType, int dimDomain>
  bool operator()(const std::pair<Dune::FieldVector<FieldType, dimDomain>, FieldType>& a,
                  const std::pair<Dune::FieldVector<FieldType, dimDomain>, FieldType>& b) const
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      if (Dune::XT::Common::FloatCmp::lt(a.first[dd], b.first[dd]))
        return true;
      else if (Dune::XT::Common::FloatCmp::gt(a.first[dd], b.first[dd]))
        return false;
    }
    return false;
  }
};

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
  Eigen::initParallel();
  using namespace Dune;
  using namespace Dune::GDT;

  // ***************** parse arguments and set up MPI and TBB
  size_t num_threads = 1;
  size_t num_save_steps = -1;
  std::string grid_size("100"), overlap_size("1");
  double t_end = 0;
  double epsilon = 1e-10;
  //  double rel_tol = 1e-2;
  //  double abs_tol = 1e-2;
  bool visualize = true;
  std::string filename;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-num_threads") {
      if (i + 1 < argc) {
        num_threads = Dune::XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "-num_threads option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-filename") {
      if (i + 1 < argc) {
        filename = std::string(argv[++i]);
      } else {
        std::cerr << "-num_save_steps option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-epsilon") {
      if (i + 1 < argc) {
        epsilon = Dune::XT::Common::from_string<double>(argv[++i]);
      } else {
        std::cerr << "-num_save_steps option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-num_save_steps") {
      if (i + 1 < argc) {
        num_save_steps = Dune::XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "-num_save_steps option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-grid_size") {
      if (i + 1 < argc) {
        grid_size = argv[++i];
      } else {
        std::cerr << "-grid_size option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-overlap_size") {
      if (i + 1 < argc) {
        overlap_size = argv[++i];
      } else {
        std::cerr << "-overlap_size option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "-t_end") {
      if (i + 1 < argc) {
        t_end = XT::Common::from_string<double>(argv[++i]);
      } else {
        std::cerr << "-t_end option requires one argument." << std::endl;
        return 1;
      }
      //    } else if (std::string(argv[i]) == "-quadrature_rel_tol") {
      //      if (i + 1 < argc) {
      //        rel_tol = XT::Common::from_string<double>(argv[++i]);
      //      } else {
      //        std::cerr << "-quadrature_rel_tol option requires one argument." << std::endl;
      //        return 1;
      //      }
      //    } else if (std::string(argv[i]) == "-quadrature_abs_tol") {
      //      if (i + 1 < argc) {
      //        abs_tol = XT::Common::from_string<double>(argv[++i]);
      //      } else {
      //        std::cerr << "-quadrature_abs_tol option requires one argument." << std::endl;
      //        return 1;
      //      }
    } else if (std::string(argv[i]) == "--no_visualization") {
      visualize = false;
    } else {
      std::cerr << "Unknown option " << std::string(argv[i]) << std::endl;
      return 1;
    }
  }

  DXTC_CONFIG.set("threading.partition_factor", 1u, true);
  Dune::XT::Common::threadManager().set_max_threads(num_threads);

#define USE_SMP_PARALLEL 1
#if HAVE_DUNE_FEM
  Dune::Fem::MPIManager::initialize(argc, argv);
#else
  MPIHelper::instance(argc, argv);
#endif

  // ********************* choose dimensions, fluxes and grid type ************************
  static const int dimDomain = 3;
  //  static const int dimDomain = 1;
  static const int momentOrder = 6;
  //  const auto numerical_flux = NumericalFluxes::kinetic;
  //  const auto numerical_flux = NumericalFluxes::godunov;
  //  const auto numerical_flux = NumericalFluxes::laxfriedrichs;
  //  const auto numerical_flux = NumericalFluxes::laxfriedrichs_with_reconstruction;
  //  const auto numerical_flux = NumericalFluxes::local_laxfriedrichs_with_reconstruction;
  //      const auto numerical_flux = NumericalFluxes::local_laxfriedrichs;
  //  const auto time_stepper_method = TimeStepperMethods::explicit_euler;
  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  //  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_third_order_ssp;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::explicit_euler;
  //      const auto rhs_time_stepper_method = TimeStepperMethods::implicit_euler;
  const auto rhs_time_stepper_method = TimeStepperMethods::matrix_exponential;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::trapezoidal_rule;
  const auto time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step;

  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  //  typedef typename XT::Grid::PeriodicGridView<GridType::LeafGridView, true> GridLayerType;
  typedef typename GridType::LeafGridView GridLayerType;
  typedef typename GridType::Codim<0>::Entity EntityType;

  //******************** choose BasisfunctionType *****************************************
  //  typedef typename Hyperbolic::Problems::LegendrePolynomials<double, dimDomain, double, momentOrder>
  //  BasisfunctionType;

  //  static const size_t refinements = 0;
  //  typedef
  //      typename Hyperbolic::Problems::HatFunctions<double,
  //                                                  3,
  //                                                  double,
  //                                                  Hyperbolic::Problems::OctaederStatistics<refinements>::num_vertices(),
  //                                                  1,
  //                                                  dimDomain>
  //          BasisfunctionType;

  static const size_t refinements = 0;
  typedef typename Hyperbolic::Problems::
      PiecewiseMonomials<double,
                         3,
                         double,
                         4 * Hyperbolic::Problems::OctaederStatistics<refinements>::num_faces(),
                         1,
                         dimDomain>
          BasisfunctionType;

  //  typedef typename Hyperbolic::Problems::RealSphericalHarmonics<double, double, momentOrder, dimDomain, false>
  //      BasisfunctionType;

  //  std::shared_ptr<const BasisfunctionType> basis_functions = std::make_shared<const BasisfunctionType>();
  std::shared_ptr<const BasisfunctionType> basis_functions = std::make_shared<const BasisfunctionType>(refinements, 4);
  static const size_t dimRange = BasisfunctionType::dimRange;
  static const size_t dimFlux = BasisfunctionType::dimFlux;
  //  static const size_t dimRange = 1.;
  static constexpr auto container_backend = Dune::XT::LA::default_sparse_backend;
  typedef FvProductSpace<GridLayerType, double, dimRange, 1> SpaceType;
  typedef typename Dune::XT::LA::Container<double, container_backend>::VectorType VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;

  //******************** choose ProblemType ***********************************************
  //  typedef typename Hyperbolic::Problems::FokkerPlanck
  //      TwoBeamsPn<BasisfunctionType, GridLayerType, EntityType, double, dimDomain, DiscreteFunctionType,
  //      double,
  //      dimRange>
  //          ProblemImp;

  //    typedef typename Hyperbolic::Problems::KineticTransport
  //      TwoBeamsMn<BasisfunctionType, GridLayerType, EntityType, double, dimDomain, DiscreteFunctionType, double,
  //      dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport
  //      SourceBeamPn<BasisfunctionType, GridLayerType, EntityType, double, dimDomain, DiscreteFunctionType, double,
  //      dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::SourceBeamMn<
  //                                                      BasisfunctionType,
  //  GridLayerType,
  //                                                      EntityType,
  //                                                      double,
  //                                                      dimDomain,
  //                                                      DiscreteFunctionType,
  //                                                      double,
  //                                                      dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::PlaneSourcePn<BasisfunctionType,
  //                                                                         GridLayerType,
  //                                                                         EntityType,
  //                                                                         double,
  //                                                                         dimDomain,
  //                                                                         DiscreteFunctionType,
  //                                                                         double,
  //                                                                         dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::PlaneSourceMn<BasisfunctionType,
  //                                                                         GridLayerType,
  //                                                                         EntityType,
  //                                                                         double,
  //                                                                         dimDomain,
  //                                                                         DiscreteFunctionType,
  //                                                                         double,
  //                                                                         dimRange>
  //      ProblemImp;

  typedef typename Hyperbolic::Problems::KineticTransport::PointSourcePn<BasisfunctionType,
                                                                         GridLayerType,
                                                                         EntityType,
                                                                         double,
                                                                         dimDomain,
                                                                         DiscreteFunctionType,
                                                                         double,
                                                                         dimRange>
      ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::PointSourceMn<
  //                                                       BasisfunctionType,
  //  GridLayerType,
  //                                                       EntityType,
  //                                                       double,
  //                                                       dimDomain,
  //                                                       DiscreteFunctionType,
  //                                                       double,
  //                                                       dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::LineSourcePn<BasisfunctionType,
  //                                                                                GridLayerType,
  //                                                                                EntityType,
  //                                                                                double,
  //                                                                                dimDomain,
  //                                                                                DiscreteFunctionType,
  //                                                                                double,
  //                                                                                dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::KineticTransport::LineSourceMn<BasisfunctionType,
  //                                                                        GridLayerType,
  //                                                                        EntityType,
  //                                                                        double,
  //                                                                        dimDomain,
  //                                                                        DiscreteFunctionType,
  //                                                                        double,
  //                                                                        dimRange>
  //      ProblemImp;

  //  typedef
  //      typename Hyperbolic::Problems::Transport<EntityType, double, dimDomain, DiscreteFunctionType, double,
  //      dimRange>
  //          ProblemType;

  //******************* get typedefs and constants from ProblemType **********************//
  typedef typename Hyperbolic::Problems::KineticEquation<ProblemImp> ProblemType;
  using DomainFieldType = typename ProblemType::DomainFieldType;
  using DomainType = typename ProblemType::DomainType;
  using RangeFieldType = typename ProblemType::RangeFieldType;
  using RangeType = typename ProblemType::RangeType;
  typedef typename ProblemType::RhsType RhsType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::BoundaryValueType BoundaryValueType;
  static const bool linear = false; // ProblemType::linear;

  //******************* create grid and FV space ***************************************
  auto grid_config = ProblemType::default_grid_cfg();
  grid_config["num_elements"] = grid_size;
  grid_config["overlap_size"] = overlap_size;
  const auto grid_ptr = Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config).grid_ptr();
  const auto& grid = *grid_ptr;
  assert(grid.comm().size() == 1 || grid.overlapSize(0) > 0);
  const GridLayerType grid_layer(grid_ptr->leafGridView());
  const SpaceType fv_space(grid_layer);

  const auto quadrature = Hyperbolic::Problems::LebedevQuadrature<DomainFieldType, true>::get(80);
  //    const auto& quadrature = basis_functions->quadrature();
  //  const auto quadrature = ProblemImp::default_quadrature(grid_config);

  //******************* create ProblemType object ***************************************
  const std::unique_ptr<ProblemImp> problem_imp =
      XT::Common::make_unique<ProblemImp>(*basis_functions, grid_layer, quadrature, grid_config);
  //    const ProblemImp problem_imp(basis_functions, grid_layer, grid_config);
  //  const std::unique_ptr<ProblemImp> problem_imp =
  //      XT::Common::make_unique<ProblemImp>(*basis_functions, grid_layer, basis_functions->quadrature(), grid_config);
  //  const ProblemImp problem_imp(
  //      basis_functions, Hyperbolic::Problems::LebedevQuadrature<DomainFieldType, true>::get(1000), grid_layer);

  //  const ProblemImp problem_imp(*basis_functions,
  //                               grid_layer,
  //                               grid_config,
  //                               ProblemImp::default_boundary_cfg(),
  //                               //                               basis_functions.quadrature());

  const ProblemType problem(*problem_imp);
  //  const ProblemType problem(grid_config);
  const InitialValueType& initial_values = problem.initial_values();
  const BoundaryValueType& boundary_values = problem.boundary_values();
  const RhsType& rhs = problem.rhs();
  const RangeFieldType CFL = problem.CFL() * 0.9;

  // ***************** project initial values to discrete function *********************
  // create a discrete function for the solution
  DiscreteFunctionType u(fv_space, "solution");
  // project initial values
  project_l2(initial_values, u);

  // ************************* create analytical flux object ***************************************
  typedef typename ProblemType::FluxType AnalyticalFluxType;
  const AnalyticalFluxType& analytical_flux = problem.flux();

  // ******************** choose Realizability limiter and eigensolver ******************************************

  //  typedef ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType, 3, dimRange>
  //      RealizabilityLimiterType;
  //  typedef LPLocalRealizabilityLimiter<DiscreteFunctionType, BasisfunctionType, 3, dimRange>
  //  RealizabilityLimiterType;
  typedef NonLimitingRealizabilityLimiter<EntityType> RealizabilityLimiterType;
  typedef DefaultEigenSolver<RangeFieldType, dimRange, dimDomain> EigenSolverType;
  auto realizability_limiter = std::make_shared<RealizabilityLimiterType>(*basis_functions, quadrature);
  //  auto realizability_limiter = std::make_shared<RealizabilityLimiterType>();

  // ******************** choose flux and rhs operator and timestepper ******************************************

  typedef typename Dune::XT::Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
      ConstantFunctionType;
  typedef AdvectionRhsOperator<RhsType> RhsOperatorType;

  //  typedef typename std::
  //      conditional<numerical_flux == NumericalFluxes::kinetic,
  //                  AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType>,
  //                  std::conditional<numerical_flux == NumericalFluxes::laxfriedrichs
  //                                       || numerical_flux == NumericalFluxes::laxfriedrichs_with_reconstruction
  //                                       || numerical_flux == NumericalFluxes::local_laxfriedrichs
  //                                       || numerical_flux ==
  //                                       NumericalFluxes::local_laxfriedrichs_with_reconstruction,
  //                                   AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
  //                                                                  BoundaryValueType,
  //                                                                  ConstantFunctionType,
  //                                                                  0,
  //                                                                  SlopeLimiters::minmod,
  //                                                                  false,
  //                                                                  BasisfunctionType>,
  //                                   AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType>>::type>::type
  //          AdvectionOperatorType;

  //  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
  //                                         BoundaryValueType,
  //                                         ConstantFunctionType,
  //                                         1,
  //                                         SlopeLimiters::minmod,
  //                                         EigenSolverType,
  //                                         RealizabilityLimiterType>
  //      AdvectionOperatorType;

  //  typedef AdvectionForceOperator<AnalyticalFluxType,
  //                                 BoundaryValueType,
  //                                 ConstantFunctionType,
  //                                 1,
  //                                 SlopeLimiters::minmod,
  //                                 EigenSolverType,
  //                                 RealizabilityLimiterType>
  //      AdvectionOperatorType;

  //  typedef AdvectionMustaOperator<AnalyticalFluxType,
  //                                 BoundaryValueType,
  //                                 ConstantFunctionType,
  //                                 1,
  //                                 SlopeLimiters::minmod,
  //                                 EigenSolverType,
  //                                 RealizabilityLimiterType>
  //      AdvectionOperatorType;

  typedef AdvectionGodunovOperator<AnalyticalFluxType,
                                   BoundaryValueType,
                                   0,
                                   SlopeLimiters::minmod,
                                   EigenSolverType,
                                   RealizabilityLimiterType>
      AdvectionOperatorType;

  //  typedef AdvectionKineticOperator<AnalyticalFluxType,
  //                                   BoundaryValueType,
  //                                   0,
  //                                   SlopeLimiters::minmod,
  //                                   EigenSolverType,
  //                                   RealizabilityLimiterType>
  //      AdvectionOperatorType;

  //  typedef AdvectionLaxWendroffOperator<AnalyticalFluxType,
  //                                       BoundaryValueType,
  //                                       ConstantFunctionType,
  //                                       0,
  //                                       SlopeLimiters::minmod,
  //                                         EigenSolverType,
  //                                         RealizabilityLimiterType>
  //      AdvectionOperatorType;


  typedef typename TimeStepperFactory<AdvectionOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType
      OperatorTimeStepperType;
  typedef typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, rhs_time_stepper_method>::TimeStepperType
      RhsOperatorTimeStepperType;
  typedef
      typename Dune::GDT::TimeStepperSplittingFactory<RhsOperatorTimeStepperType,
                                                      OperatorTimeStepperType,
                                                      time_stepper_splitting_method>::TimeStepperType TimeStepperType;


  // *************** choose t_end and initial dt **************************************
  // calculate dx and choose initial dt
  Dune::XT::Grid::Dimensions<typename SpaceType::GridLayerType> dimensions(grid_layer);
  RangeFieldType dx = dimensions.entity_width.max();
  if (dimDomain == 2)
    dx /= std::sqrt(2);
  if (dimDomain == 3)
    dx /= std::sqrt(3);
  RangeFieldType dt = CFL * dx;
  t_end = XT::Common::FloatCmp::eq(t_end, 0.) ? problem.t_end() : t_end;


  // *********************** create operators and timesteppers ************************************
  const ConstantFunctionType dx_function(dx);

  //  AdvectionOperatorType advection_operator =
  //      internal::AdvectionOperatorCreator<AdvectionOperatorType, numerical_flux>::create(
  //          analytical_flux, boundary_values, dx_function, linear);

  //  AdvectionOperatorType advection_operator(analytical_flux, boundary_values, dx_function, false, linear);

  //  AdvectionOperatorType advection_operator(analytical_flux,
  //                                           boundary_values,
  //                                           dx_function,
  //                                           AdvectionOperatorType::default_quadrature(),
  //                                           realizability_limiter,
  //                                           false,
  //                                           linear);

  //  AdvectionOperatorType advection_operator(analytical_flux,
  //                                           boundary_values,
  //                                           dx_function,
  //                                           AdvectionOperatorType::default_quadrature(),
  //                                           realizability_limiter,
  //                                           linear);

  //  AdvectionOperatorType advection_operator(analytical_flux,
  //                                           boundary_values,
  //                                           dx_function,
  //                                           AdvectionOperatorType::default_quadrature(),
  //                                           realizability_limiter,
  //                                           linear,
  //                                           10);

  AdvectionOperatorType advection_operator(
      analytical_flux, boundary_values, AdvectionOperatorType::default_1d_quadrature(), realizability_limiter, linear);

  //  AdvectionOperatorType advection_operator(analytical_flux, boundary_values, dx_function, linear);
  //  AdvectionOperatorType advection_operator(analytical_flux, boundary_values, linear);
  //  advection_operator.set_realizability_limiter(realizability_limiter);
  //  advection_operator.set_basisfunctions(basis_functions);
  //  advection_operator.set_quadrature(problem_imp.quadrature());
  //  advection_operator.set_quadrature(quadrature);
  //  advection_operator.set_epsilon(epsilon);
  //  AdvectionOperatorType advection_operator(*analytical_flux,
  //                                           *boundary_values,
  //                                           grid_layer,
  //                                           grid_sizes,
  //                                           plane_coefficients,
  //                                           linear,
  //                                           true,
  //                                           space_quadrature_rules);

  RhsOperatorType rhs_operator(rhs);


  // ******************************** do the time steps ***********************************************************
  OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);
  RhsOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
  TimeStepperType timestepper(timestepper_rhs, timestepper_op);
  //  TimeStepperType timestepper(timestepper_op, timestepper_rhs);
  filename += "_" + ProblemType::static_id();
  filename += Dune::XT::Common::to_string(dimRange);
  filename += rhs_time_stepper_method == TimeStepperMethods::implicit_euler
                  ? "_implicit"
                  : (rhs_time_stepper_method == TimeStepperMethods::matrix_exponential ? "_matexp" : "_explicit");

  timestepper.solve(t_end,
                    dt,
                    num_save_steps,
                    /*save_solution = */ false,
                    /*output_progress = */ true,
                    visualize,
                    filename,
                    basis_functions->visualizer<DiscreteFunctionType>());

  const auto& sol = timestepper.current_solution();
  std::vector<std::pair<DomainType, RangeFieldType>> values;

  for (const auto& entity : Dune::elements(grid_layer)) {
    const auto& local_sol = sol.local_function(entity);
    values.push_back(std::make_pair(entity.geometry().center(),
                                    local_sol->evaluate(entity.geometry().local(entity.geometry().center()))[0]));
  }
  std::sort(values.begin(), values.end(), CmpStruct());
  std::ofstream valuesfile(filename + ".txt");
  for (const auto& pair : values)
    valuesfile << XT::Common::to_string(pair.first, 15) << "\t" << XT::Common::to_string(pair.second, 15) << std::endl;
  valuesfile.close();

  // normalize solution
  RangeFieldType l1norm = 0;
  RangeFieldType l2norm = 0;
  RangeFieldType linfnorm = 0;
  for (const auto& entity : elements(grid_layer)) {
    const auto local_sol = sol.local_function(entity);
    const auto val = local_sol->evaluate(entity.geometry().local(entity.geometry().center()));
    RangeFieldType psi(0);
    //    for (const auto& entry : val) // for hatfunctions
    //      psi += entry;
    for (size_t rr = 0; rr < dimRange; rr += 4) // for piecewise
      psi += val[rr];
    //    psi = val[0] * std::sqrt(4 * M_PI); // for real spherical harmonics

    l1norm += std::abs(psi) * entity.geometry().volume();
    l2norm += std::pow(psi, 2) * entity.geometry().volume();
    linfnorm = std::max(std::abs(psi), linfnorm);
  }

  l1norm = grid_layer.comm().sum(l1norm);
  l2norm = grid_layer.comm().sum(l2norm);
  linfnorm = grid_layer.comm().max(linfnorm);
  l2norm = std::sqrt(l2norm);
  if (grid_layer.comm().rank() == 0) {
    std::cout << "l1norm: " << l1norm << std::endl;
    std::cout << "l2norm: " << l2norm << std::endl;
    std::cout << "linfnorm: " << linfnorm << std::endl;
  }

  std::ifstream matlabvaluesfile("values_matlab.txt");
  std::string line;
  std::vector<DomainType> x_matlab;
  std::vector<RangeFieldType> values_matlab;
  while (std::getline(matlabvaluesfile, line)) {
    auto tokens = XT::Common::tokenize(line, "\t", boost::algorithm::token_compress_on);
    trim(tokens);
    assert(tokens.size() == 2);
    auto x = XT::Common::from_string<DomainType>(tokens[0]);
    x_matlab.push_back(x);
    auto val_matlab = XT::Common::from_string<RangeFieldType>(tokens[1]);
    values_matlab.push_back(val_matlab);
  }

  const size_t grid_size_ns = XT::Common::from_string<size_t>(grid_size);
  Dune::XT::Grid::EntityInlevelSearch<GridLayerType> entity_search(grid_layer);
  const auto entities = entity_search(x_matlab);
  assert(entities.size() == grid_size_ns * grid_size_ns * grid_size_ns);
  RangeFieldType l2error = 0;
  RangeFieldType l1error = 0;
  RangeFieldType linferror = 0;

  for (size_t ii = 0; ii < entities.size(); ++ii) {
    const auto& entity = entities[ii];
    const auto& point = x_matlab[ii];
    if (entity) {
      const auto local_sol = sol.local_function(*entity);
      const auto val = local_sol->evaluate(entity->geometry().local(point));
      RangeFieldType psi(0);
      //    for (const auto& entry : val) // for hatfunctions
      //      psi += entry;
      for (size_t rr = 0; rr < dimRange; rr += 4) // for piecewise
        psi += val[rr];
      //      psi = val[0] * std::sqrt(4 * M_PI); // for real spherical harmonics

      //      psi /= l1norm; // normalize

      const auto& val_matlab = values_matlab[ii];
      l2error += std::pow(psi - val_matlab, 2) * entity->geometry().volume();
      l1error += std::abs(psi - val_matlab) * entity->geometry().volume();
      linferror = std::max(std::abs(psi - val_matlab), linferror);
    }
  }

  l1error = grid_layer.comm().sum(l1error);
  l2error = grid_layer.comm().sum(l2error);
  linferror = grid_layer.comm().max(linferror);
  l2error = std::sqrt(l2error);
  if (grid_layer.comm().rank() == 0) {
    std::cout << "l2error: " << XT::Common::to_string(l2error) << std::endl;
    std::cout << "l1error: " << XT::Common::to_string(l1error) << std::endl;
    std::cout << "linferror: " << XT::Common::to_string(linferror) << std::endl;
  }

  return 0;
}
