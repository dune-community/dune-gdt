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

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/istl.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshwriter.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
//#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/checkerboard3d.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/basisfunctions.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/twobeams.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/planesource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/pointsource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/linesource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/lebedevquadrature.hh>

#include <Eigen/Core>

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
  const auto numerical_flux = NumericalFluxes::laxfriedrichs;
  //  const auto numerical_flux = NumericalFluxes::laxfriedrichs_with_reconstruction;
  //  const auto numerical_flux = NumericalFluxes::local_laxfriedrichs_with_reconstruction;
  //      const auto numerical_flux = NumericalFluxes::local_laxfriedrichs;
  // const auto time_stepper_method = TimeStepperMethods::explicit_euler;
  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  //  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_third_order_ssp;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::explicit_euler;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::implicit_euler;
  const auto rhs_time_stepper_method = TimeStepperMethods::matrix_exponential;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::trapezoidal_rule;
  const auto container_backend = Dune::XT::LA::default_sparse_backend;

  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::Codim<0>::Entity EntityType;

  //******************** choose BasisfunctionType *****************************************
  //  typedef typename Hyperbolic::Problems::LegendrePolynomials<double, dimDomain, double, momentOrder>
  //  BasisfunctionType;
  //  typedef typename Hyperbolic::Problems::HatFunctions<double, dimDomain, double, momentOrder> BasisfunctionType;

  static const size_t refinements = 0;
  typedef
      typename Hyperbolic::Problems::HatFunctions<double,
                                                  dimDomain,
                                                  double,
                                                  Hyperbolic::Problems::OctaederStatistics<refinements>::num_vertices(),
                                                  1,
                                                  3>
          BasisfunctionType;

  //  typedef typename Hyperbolic::Problems::
  //      HatFunctions<double, 3, double, Hyperbolic::Problems::OctaederStatistics<refinements>::num_vertices(), 1, 2>
  //          BasisfunctionType;

  //  typedef typename Hyperbolic::Problems::RealSphericalHarmonics<double, double, momentOrder, dimDomain, true>
  //      BasisfunctionType;

  //  typedef typename Hyperbolic::Problems::
  //      PiecewiseMonomials<double,
  //                         dimDomain,
  //                         double,
  //                         4 * Hyperbolic::Problems::OctaederStatistics<refinements>::num_faces()>
  //          BasisfunctionType;
  //  BasisfunctionType basis_functions;
  std::shared_ptr<const BasisfunctionType> basis_functions =
      std::make_shared<const BasisfunctionType>(refinements, refinements + 4);
  static const size_t dimRange = BasisfunctionType::dimRange;
  typedef FvProductSpace<GridViewType, double, dimRange, 1> SpaceType;
  typedef typename Dune::XT::LA::Container<double, container_backend>::VectorType VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;

  //******************** choose ProblemType ***********************************************
  //  typedef typename Hyperbolic::Problems::
  //      TwoBeamsFokkerPlanckPn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double,
  //      dimRange>
  //          ProblemImp;

  //    typedef typename Hyperbolic::Problems::
  //      TwoBeamsMn<GridViewType, BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double,
  //      dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::
  //      SourceBeamPn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double, dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::SourceBeamMn<GridViewType,
  //                                                      BasisfunctionType,
  //                                                      EntityType,
  //                                                      double,
  //                                                      dimDomain,
  //                                                      DiscreteFunctionType,
  //                                                      double,
  //                                                      dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::
  //      PlaneSourcePn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double, dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::PlaneSourceMn<GridViewType,
  //                                                       BasisfunctionType,
  //                                                       EntityType,
  //                                                       double,
  //                                                       dimDomain,
  //                                                       DiscreteFunctionType,
  //                                                       double,
  //                                                       dimRange>
  //      ProblemImp;

  //  typedef typename Hyperbolic::Problems::
  //      PointSourcePn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double, dimRange>
  //          ProblemImp;

  typedef typename Hyperbolic::Problems::PointSourceMn<GridViewType,
                                                       BasisfunctionType,
                                                       EntityType,
                                                       double,
                                                       dimDomain,
                                                       DiscreteFunctionType,
                                                       double,
                                                       dimRange>
      ProblemImp;

  //  typedef typename Hyperbolic::Problems::
  //      ModifiedLineSourcePn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double, dimRange>
  //          ProblemImp;

  //  typedef typename Hyperbolic::Problems::ModifiedLineSourceMn<GridViewType,
  //                                                              BasisfunctionType,
  //                                                              EntityType,
  //                                                              double,
  //                                                              dimDomain,
  //                                                              DiscreteFunctionType,
  //                                                              double,
  //                                                              dimRange>
  //      ProblemImp;


  //******************* get typedefs and constants from ProblemType **********************//
  typedef typename Hyperbolic::Problems::KineticEquation<ProblemImp> ProblemType;
  using DomainFieldType = typename ProblemType::DomainFieldType;
  using DomainType = typename ProblemType::DomainType;
  using RangeFieldType = typename ProblemType::RangeFieldType;
  using RangeType = typename ProblemType::RangeType;
  //  static const size_t dimRange = ProblemType::dimRange;
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
  const GridViewType& grid_view = grid_ptr->leafGridView();
  const SpaceType fv_space(grid_view);

  //******************* create ProblemType object ***************************************
  //  const ProblemImp problem_imp(*basis_functions, grid_config);
  //  const ProblemImp problem_imp(basis_functions, grid_view, grid_config);
  const ProblemImp problem_imp(*basis_functions, basis_functions->quadrature(), grid_view, grid_config);
  //  const ProblemImp problem_imp(
  //      basis_functions, Hyperbolic::Problems::LebedevQuadrature<DomainFieldType, true>::get(1000), grid_view);

  //  const ProblemImp problem_imp(basis_functions,
  //                               grid_view,
  //                               ProblemImp::default_grid_cfg(),
  //                               ProblemImp::default_boundary_cfg(),
  //                               //                               basis_functions.quadrature());
  //                               Hyperbolic::Problems::LebedevQuadrature<DomainFieldType, true>::get(100));

  const ProblemType problem(problem_imp);
  const InitialValueType& initial_values = problem.initial_values();
  const BoundaryValueType& boundary_values = problem.boundary_values();
  const RhsType& rhs = problem.rhs();
  const RangeFieldType CFL = problem.CFL();

  // ***************** project initial values to discrete function *********************
  // create a discrete function for the solution
  DiscreteFunctionType u(fv_space, "solution");
  // project initial values
  project_l2(initial_values, u);


  // ************************* create analytical flux object ***************************************
  typedef typename ProblemType::FluxType AnalyticalFluxType;
  const AnalyticalFluxType& analytical_flux = problem.flux();


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
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                         BoundaryValueType,
                                         ConstantFunctionType,
                                         1,
                                         SlopeLimiters::minmod,
                                         true,
                                         BasisfunctionType>
      AdvectionOperatorType;

  typedef
      typename TimeStepperFactory<AdvectionOperatorType, DiscreteFunctionType, RangeFieldType, time_stepper_method>::
          TimeStepperType OperatorTimeStepperType;
  typedef typename TimeStepperFactory<RhsOperatorType,
                                      DiscreteFunctionType,
                                      RangeFieldType,
                                      rhs_time_stepper_method,
                                      Dune::XT::LA::default_sparse_backend>::TimeStepperType RhsOperatorTimeStepperType;
  //    typedef FractionalTimeStepper<OperatorTimeStepperType, RhsOperatorTimeStepperType> TimeStepperType;
  typedef StrangSplittingTimeStepper<RhsOperatorTimeStepperType, OperatorTimeStepperType> TimeStepperType;


  // *************** choose t_end and initial dt **************************************
  // calculate dx and choose initial dt
  Dune::XT::Grid::Dimensions<typename SpaceType::GridViewType> dimensions(grid_view);
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

  FieldVector<size_t, dimDomain> grid_sizes;
  std::fill(grid_sizes.begin(), grid_sizes.end(), XT::Common::from_string<size_t>(grid_size));
  AdvectionOperatorType advection_operator(analytical_flux, boundary_values, dx_function);
  advection_operator.set_basisfunctions(basis_functions);
  //  advection_operator.set_quadrature(problem_imp.quadrature());
  advection_operator.set_quadrature(basis_functions->quadrature());

  //  AdvectionOperatorType advection_operator(*analytical_flux,
  //                                           *boundary_values,
  //                                           grid_view,
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
  filename += "_" + ProblemType::static_id();
  filename += Dune::XT::Common::to_string(dimRange);
  filename += rhs_time_stepper_method == TimeStepperMethods::implicit_euler
                  ? "_implicit"
                  : (rhs_time_stepper_method == TimeStepperMethods::matrix_exponential ? "_matexp" : "_explicit");

  timestepper.solve(t_end, dt, num_save_steps, false, true, visualize, filename, 1);
  return 0;
}
