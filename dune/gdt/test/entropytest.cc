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

#include <dune/gdt/discretefunction/datahandle.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/fractional-step.hh>
#include <dune/gdt/timestepper/implicit-rungekutta.hh>
#include <dune/gdt/timestepper/implicit-rungekutta-parallel.hh>
#include <dune/gdt/timestepper/matrix_exponential.hh>
//#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>
//#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/planesource.hh>
//#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/pointsource.hh>
//#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/checkerboard3d.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/basisfunctions.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/twobeams.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>

#include <libqhullcpp/RboxPoints.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/QhullQh.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullLinkedList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/Qhull.h>

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
  static const int dimDomain = 1;
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
  const auto rhs_time_stepper_method = TimeStepperMethods::implicit_euler;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::matrix_exponential;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::trapezoidal_rule;
  const auto container_backend = Dune::XT::LA::default_sparse_backend;

  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::Codim<0>::Entity EntityType;

  //******************** choose BasisfunctionType *****************************************
  //  typedef typename Hyperbolic::Problems::LegendrePolynomials<double, dimDomain, double, momentOrder>
  //  BasisfunctionType;
  typedef typename Hyperbolic::Problems::HatFunctions<double, dimDomain, double, momentOrder> BasisfunctionType;
  BasisfunctionType basis_functions;
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
  typedef typename Hyperbolic::Problems::
      SourceBeamPn<BasisfunctionType, EntityType, double, dimDomain, DiscreteFunctionType, double, dimRange>
          ProblemImp;
  //  typedef typename Hyperbolic::Problems::SourceBeamMn<GridViewType,
  //                                                      BasisfunctionType,
  //                                                      EntityType,
  //                                                      double,
  //                                                      dimDomain,
  //                                                      DiscreteFunctionType,
  //                                                      double,
  //                                                      dimRange>
  //      ProblemImp;
  typedef typename Hyperbolic::Problems::KineticEquation<ProblemImp> ProblemType;

  //  typedef typename Hyperbolic::Problems::SourceBeamPnHatFunctions<EntityType, double, dimDomain, double,
  //  momentOrder>
  //      ProblemType;
  //  typedef
  //      typename Hyperbolic::Problems::SourceBeamPnLegendre<EntityType, double, dimDomain, double,
  //      momentOrder>
  //          ProblemType;
  //  typedef typename Hyperbolic::Problems::
  //      SourceBeamPnLegendreLaplaceBeltrami<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //    typedef typename Hyperbolic::Problems::SourceBeamPnFirstOrderDG<EntityType, double, dimDomain, double,
  //    momentOrder>
  //      ProblemType;
  //  typedef typename Hyperbolic::Problems::
  //      PlaneSourcePnLegendre<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //  typedef typename Hyperbolic::Problems::PlaneSourcePnHatFunctions<EntityType, double, dimDomain, double,
  //  momentOrder>
  //      ProblemType;
  //  typedef typename Hyperbolic::Problems::PlaneSourcePnFirstOrderDG<EntityType, double, dimDomain, double,
  //  momentOrder>
  //      ProblemType;
  //    typedef typename Hyperbolic::Problems::
  //        PointSourcePnLegendre<EntityType, double, dimDomain, double, momentOrder>
  //            ProblemType;
  // typedef typename Hyperbolic::Problems::PointSourcePnHatFunctions<EntityType, double, dimDomain, double, 6>
  //    ProblemType;
  //  typedef typename Hyperbolic::Problems::PointSourcePnPartialMoments<EntityType, double, dimDomain, double, 8>
  //      ProblemType;
  //  typedef typename Hyperbolic::Problems::CheckerboardPnHatFunctions<EntityType, double, dimDomain, double, 6>
  //      ProblemType;

  // typedef typename Hyperbolic::Problems::CheckerboardPnPartialMoments<EntityType, double, dimDomain, double, 8>
  //    ProblemType;

  //******************* get typedefs and constants from ProblemType **********************//
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

  // ***************** get quadrature rule *********************************************


  //    // Lebedev quadrature on unit sphere (in polar coordinates)
  //    const size_t quadrature_order = 20;
  //    const auto quadrature_rule = Hyperbolic::Problems::get_lebedev_quadrature(quadrature_order);

  // 3D quadrature on sphere (from http://www.unizar.es/galdeano/actas_pau/PDFVIII/pp61-69.pdf)

  // 3d adaptive quadrature on sphere (from http://www.unizar.es/galdeano/actas_pau/PDFVIII/pp61-69.pdf)
  //  typedef typename GDT::Hyperbolic::Problems::AdaptiveQuadrature<DomainType, RangeType, RangeType>
  //      AdaptiveQuadratureType;
  //  typedef typename AdaptiveQuadratureType::QuadraturePointType QuadraturePointType;
  //  typedef std::function<RangeType(const QuadraturePointType&)> BasisfunctionsType;
  //  std::function<RangeType(const DomainType&, const CGALWrapper::Polyhedron_3&)> basisevaluation =
  //      GDT::Hyperbolic::Problems::evaluate_linear_partial_basis<RangeType, DomainType, CGALWrapper::Polyhedron_3>;
  //  std::function<RangeType(const DomainType&, const CGALWrapper::Polyhedron_3&)> basisevaluation =
  //      GDT::Hyperbolic::Problems::evaluate_spherical_barycentric_coordinates<RangeType,
  //                                                                            DomainType,
  //                                                                            CGALWrapper::Polyhedron_3>;
  //  BasisfunctionsType basisfunctions(
  //      [&](const QuadraturePointType& quadpoint) { return basisevaluation(quadpoint.position(), poly); });
  //  AdaptiveQuadratureType adaptive_quadrature(poly, basisfunctions, rel_tol, abs_tol);

  //******************* create ProblemType object ***************************************
  //  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config));
  //  const auto problem_ptr = ProblemType::create(
  //      quadrature_rule, triangulation, ProblemType::default_config(grid_config, quadrature_rule, triangulation));
  //  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config, true));
  const ProblemImp problem_imp(basis_functions);
  //  const ProblemImp problem_imp(basis_functions, grid_view);
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

  // ********************** store evaluation of basisfunctions at quadrature points in matrix **********************
  // ********************** (for non-adaptive quadratures)                                    **********************
  //  const auto& basis_values_matrix = basis_values_matrices.back();

  // ********************* calculate half-space representation of realizable set **********************
  //  using orgQhull::Qhull;
  //  Qhull qhull;
  //  void     runQhull(const char *inputComment2, int pointDimension, int pointCount, const realT *pointCoordinates,
  //  const char *qhullCommand2);
  // realT is double if the constant REALfloat is not set to 1 manually
  //  std::vector<FieldVector<double, dimRange>> points(quadrature_rule.size() + 1);
  //  points[0] = FieldVector<double, dimRange>(0);
  //  size_t ii = 1;
  //  for (const auto& basis_value : basis_values_matrix)
  //    points[ii++] = basis_value;

  //  qhull.runQhull("Realizable set", int(dimRange), int(points.size()), &(points[0][0]), "Qt");
  //  qhull.outputQhull("n");
  //  const auto facet_end = qhull.endFacet();
  //  std::vector<std::pair<RangeType, RangeFieldType>> plane_coefficients(qhull.facetList().count());
  //  ii = 0;
  //  for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next()) {
  //    for (size_t jj = 0; jj < dimRange; ++jj)
  //      plane_coefficients[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
  //    plane_coefficients[ii].second = -facet.hyperplane().offset();
  //  }

  //*********************** choose analytical flux *************************************************************

  //  typedef EntropyBasedLocalFlux<GridViewType, EntityType, double, dimDomain, double, dimRange, 1>
  //  AnalyticalFluxType;

  //  typedef AdaptiveEntropyBasedLocalFlux<GridViewType, EntityType, double, dimDomain, double, dimRange, 1>
  //      AnalyticalFluxType;

  //  typedef AdaptiveEntropyBasedLocalFlux<GridViewType, EntityType, double, dimDomain, double, dimRange, 1>
  //      AnalyticalFluxType;

  //  typedef typename EntropyBasedLocalFlux3D<GridViewType,
  //                                                      EntityType,
  //                                                      double,
  //                                                      dimDomain,
  //                                                      double,
  //                                                      dimRange,
  //                                                      1,
  //                                                      container_backend>
  //      AnalyticalFluxType;


  typedef typename ProblemType::FluxType AnalyticalFluxType;

  //  typedef typename EntropyBasedLocalFluxHatFunctions<GridViewType,
  //                                                                typename SpaceType::EntityType,
  //                                                                double,
  //                                                                dimDomain,
  //                                                                double,
  //                                                                dimRange,
  //                                                                1>
  //      AnalyticalFluxType;


  // ************************* create analytical flux object ***************************************

  const AnalyticalFluxType& analytical_flux = problem.flux();

  //  const auto analytical_flux =
  //      std::make_shared<const AnalyticalFluxType>(grid_view, ProblemType::create_equidistant_points());

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, ProblemType::create_equidistant_points());

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_1d_hatfunctions);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_1d_firstorderdg);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_3d_hatfunctions);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rules, basis_values_matrix, isotropic_dist_calculator_3d_hatfunctions);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_3d_partialbasis);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(grid_view,
  // isotropic_dist_calculator_3d_partialbasis,
  //                                                                          isotropic_dist_calculator_3d_hatfunctions,
  //                                                                          adaptive_quadrature);


  // ******************** choose flux and rhs operator and timestepper
  // *************************************************

  typedef typename Dune::XT::Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
      ConstantFunctionType;
  typedef AdvectionRHSOperator<RhsType> RHSOperatorType;

  typedef typename std::
      conditional<numerical_flux == NumericalFluxes::kinetic,
                  AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType>,
                  std::conditional<numerical_flux == NumericalFluxes::laxfriedrichs
                                       || numerical_flux == NumericalFluxes::laxfriedrichs_with_reconstruction
                                       || numerical_flux == NumericalFluxes::local_laxfriedrichs
                                       || numerical_flux == NumericalFluxes::local_laxfriedrichs_with_reconstruction,
                                   AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                                                  BoundaryValueType,
                                                                  ConstantFunctionType,
                                                                  SlopeLimiters::minmod>,
                                   AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType>>::type>::type
          AdvectionOperatorType;

  //  typedef AdvectionLaxFriedrichsWENOOperator<AnalyticalFluxType,
  //                                             BoundaryValueType,
  //                                             ConstantFunctionType,
  //                                             GridViewType,
  //                                             BasisFunction::hat_functions,
  //                                             1,
  //                                             SlopeLimiters::minmod>
  //      AdvectionOperatorType;

  //  typedef AdvectionGodunovWENOOperator<AnalyticalFluxType,
  //                                       BoundaryValueType,
  //                                       GridViewType,
  //                                       BasisFunctionType::hat_functions,
  //                                       1,
  //                                       SlopeLimiters::minmod>
  //      AdvectionOperatorType;

  //    typedef AdvectionKineticWENOOperator<AnalyticalFluxType,
  //                                         BoundaryValueType,
  //                                         GridViewType,
  //                                         BasisFunctionType::hat_functions,
  //                                         1,
  //                                         SlopeLimiters::minmod>
  //        AdvectionOperatorType;


  typedef
      typename TimeStepperFactory<AdvectionOperatorType, DiscreteFunctionType, RangeFieldType, time_stepper_method>::
          TimeStepperType OperatorTimeStepperType;
  typedef typename TimeStepperFactory<RHSOperatorType,
                                      DiscreteFunctionType,
                                      RangeFieldType,
                                      rhs_time_stepper_method,
                                      Dune::XT::LA::default_sparse_backend>::TimeStepperType RHSOperatorTimeStepperType;
  //  typedef MatrixExponentialTimeStepper<RHSOperatorType, DiscreteFunctionType, RangeFieldType>
  //      RHSOperatorTimeStepperType;
  //  typedef FractionalTimeStepper<OperatorTimeStepperType, RHSOperatorTimeStepperType> TimeStepperType;
  typedef StrangSplittingTimeStepper<RHSOperatorTimeStepperType, OperatorTimeStepperType> TimeStepperType;


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

  // *********************** get quadrature rules for the FV reconstruction ***********************

  // get 1D quadrature rules
  Dune::QuadratureRule<double, 1> space_quadrature_rule;
  space_quadrature_rule.push_back(Dune::QuadraturePoint<double, 1>(0.5 * (1. - 1. / std::sqrt(3)), 0.5));
  space_quadrature_rule.push_back(Dune::QuadraturePoint<double, 1>(0.5 * (1. + 1. / std::sqrt(3)), 0.5));
  FieldVector<Dune::QuadratureRule<double, 1>, dimDomain> space_quadrature_rules;
  std::fill(space_quadrature_rules.begin(), space_quadrature_rules.end(), space_quadrature_rule);


  // *********************** create operators and timesteppers ************************************
  const ConstantFunctionType dx_function(dx);

  AdvectionOperatorType advection_operator =
      internal::AdvectionOperatorCreator<AdvectionOperatorType, numerical_flux>::create(
          analytical_flux, boundary_values, dx_function, linear);

  FieldVector<size_t, dimDomain> grid_sizes;
  std::fill(grid_sizes.begin(), grid_sizes.end(), XT::Common::from_string<size_t>(grid_size));
  //  AdvectionOperatorType advection_operator(*analytical_flux,
  //                                           *boundary_values,
  //                                           dx_function,
  //                                           grid_view,
  //                                           grid_sizes,
  //                                           plane_coefficients,
  //                                           linear,
  //                                           true,
  //                                           space_quadrature_rules,
  //                                           epsilon,
  //                                           false,
  //                                           DomainType(0));

  //  AdvectionOperatorType advection_operator(*analytical_flux,
  //                                           *boundary_values,
  //                                           grid_view,
  //                                           grid_sizes,
  //                                           plane_coefficients,
  //                                           linear,
  //                                           true,
  //                                           space_quadrature_rules);

  RHSOperatorType rhs_operator(rhs);

  // create timestepper
  OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

  // ******************************** do the time steps ***********************************************************
  // use fractional step method
  RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
  TimeStepperType timestepper(timestepper_rhs, timestepper_op);
  //    std::string filename = ProblemType::static_id() + "WENO";
  //    std::string filename = ProblemType::static_id() + "MinMod";
  //    std::string filename = ProblemType::static_id();
  filename += "_" + ProblemType::static_id();
  filename += Dune::XT::Common::to_string(momentOrder);
  filename += rhs_time_stepper_method == TimeStepperMethods::implicit_euler
                  ? "_implicit"
                  : (rhs_time_stepper_method == TimeStepperMethods::matrix_exponential ? "_matexp" : "_explicit");

  timestepper.solve(t_end, dt, num_save_steps, false, true, visualize, filename, 0);
  return 0;
}
