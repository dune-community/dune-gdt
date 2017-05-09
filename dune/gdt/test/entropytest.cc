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
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/fractional-step.hh>
#include <dune/gdt/timestepper/implicit-rungekutta.hh>
#include <dune/gdt/timestepper/implicit-rungekutta-parallel.hh>
#include <dune/gdt/timestepper/matrix_exponential.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/planesource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/pointsource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/checkerboard3d.hh>

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
  //  static const int dimDomain = 3;
  static const int dimDomain = 1;
  static const int momentOrder = 7;
  const auto numerical_flux = NumericalFluxes::kinetic;
  //  const auto numerical_flux = NumericalFluxes::godunov;
  //  const auto numerical_flux = NumericalFluxes::laxfriedrichs;
  //  const auto numerical_flux = NumericalFluxes::laxfriedrichs_with_reconstruction;
  //  const auto numerical_flux = NumericalFluxes::local_laxfriedrichs_with_reconstruction;
  //      const auto numerical_flux = NumericalFluxes::local_laxfriedrichs;
  const auto time_stepper_method = TimeStepperMethods::explicit_euler;
  //  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  //  const auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_third_order_ssp;
  const auto rhs_time_stepper_method = TimeStepperMethods::implicit_euler;
  //  const auto rhs_time_stepper_method = TimeStepperMethods::trapezoidal_rule;
  const auto container_backend = Dune::XT::LA::default_sparse_backend;

  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::Codim<0>::Entity EntityType;

  //******************** choose ProblemType ***********************************************

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
  typedef typename Hyperbolic::Problems::PlaneSourcePnHatFunctions<EntityType, double, dimDomain, double, momentOrder>
      ProblemType;
  //  typedef typename Hyperbolic::Problems::PlaneSourcePnFirstOrderDG<EntityType, double, dimDomain, double,
  //  momentOrder>
  //      ProblemType;
  //    typedef typename Hyperbolic::Problems::
  //        PointSourcePnLegendre<EntityType, double, dimDomain, double, momentOrder>
  //            ProblemType;
  //    typedef typename Hyperbolic::Problems::PointSourcePnHatFunctions<EntityType, double, dimDomain, double, 6>
  //        ProblemType;
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
  static const size_t dimRange = ProblemType::dimRange;
  typedef typename ProblemType::RHSType RHSType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::BoundaryValueType BoundaryValueType;
  static const bool linear = ProblemType::linear;

  //******************* create grid and FV space ***************************************
  auto grid_config = ProblemType::default_grid_config();
  grid_config["num_elements"] = grid_size;
  grid_config["overlap_size"] = overlap_size;
  const auto grid_ptr = Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config).grid_ptr();
  const auto& grid = *grid_ptr;
  assert(grid.comm().size() == 1 || grid.overlapSize(0) > 0);
  const GridViewType& grid_view = grid_ptr->leafGridView();
  typedef FvProductSpace<GridViewType, RangeFieldType, dimRange, 1> SpaceType;
  const SpaceType fv_space(grid_view);

  // ***************** get quadrature rule *********************************************

  //  // 1D quadrature that consists of a Gauss-Legendre quadrature on each cell of the velocity grid
  Dune::QuadratureRule<double, dimDomain> quadrature_rule;
  static const int num_cells = 100;
  Dune::FieldVector<double, dimDomain> lower_left(-1);
  Dune::FieldVector<double, dimDomain> upper_right(1);
  static const std::array<int, dimDomain> s({num_cells});
  GridType velocity_grid(lower_left, upper_right, s);
  const auto velocity_grid_view = velocity_grid.leafGridView();
  const size_t quadrature_order = 20;
  for (const auto& entity : elements(velocity_grid_view)) {
    const auto local_quadrature_rule = Dune::QuadratureRules<double, dimDomain>::rule(
        entity.type(), quadrature_order, Dune::QuadratureType::GaussLegendre);
    for (const auto& quad : local_quadrature_rule) {
      quadrature_rule.push_back(Dune::QuadraturePoint<double, dimDomain>(
          entity.geometry().global(quad.position()),
          quad.weight() * entity.geometry().integrationElement(quad.position())));
    }
  }


  //    // Lebedev quadrature on unit sphere (in polar coordinates)
  //    const size_t quadrature_order = 20;
  //    const auto quadrature_rule = Hyperbolic::Problems::get_lebedev_quadrature(quadrature_order);

  //    // 3D quadrature on sphere (from http://www.unizar.es/galdeano/actas_pau/PDFVIII/pp61-69.pdf)
  //    const size_t octaeder_refinements = 0;
  //    std::vector<Dune::XT::Common::FieldVector<double, dimDomain>> initial_points{
  //        {1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}};
  //    const Dune::GDT::Hyperbolic::Problems::SphericalTriangulation<double> triangulation(initial_points,
  //                                                                                        octaeder_refinements);
  //    const size_t max_quadrature_refinements = 4;
  //    Dune::GDT::Hyperbolic::Problems::SphericalTriangulation<double> quadrature_triangulation(initial_points, 0);
  //    std::vector<Dune::QuadratureRule<double, dimDomain>> quadrature_rules(max_quadrature_refinements);
  //    for (size_t ii = 0; ii < max_quadrature_refinements; ++ii) {
  //      quadrature_triangulation.refine();
  //      quadrature_rules[ii] = quadrature_triangulation.quadrature_rule();
  //    }
  //    const auto& quadrature_rule = quadrature_rules.back();

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
  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config));
  //    const auto problem_ptr = ProblemType::create(
  //        quadrature_rule, triangulation, ProblemType::default_config(grid_config, quadrature_rule, triangulation));
  //  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config, true));
  const ProblemType& problem = *problem_ptr;
  const std::shared_ptr<const InitialValueType> initial_values = problem.initial_values();
  const std::shared_ptr<const BoundaryValueType> boundary_values = problem.boundary_values();
  const std::shared_ptr<const RHSType> rhs = problem.rhs();
  const RangeFieldType CFL = problem.CFL();

  // ***************** project initial values to discrete function *********************
  // create a discrete function for the solution
  typedef typename Dune::XT::LA::Container<double, container_backend>::VectorType VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  DiscreteFunctionType u(fv_space, "solution");
  // project initial values
  project(*initial_values, u);


  // *************************** get function to calculate u_iso and alpha_iso from u **********************

  // 1D legendre or 3D spherical harmonics
  auto isotropic_dist_calculator_legendre = [](const typename ProblemType::RangeType& uu) {
    typename ProblemType::RangeType u_iso(0), alpha_iso(0);
    u_iso[0] = uu[0];
    alpha_iso[0] = std::log(uu[0] / (ProblemType::dimDomain == 1 ? 2. : 4. * M_PI));
    return std::make_pair(u_iso, alpha_iso);
  };

  const auto v_points = ProblemType::create_equidistant_points();
  auto isotropic_dist_calculator_1d_firstorderdg = [v_points](const typename ProblemType::RangeType& uu) {
    typename ProblemType::RangeType alpha_iso(0);
    typename ProblemType::RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
      psi_iso += uu[ii];
      alpha_iso[ii] = 1;
    }
    psi_iso /= 2.;
    alpha_iso *= std::log(psi_iso);
    typename ProblemType::RangeType u_iso(0);
    for (size_t ii = 0; ii < ProblemType::dimRange / 2; ++ii) {
      u_iso[2 * ii] = v_points[ii + 1] - v_points[ii];
      u_iso[2 * ii + 1] = (std::pow(v_points[ii + 1], 2) - std::pow(v_points[ii], 2)) / 2.;
    }
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  };

  //  const auto v_points = ProblemType::create_equidistant_points();
  auto isotropic_dist_calculator_1d_hatfunctions = [v_points](const typename ProblemType::RangeType& uu) {
    typename ProblemType::RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < ProblemType::dimRange; ++ii)
      psi_iso += uu[ii];
    psi_iso /= 2.;
    typename ProblemType::RangeType alpha_iso(std::log(psi_iso)), u_iso;
    u_iso[0] = v_points[1] - v_points[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      u_iso[ii] = v_points[ii + 1] - v_points[ii - 1];
    u_iso[dimRange - 1] = v_points[dimRange - 1] - v_points[dimRange - 2];
    u_iso *= psi_iso / 2.;
    return std::make_pair(u_iso, alpha_iso);
  };

  //    const auto basis_integrated = ProblemType::basisfunctions_integrated(quadrature_rule, triangulation);
  //    auto isotropic_dist_calculator_3d_hatfunctions = [basis_integrated](const typename ProblemType::RangeType& uu) {
  //      typename ProblemType::RangeFieldType psi_iso(0);
  //      for (size_t ii = 0; ii < ProblemType::dimRange; ++ii)
  //        psi_iso += uu[ii];
  //      psi_iso /= 4. * M_PI;
  //      ProblemType::RangeType alpha_iso(std::log(psi_iso));
  //      auto u_iso = basis_integrated;
  //      u_iso *= psi_iso;
  //      return std::make_pair(u_iso, alpha_iso);
  //    };

  //    auto isotropic_dist_calculator_3d_partialbasis = [basis_integrated](const typename ProblemType::RangeType& uu) {
  //      typename ProblemType::RangeFieldType psi_iso(0);
  //      ProblemType::RangeType alpha_iso(0);
  //      for (size_t ii = 0; ii < ProblemType::dimRange; ii += 4) {
  //        psi_iso += uu[ii];
  //        alpha_iso[ii] = 1.;
  //      }
  //      psi_iso /= 4. * M_PI;
  //      alpha_iso *= std::log(psi_iso);
  //      auto u_iso = basis_integrated;
  //      u_iso *= psi_iso;
  //      return std::make_pair(u_iso, alpha_iso);
  //    };

  // ********************** store evaluation of basisfunctions at quadrature points in matrix **********************
  // ********************** (for non-adaptive quadratures)                                    **********************
  using BasisValuesMatrixType = std::vector<Dune::FieldVector<double, dimRange>>;
  BasisValuesMatrixType basis_values_matrix(quadrature_rule.size());
  //  std::vector<BasisValuesMatrixType> basis_values_matrices(max_quadrature_refinements);
  //  using BasisValuesMatrixType = std::vector<VectorType>;
  //    BasisValuesMatrixType basis_values_matrix(quadrature_rule.size(), VectorType(dimRange));
  //  for (size_t qq = 0; qq < max_quadrature_refinements; ++qq) {
  //    const auto& current_quadrature = quadrature_rules[qq];
  //    basis_values_matrices[qq].resize(current_quadrature.size());
  for (size_t ii = 0; ii < quadrature_rule.size(); ++ii) {
    //    3D hatfunctions on sphere
    //        const auto hatfunctions_evaluated =
    //            Hyperbolic::Problems::evaluate_spherical_barycentric_coordinates<RangeType, DomainType>(
    //                quadrature_rule[ii].position(), triangulation);

    // 3D partial moments
    //    const auto partial_basis_evaluated =
    //        GDT::Hyperbolic::Problems::evaluate_linear_partial_basis<RangeType,
    //        DomainType>(quadrature_rule[ii].position(),
    //        triangulation);

    for (size_t nn = 0; nn < dimRange; ++nn) {
      //            basis_values_matrix[ii][nn] = hatfunctions_evaluated[nn];
      //        basis_values_matrices[qq][ii][nn]
      //        =
      //        hatfunctions_evaluated[nn];
      //      basis_values_matrix[ii][nn]
      //      =
      //      partial_basis_evaluated[nn];
      //      basis_values_matrix[ii][nn]
      //      =
      //          Hyperbolic::Problems::evaluate_legendre_polynomial(quadrature_rule[ii].position(),
      //          nn);
      basis_values_matrix[ii][nn] = Hyperbolic::Problems::evaluate_hat_function(
          quadrature_rule[ii].position()[0], nn, ProblemType::create_equidistant_points());
      //      basis_values_matrix[ii][nn]
      //      =
      //      Hyperbolic::Problems::evaluate_first_order_dg(
      //          quadrature_rule[ii].position()[0],
      //          nn,
      //          ProblemType::create_equidistant_points());
      //      basis_values_matrix[ii][nn]
      //      =
      //      Hyperbolic::Problems::evaluate_real_spherical_harmonics(
      //          quadrature_rule[ii].position()[0],
      //          quadrature_rule[ii].position()[1],
      //          Hyperbolic::Problems::get_l_and_m(nn).first,
      //          Hyperbolic::Problems::get_l_and_m(nn).second);
    }
  }
  //  }

  //  const auto& basis_values_matrix = basis_values_matrices.back();

  // ********************* calculate half-space representation of realizable set **********************
  using orgQhull::Qhull;
  Qhull qhull;
  //  void     runQhull(const char *inputComment2, int pointDimension, int pointCount, const realT *pointCoordinates,
  //  const char *qhullCommand2);
  // realT is double if the constant REALfloat is not set to 1 manually
  std::vector<FieldVector<double, dimRange>> points(quadrature_rule.size() + 1);
  points[0] = FieldVector<double, dimRange>(0);
  size_t ii = 1;
  for (const auto& basis_value : basis_values_matrix)
    points[ii++] = basis_value;

  qhull.runQhull("Realizable set", int(dimRange), int(points.size()), &(points[0][0]), "Qt");
  //  qhull.outputQhull("n");
  const auto facet_end = qhull.endFacet();
  std::vector<std::pair<RangeType, RangeFieldType>> plane_coefficients(qhull.facetList().count());
  ii = 0;
  for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next()) {
    for (size_t jj = 0; jj < dimRange; ++jj)
      plane_coefficients[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
    plane_coefficients[ii].second = -facet.hyperplane().offset();
  }

  //*********************** choose analytical flux *************************************************************

  typedef EntropyBasedLocalFlux<GridViewType, EntityType, double, dimDomain, double, dimRange, 1> AnalyticalFluxType;

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


  //  typedef typename ProblemType::FluxType AnalyticalFluxType;

  //  typedef typename EntropyBasedLocalFluxHatFunctions<GridViewType,
  //                                                                typename SpaceType::EntityType,
  //                                                                double,
  //                                                                dimDomain,
  //                                                                double,
  //                                                                dimRange,
  //                                                                1>
  //      AnalyticalFluxType;


  // ************************* create analytical flux object ***************************************

  //  const std::shared_ptr<const AnalyticalFluxType> analytical_flux = problem.flux();

  //  const auto analytical_flux =
  //      std::make_shared<const AnalyticalFluxType>(grid_view, ProblemType::create_equidistant_points());

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, ProblemType::create_equidistant_points());

  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_1d_hatfunctions);

  //  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //      grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_1d_firstorderdg);

  //    const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
  //        grid_view, quadrature_rule, basis_values_matrix, isotropic_dist_calculator_3d_hatfunctions);

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
  typedef AdvectionRHSOperator<RHSType> RHSOperatorType;

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
  //                                                                  SlopeLimiters::minmod>,
  //                                   AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType>>::type>::type
  //          AdvectionOperatorType;

  typedef AdvectionLaxFriedrichsWENOOperator<AnalyticalFluxType,
                                             BoundaryValueType,
                                             ConstantFunctionType,
                                             GridViewType,
                                             BasisFunctionType::hat_functions,
                                             1,
                                             SlopeLimiters::minmod>
      AdvectionOperatorType;

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
  //  typedef typename TimeStepperFactory<RHSOperatorType,
  //                                      DiscreteFunctionType,
  //                                      RangeFieldType,
  //                                      rhs_time_stepper_method,
  //                                      Dune::XT::LA::default_sparse_backend>::TimeStepperType
  //                                      RHSOperatorTimeStepperType;
  typedef MatrixExponentialTimeStepper<RHSOperatorType, DiscreteFunctionType, RangeFieldType>
      RHSOperatorTimeStepperType;
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

  //    AdvectionOperatorType advection_operator =
  //        internal::AdvectionOperatorCreator<AdvectionOperatorType, numerical_flux>::create(
  //            *analytical_flux, *boundary_values, dx_function, linear);

  FieldVector<size_t, dimDomain> grid_sizes;
  std::fill(grid_sizes.begin(), grid_sizes.end(), XT::Common::from_string<size_t>(grid_size));
  AdvectionOperatorType advection_operator(*analytical_flux,
                                           *boundary_values,
                                           dx_function,
                                           grid_view,
                                           grid_sizes,
                                           plane_coefficients,
                                           linear,
                                           true,
                                           space_quadrature_rules,
                                           epsilon,
                                           true,
                                           DomainType(1));

  //    AdvectionOperatorType advection_operator(*analytical_flux,
  //                                             *boundary_values,
  //                                             grid_view,
  //                                             grid_sizes,
  //                                             plane_coefficients,
  //                                             linear,
  //                                             true,
  //                                             space_quadrature_rules);

  RHSOperatorType rhs_operator(*rhs);

  // create timestepper
  OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

  // ******************************** do the time steps ***********************************************************
  if (problem.has_non_zero_rhs()) {
    // use fractional step method
    RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
    TimeStepperType timestepper(timestepper_rhs, timestepper_op);
    //    std::string filename = ProblemType::static_id() + "WENO";
    //    std::string filename = ProblemType::static_id() + "MinMod";
    //    std::string filename = ProblemType::static_id();
    filename += "_" + ProblemType::static_id();
    filename +=
        std::string("_") + (std::is_same<typename ProblemType::FluxType, AnalyticalFluxType>::value ? "p" : "m");
    filename += Dune::XT::Common::to_string(dimRange);
    filename += rhs_time_stepper_method == TimeStepperMethods::implicit_euler ? "_implicit" : "_explicit";

    timestepper.solve(t_end, dt, num_save_steps, false, true, visualize, filename, 1);
  } else {
    timestepper_op.solve(t_end, dt, num_save_steps, false, true, visualize, "entropy_implicit_trapezoidal");
  }
  return 0;
}
