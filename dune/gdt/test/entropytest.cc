// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#include <stdio.h>
#include <stdlib.h>
#include <utility>
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
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/planesource.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/pointsource.hh>

struct cmp_pair_struct
{
  template <class Pair1, class Pair2>
  bool operator()(const Pair1& a, const Pair2& b)
  {
    return a.position()[0] < b.position()[0];
  }
};

int main(int argc, char** argv)
{

  size_t num_threads = 1;
  size_t num_save_steps = 1000;
  std::string grid_size("100"), overlap_size("1");
  bool visualize = true;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-num_threads") {
      if (i + 1 < argc) {
        num_threads = Dune::XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "-num_threads option requires one argument." << std::endl;
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
    } else if (std::string(argv[i]) == "--no_visualization") {
      visualize = false;
    }
  }

  DXTC_CONFIG.set("threading.partition_factor", 1u, true);
  Dune::XT::Common::threadManager().set_max_threads(num_threads);

  using namespace Dune::GDT;
#if HAVE_DUNE_FEM
  Dune::Fem::MPIManager::initialize(argc, argv);
#else
  MPIHelper::instance(argc, argv);
#endif

  const int dimDomain = 1;
  const int momentOrder = 3;
  const auto numerical_flux = Dune::GDT::NumericalFluxes::kinetic;
  const auto time_stepper_method = Dune::GDT::TimeStepperMethods::explicit_euler;
  const auto rhs_time_stepper_method = Dune::GDT::TimeStepperMethods::explicit_euler;
  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::Codim<0>::Entity EntityType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      SourceBeamPnHatFunctions<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //  typedef
  //      typename Dune::GDT::Hyperbolic::Problems::SourceBeamPnLegendre<EntityType, double, dimDomain, double,
  //      momentOrder>
  //          ProblemType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      SourceBeamPnLegendreLaplaceBeltrami<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  typedef typename Dune::GDT::Hyperbolic::Problems::
      SourceBeamPnFirstOrderDG<EntityType, double, dimDomain, double, momentOrder>
          ProblemType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      PlaneSourcePnLegendre<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      PlaneSourcePnHatFunctions<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      PlaneSourcePnFirstOrderDG<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  //  typedef typename Dune::GDT::Hyperbolic::Problems::
  //      PointSourcePnLegendre<EntityType, double, dimDomain, double, momentOrder>
  //          ProblemType;
  auto grid_config = ProblemType::default_grid_config();
  grid_config["num_elements"] = grid_size;
  grid_config["overlap_size"] = overlap_size;
  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config));
  //  const auto problem_ptr = ProblemType::create(ProblemType::default_config(grid_config, true));
  const ProblemType& problem = *problem_ptr;
  const auto grid_ptr = Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config).grid_ptr();
  const auto& grid = *grid_ptr;
  assert(grid.comm().size() == 1 || grid.overlapSize(0) > 0);
  const GridViewType& grid_view = grid_ptr->leafGridView();
  const auto dimRange = ProblemType::dimRange;
  typedef typename Dune::GDT::FvProductSpace<GridViewType, double, dimRange, 1> SpaceType;
  const SpaceType fv_space(grid_view);

  Dune::FieldVector<double, dimDomain> lower_left(-1);
  Dune::FieldVector<double, dimDomain> upper_right(1);
  static const int num_cells = 10;
  static const std::array<int, dimDomain> s{num_cells};
  GridType velocity_grid(lower_left, upper_right, s);
  const auto velocity_grid_view = velocity_grid.leafGridView();
  const auto quadrature_order = 60;

  Dune::QuadratureRule<double, dimDomain> quadrature_rule;
  for (const auto& entity : elements(velocity_grid_view)) {
    const auto local_quadrature_rule = Dune::QuadratureRules<double, dimDomain>::rule(
        entity.type(), quadrature_order, Dune::QuadratureType::GaussLegendre);
    for (const auto& quad : local_quadrature_rule) {
      quadrature_rule.push_back(Dune::QuadraturePoint<double, dimDomain>(
          entity.geometry().global(quad.position()),
          quad.weight() * entity.geometry().integrationElement(quad.position())));
      std::cout << Dune::XT::Common::to_string(quadrature_rule.back().position()) << " and "
                << quadrature_rule.back().weight() << std::endl;
    }
  }

  //  const auto quadrature_rule = Dune::GDT::Hyperbolic::Problems::get_lebedev_quadrature(quadrature_order);
  //  std::cout << quadrature_rule.size() << std::endl;
  //  return 1;
  //    std::sort(quadrature_rule.begin(), quadrature_rule.end(), cmp_pair_struct{});

  const size_t num_quad_points = num_cells * 31;
  using BasisValuesMatrixType = Dune::FieldMatrix<double, num_quad_points, dimRange>;
  BasisValuesMatrixType basis_values_matrix(0);
  assert(num_quad_points == quadrature_rule.size());
  for (size_t ii = 0; ii < num_quad_points; ++ii) {
    for (size_t nn = 0; nn < dimRange; ++nn) {
      //      basis_values_matrix[ii][nn] =
      //          Dune::GDT::Hyperbolic::Problems::evaluate_legendre_polynomial(quadrature_rule[ii].position(), nn);
      //      basis_values_matrix[ii][nn] = Dune::GDT::Hyperbolic::Problems::evaluate_hat_function(
      //          quadrature_rule[ii].position()[0], nn, ProblemType::create_equidistant_points());
      basis_values_matrix[ii][nn] = Dune::GDT::Hyperbolic::Problems::evaluate_first_order_dg(
          quadrature_rule[ii].position()[0], nn, ProblemType::create_equidistant_points());
      //      basis_values_matrix[ii][nn] = Dune::GDT::Hyperbolic::Problems::evaluate_real_spherical_harmonics(
      //          quadrature_rule[ii].position()[0],
      //          quadrature_rule[ii].position()[1],
      //          Dune::GDT::Hyperbolic::Problems::get_l_and_m(nn).first,
      //          Dune::GDT::Hyperbolic::Problems::get_l_and_m(nn).second);
    }
  }

  typedef typename Dune::XT::LA::Container<double, Dune::XT::LA::default_backend>::VectorType VectorType;
  typedef typename Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;

  static const bool linear = ProblemType::linear;
  //  typedef typename ProblemType::FluxType AnalyticalFluxType;
  typedef typename Dune::GDT::EntropyBasedLocalFlux<GridViewType,
                                                    typename SpaceType::EntityType,
                                                    double,
                                                    dimDomain,
                                                    double,
                                                    dimRange,
                                                    1,
                                                    num_quad_points,
                                                    decltype(ProblemType::create_equidistant_points())>
      AnalyticalFluxType;
  //  typedef typename Dune::GDT::EntropyBasedLocalFluxHatFunctions<GridViewType,
  //                                                                typename SpaceType::EntityType,
  //                                                                double,
  //                                                                dimDomain,
  //                                                                double,
  //                                                                dimRange,
  //                                                                1>
  //      AnalyticalFluxType;
  typedef typename ProblemType::RHSType RHSType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::BoundaryValueType BoundaryValueType;
  typedef typename ProblemType::DomainFieldType DomainFieldType;
  typedef typename ProblemType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::Functions::
      ConstantFunction<typename SpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
          ConstantFunctionType;
  typedef typename Dune::GDT::AdvectionRHSOperator<RHSType> RHSOperatorType;
  typedef typename std::
      conditional<numerical_flux == NumericalFluxes::kinetic,
                  typename Dune::GDT::AdvectionKineticOperator<AnalyticalFluxType, BoundaryValueType>,
                  std::conditional<numerical_flux == NumericalFluxes::laxfriedrichs
                                       || numerical_flux == NumericalFluxes::laxfriedrichs_with_reconstruction
                                       || numerical_flux == NumericalFluxes::local_laxfriedrichs
                                       || numerical_flux == NumericalFluxes::local_laxfriedrichs_with_reconstruction,
                                   typename Dune::GDT::AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                                                                      BoundaryValueType,
                                                                                      ConstantFunctionType>,
                                   typename Dune::GDT::AdvectionGodunovOperator<AnalyticalFluxType,
                                                                                BoundaryValueType>>::type>::type
          AdvectionOperatorType;
  typedef typename Dune::GDT::TimeStepperFactory<AdvectionOperatorType,
                                                 DiscreteFunctionType,
                                                 RangeFieldType,
                                                 time_stepper_method>::TimeStepperType OperatorTimeStepperType;
  typedef typename Dune::GDT::TimeStepperFactory<RHSOperatorType,
                                                 DiscreteFunctionType,
                                                 RangeFieldType,
                                                 rhs_time_stepper_method>::TimeStepperType RHSOperatorTimeStepperType;
  typedef typename Dune::GDT::FractionalTimeStepper<OperatorTimeStepperType, RHSOperatorTimeStepperType>
      TimeStepperType;

  // get analytical flux, initial and boundary values
  //  const std::shared_ptr<const AnalyticalFluxType> analytical_flux = problem.flux();
  //    const auto analytical_flux =
  //        std::make_shared<const AnalyticalFluxType>(grid_view, ProblemType::create_equidistant_points());
  const auto analytical_flux = std::make_shared<const AnalyticalFluxType>(
      grid_view, quadrature_rule, basis_values_matrix, ProblemType::create_equidistant_points());
  const std::shared_ptr<const InitialValueType> initial_values = problem.initial_values();
  const std::shared_ptr<const BoundaryValueType> boundary_values = problem.boundary_values();
  const std::shared_ptr<const RHSType> rhs = problem.rhs();

  // create a discrete function for the solution
  DiscreteFunctionType u(fv_space, "solution");

  // project initial values
  project(*initial_values, u);

  RangeFieldType t_end = problem.t_end();
  const RangeFieldType CFL = problem.CFL();

  // calculate dx and choose initial dt
  Dune::XT::Grid::Dimensions<typename SpaceType::GridViewType> dimensions(grid_view);
  RangeFieldType dx = dimensions.entity_width.max();
  if (dimDomain == 2)
    dx /= std::sqrt(2);
  RangeFieldType dt = 0.5 * CFL * dx;

  // create operators
  const ConstantFunctionType dx_function(dx);
  AdvectionOperatorType advection_operator =
      internal::AdvectionOperatorCreator<AdvectionOperatorType, numerical_flux>::create(
          *analytical_flux, *boundary_values, dx_function, dt, linear);
  RHSOperatorType rhs_operator(*rhs);

  // create timestepper
  OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

  // do the time steps
  if (problem.has_non_zero_rhs()) {
    // use fractional step method
    RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
    TimeStepperType timestepper(timestepper_op, timestepper_rhs);
    std::string filename = ProblemType::static_id();
    filename +=
        std::string("_") + (std::is_same<typename ProblemType::FluxType, AnalyticalFluxType>::value ? "p" : "m");
    filename += Dune::XT::Common::to_string(momentOrder);
    filename += rhs_time_stepper_method == Dune::GDT::TimeStepperMethods::implicit_euler ? "_implicit" : "_explicit";

    timestepper.solve(t_end, dt, num_save_steps, false, true, visualize, filename, 2);
  } else {
    timestepper_op.solve(t_end, dt, num_save_steps, false, true, visualize, "entropy_implicit_trapezoidal");
  }
}
