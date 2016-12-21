// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#include <utility>
#include <vector>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/la/container/common.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/spaces/fv/product.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/fractional-step.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/sourcebeam.hh>

int main()
{

  using namespace Dune::GDT;

  const int dimDomain = 1;
  const int momentOrder = 5;
  const auto numerical_flux = Dune::GDT::NumericalFluxes::laxfriedrichs;
  const auto time_stepper_method = Dune::GDT::TimeStepperMethods::explicit_euler;
  const auto rhs_time_stepper_method = Dune::GDT::TimeStepperMethods::explicit_euler;
  typedef typename Dune::YaspGrid<dimDomain, Dune::EquidistantOffsetCoordinates<double, dimDomain>> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::Codim<0>::Entity EntityType;
  typedef typename Dune::GDT::Hyperbolic::Problems::SourceBeam<EntityType, double, dimDomain, double, momentOrder>
      ProblemType;
  const auto problem_ptr = ProblemType::create();
  const ProblemType& problem = *problem_ptr;
  const auto grid_ptr =
      Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(ProblemType::default_grid_config()).grid_ptr();
  const GridViewType& grid_view = grid_ptr->leafGridView();
  const auto dimRange = ProblemType::dimRange;
  typedef typename Dune::GDT::FvProductSpace<GridViewType, double, dimRange, 1> SpaceType;
  const SpaceType fv_space(grid_view);

  typedef typename Dune::XT::LA::CommonDenseVector<double> VectorType;
  typedef typename Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;

  static const bool linear = ProblemType::linear;
  typedef typename ProblemType::FluxType AnalyticalFluxType;
  typedef typename ProblemType::RHSType RHSType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::BoundaryValueType BoundaryValueType;
  typedef typename ProblemType::DomainFieldType DomainFieldType;
  typedef typename ProblemType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::Functions::
      ConstantFunction<typename SpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
          ConstantFunctionType;
  typedef typename Dune::GDT::AdvectionRHSOperator<RHSType> RHSOperatorType;
  typedef
      typename std::conditional<numerical_flux == NumericalFluxes::laxfriedrichs
                                    || numerical_flux == NumericalFluxes::laxfriedrichs_with_reconstruction
                                    || numerical_flux == NumericalFluxes::local_laxfriedrichs
                                    || numerical_flux == NumericalFluxes::local_laxfriedrichs_with_reconstruction,
                                typename Dune::GDT::AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                                                                   BoundaryValueType,
                                                                                   ConstantFunctionType>,
                                typename Dune::GDT::AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType>>::
          type AdvectionOperatorType;
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
  const std::shared_ptr<const AnalyticalFluxType> analytical_flux = problem.flux();
  const std::shared_ptr<const InitialValueType> initial_values = problem.initial_values();
  const std::shared_ptr<const BoundaryValueType> boundary_values = problem.boundary_values();
  const std::shared_ptr<const RHSType> rhs = problem.rhs();

  // create a discrete function for the solution
  DiscreteFunctionType u(fv_space, "solution");

  // project initial values
  project(*initial_values, u);

  RangeFieldType t_end = 1;
  const RangeFieldType CFL = problem.CFL();

  // calculate dx and choose initial dt
  Dune::XT::Grid::Dimensions<typename SpaceType::GridViewType> dimensions(grid_view);
  RangeFieldType dx = dimensions.entity_width.max();
  if (dimDomain == 2)
    dx /= std::sqrt(2);
  RangeFieldType dt = CFL * dx;

  // create operators
  const ConstantFunctionType dx_function(dx);
  AdvectionOperatorType advection_operator =
      internal::AdvectionOperatorCreator<AdvectionOperatorType, numerical_flux>::create(
          *analytical_flux, *boundary_values, dx_function, dt, linear);
  RHSOperatorType rhs_operator(*rhs);

  // create timestepper
  OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

  // do the time steps
  const size_t num_save_steps = 100;
  if (problem.has_non_zero_rhs()) {
    // use fractional step method
    RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
    TimeStepperType timestepper(timestepper_op, timestepper_rhs);
    timestepper.solve(t_end, dt, num_save_steps);
  } else {
    timestepper_op.solve(t_end, dt, num_save_steps);
  }
}
