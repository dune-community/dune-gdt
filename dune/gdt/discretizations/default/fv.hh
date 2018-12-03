// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_HH

#include <dune/xt/grid/information.hh>

#include <dune/xt/la/container/common.hh>

#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/projections/l2.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/fractional-step.hh>

#include "fv-internal.hh"
#include "../interfaces.hh"

namespace Dune {
namespace GDT {


template <class ProblemImp,
          class FvSpaceImp,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method,
          TimeStepperSplittingMethods time_stepper_splitting_method,
          class Traits>
class HyperbolicFvDefaultDiscretization;


namespace internal {


template <class TestCaseImp,
          class FvSpaceImp,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method,
          TimeStepperSplittingMethods time_stepper_splitting_method>
class HyperbolicFvDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  using derived_type = HyperbolicFvDefaultDiscretization<TestCaseImp,
                                                         FvSpaceImp,
                                                         numerical_flux,
                                                         reconstruction_order,
                                                         time_stepper_method,
                                                         rhs_time_stepper_method,
                                                         time_stepper_splitting_method,
                                                         HyperbolicFvDefaultDiscretizationTraits>;
  using ProblemType = typename TestCaseImp::ProblemType;
  using SpaceType = FvSpaceImp;
  using RangeFieldType = typename SpaceType::RangeFieldType;
  using VectorType = typename Dune::XT::LA::CommonDenseVector<RangeFieldType>;
  using DiscreteFunctionType = DiscreteFunction<SpaceType, VectorType>;
  using DiscreteSolutionType = std::map<double, DiscreteFunctionType, Dune::GDT::internal::FloatCmpLt>;
}; // class HyperbolicFvDefaultDiscretizationTraits


} // namespace internal


template <class TestCaseImp,
          class FvSpaceImp,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method = time_stepper_method,
          TimeStepperSplittingMethods time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step,
          class Traits = internal::HyperbolicFvDefaultDiscretizationTraits<TestCaseImp,
                                                                           FvSpaceImp,
                                                                           numerical_flux,
                                                                           reconstruction_order,
                                                                           time_stepper_method,
                                                                           rhs_time_stepper_method,
                                                                           time_stepper_splitting_method>>
class HyperbolicFvDefaultDiscretization : public FvDiscretizationInterface<Traits>
{
  static_assert(reconstruction_order <= 1, "Not yet implemented for higher reconstruction orders!");
  using BaseType = FvDiscretizationInterface<Traits>;
  using ThisType = HyperbolicFvDefaultDiscretization;

public:
  using TestCaseType = TestCaseImp;
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::ProblemType;
  using typename BaseType::SpaceType;
  using typename BaseType::VectorType;

private:
  using GridLayerType = typename SpaceType::GridLayerType;
  using IntersectionType = typename GridLayerType::Intersection;

  static const size_t dimDomain = ProblemType::dimDomain;
  using AnalyticalFluxType = typename ProblemType::FluxType;
  using RhsType = typename ProblemType::RhsType;
  using InitialValueType = typename ProblemType::InitialValueType;
  using BoundaryValueType = typename ProblemType::BoundaryValueType;
  using DomainFieldType = typename ProblemType::DomainFieldType;
  using RangeFieldType = typename ProblemType::RangeFieldType;
  using RangeType = typename DiscreteFunctionType::RangeType;
  using ConstantFunctionType = typename Dune::XT::Functions::
      ConstantFunction<typename SpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>;
  using RhsOperatorType = typename Dune::GDT::AdvectionRhsOperator<RhsType>;

  using AdvectionOperatorCreatorType =
      internal::AdvectionOperatorCreator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType, numerical_flux>;
  using AdvectionOperatorType = typename AdvectionOperatorCreatorType::type;
  using NumericalCouplingFluxType = typename AdvectionOperatorType::NumericalCouplingFluxType;
  using NumericalBoundaryFluxType = typename AdvectionOperatorType::NumericalBoundaryFluxType;
  using ReconstructionOperatorType = LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType>;
  using AdvectionWithReconstructionOperatorType =
      AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
  using FvOperatorType =
      std::conditional_t<reconstruction_order == 0, AdvectionOperatorType, AdvectionWithReconstructionOperatorType>;

  using OperatorTimeStepperType =
      typename TimeStepperFactory<FvOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType;
  using RhsOperatorTimeStepperType =
      typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, rhs_time_stepper_method>::TimeStepperType;
  using TimeStepperType =
      typename Dune::GDT::TimeStepperSplittingFactory<OperatorTimeStepperType,
                                                      RhsOperatorTimeStepperType,
                                                      time_stepper_splitting_method>::TimeStepperType;

public:
  HyperbolicFvDefaultDiscretization(const TestCaseImp& tst_cs, const std::shared_ptr<const SpaceType> fv_space_ptr)
    : test_case_(tst_cs)
    , fv_space_(fv_space_ptr)
  {}

  /// \name Required by FvDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return test_case_.problem();
  }

  const SpaceType& space() const
  {
    return *fv_space_;
  }

  using BaseType::solve;

  void solve(DiscreteSolutionType& solution) const
  {
    try {
      const auto& problem = test_case_.problem();
      // get analytical flux, initial and boundary values
      const AnalyticalFluxType& analytical_flux = problem.flux();
      const InitialValueType& initial_values = problem.initial_values();
      const auto& boundary_values = problem.boundary_values();
      const RhsType& rhs = problem.rhs();

      // create a discrete function for the solution
      DiscreteFunctionType u(*fv_space_, "solution");

      // project initial values
      project_l2(initial_values, u, true);

      RangeFieldType t_end = test_case_.t_end();
      const RangeFieldType CFL = problem.CFL();

      // calculate dx and choose initial dt
      Dune::XT::Grid::Dimensions<typename SpaceType::GridLayerType> dimensions(fv_space_->grid_layer());
      RangeFieldType dx = dimensions.entity_width.max();
      if (dimDomain == 2)
        dx /= std::sqrt(2);
      if (dimDomain == 3)
        dx /= std::sqrt(3);
      RangeFieldType dt = CFL * dx;

      // create operators
      const ConstantFunctionType dx_function(dx);
      std::unique_ptr<AdvectionOperatorType> advection_op =
          AdvectionOperatorCreatorType::create(analytical_flux, boundary_values, dx_function);

      MinmodSlope<RangeType, typename ReconstructionOperatorType::MatrixType> slope;
      ReconstructionOperatorType reconstruction_op(analytical_flux, boundary_values, slope);

      AdvectionWithReconstructionOperatorType advection_with_reconstruction_op(*advection_op, reconstruction_op);

      const auto& fv_op = choose_fv_operator(*advection_op, advection_with_reconstruction_op);

      RhsOperatorType rhs_op(rhs);

      // create timestepper
      OperatorTimeStepperType timestepper_op(fv_op, u, -1.0);

      // do the time steps
      const size_t num_save_steps = 100;
      solution.clear();
      if (problem.has_non_zero_rhs()) {
        // use fractional step method
        RhsOperatorTimeStepperType timestepper_rhs(rhs_op, u);
        TimeStepperType timestepper(timestepper_op, timestepper_rhs);
        timestepper.solve(t_end, dt, num_save_steps, solution);
      } else {
        timestepper_op.solve(t_end, dt, num_save_steps, solution);
      }

    } catch (Dune::Exception& e) {
      std::cerr << "Dune reported: " << e.what() << std::endl;
      std::abort();
    }
  } // void solve(...)

  /// \}

  void visualize(const DiscreteSolutionType& solution, const std::string filename) const
  {
    for (size_t ii = 0; ii < solution.size(); ++ii)
      solution[ii].second.visualize(filename, Dune::XT::Common::to_string(ii));
  }

private:
  template <size_t order = reconstruction_order>
  static std::enable_if_t<order == 1, const FvOperatorType&>
  choose_fv_operator(const AdvectionOperatorType& /*advection_op*/,
                     const AdvectionWithReconstructionOperatorType& advection_with_reconstruction_op)
  {
    return advection_with_reconstruction_op;
  }

  template <size_t order = reconstruction_order>
  static std::enable_if_t<order == 0, const FvOperatorType&>
  choose_fv_operator(const AdvectionOperatorType& advection_op,
                     const AdvectionWithReconstructionOperatorType& /*advection_with_reconstruction_op*/)
  {
    return advection_op;
  }

  const TestCaseType& test_case_;
  const std::shared_ptr<const SpaceType> fv_space_;
}; // class HyperbolicFvDefaultDiscretization


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_HH
