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
          SlopeLimiters slope_limiter,
          class Traits>
class HyperbolicFvDefaultDiscretization;


namespace internal {


template <class TestCaseImp,
          class FvSpaceImp,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method,
          TimeStepperSplittingMethods time_stepper_splitting_method,
          SlopeLimiters slope_limiter>
class HyperbolicFvDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef HyperbolicFvDefaultDiscretization<TestCaseImp,
                                            FvSpaceImp,
                                            numerical_flux,
                                            reconstruction_order,
                                            time_stepper_method,
                                            rhs_time_stepper_method,
                                            time_stepper_splitting_method,
                                            slope_limiter,
                                            HyperbolicFvDefaultDiscretizationTraits>
      derived_type;
  typedef typename TestCaseImp::ProblemType ProblemType;
  typedef FvSpaceImp SpaceType;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef std::map<double, DiscreteFunctionType, Dune::GDT::internal::FloatCmpLt> DiscreteSolutionType;
}; // class HyperbolicFvDefaultDiscretizationTraits


} // namespace internal


template <class TestCaseImp,
          class FvSpaceImp,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method = time_stepper_method,
          TimeStepperSplittingMethods time_stepper_splitting_method = TimeStepperSplittingMethods::fractional_step,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod,
          class Traits = internal::HyperbolicFvDefaultDiscretizationTraits<TestCaseImp,
                                                                           FvSpaceImp,
                                                                           numerical_flux,
                                                                           reconstruction_order,
                                                                           time_stepper_method,
                                                                           rhs_time_stepper_method,
                                                                           time_stepper_splitting_method,
                                                                           slope_limiter>>
class HyperbolicFvDefaultDiscretization : public FvDiscretizationInterface<Traits>
{
  static_assert(reconstruction_order <= 1, "Not yet implemented for higher reconstruction orders!");
  typedef FvDiscretizationInterface<Traits> BaseType;
  typedef HyperbolicFvDefaultDiscretization ThisType;

public:
  typedef TestCaseImp TestCaseType;
  using typename BaseType::SpaceType;
  using typename BaseType::ProblemType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::VectorType;
  using typename BaseType::DiscreteFunctionType;

private:
  using GridLayerType = typename SpaceType::GridLayerType;
  using IntersectionType = typename GridLayerType::Intersection;

  static const size_t dimDomain = ProblemType::dimDomain;
  typedef typename ProblemType::FluxType AnalyticalFluxType;
  typedef typename ProblemType::RhsType RhsType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  using BoundaryValueType =
      LocalizableFunctionBasedLocalizableDirichletBoundaryValue<GridLayerType, typename ProblemType::BoundaryValueType>;
  typedef typename ProblemType::DomainFieldType DomainFieldType;
  typedef typename ProblemType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::Functions::
      ConstantFunction<typename SpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
          ConstantFunctionType;
  typedef typename Dune::GDT::AdvectionRhsOperator<RhsType> RhsOperatorType;

  typedef internal::
      AdvectionOperatorCreator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType, numerical_flux>
          AdvectionOperatorCreatorType;
  typedef typename AdvectionOperatorCreatorType::type AdvectionOperatorType;
  typedef typename AdvectionOperatorType::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename AdvectionOperatorType::NumericalBoundaryFluxType NumericalBoundaryFluxType;
  using ReconstructionOperatorType = LinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType>;
  using AdvectionWithReconstructionOperatorType =
      AdvectionWithReconstructionOperator<AdvectionOperatorType, ReconstructionOperatorType>;
  using FvOperatorType =
      std::conditional_t<reconstruction_order == 0, AdvectionOperatorType, AdvectionWithReconstructionOperatorType>;

  typedef typename TimeStepperFactory<FvOperatorType, DiscreteFunctionType, time_stepper_method>::TimeStepperType
      OperatorTimeStepperType;
  typedef typename TimeStepperFactory<RhsOperatorType, DiscreteFunctionType, rhs_time_stepper_method>::TimeStepperType
      RhsOperatorTimeStepperType;
  typedef
      typename Dune::GDT::TimeStepperSplittingFactory<OperatorTimeStepperType,
                                                      RhsOperatorTimeStepperType,
                                                      time_stepper_splitting_method>::TimeStepperType TimeStepperType;

public:
  HyperbolicFvDefaultDiscretization(const TestCaseImp& tst_cs, const std::shared_ptr<const SpaceType> fv_space_ptr)
    : test_case_(tst_cs)
    , fv_space_(fv_space_ptr)
  {
  }

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
      const auto& dirichlet_boundary_values = problem.boundary_values();
      const auto boundary_info = XT::Grid::AllDirichletBoundaryInfo<IntersectionType>();
      const BoundaryValueType boundary_values(boundary_info, dirichlet_boundary_values);
      const RhsType& rhs = problem.rhs();

      // create a discrete function for the solution
      DiscreteFunctionType u(*fv_space_, "solution");

      // project initial values
      project_l2(initial_values, u);

      RangeFieldType t_end = problem.t_end();
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

      ReconstructionOperatorType reconstruction_op(analytical_flux, boundary_values);

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
