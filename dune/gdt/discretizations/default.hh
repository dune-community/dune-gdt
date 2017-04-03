// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH

#include <utility>
#include <vector>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/la/container/common.hh>

#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/timestepper/factory.hh>
#include <dune/gdt/timestepper/fractional-step.hh>
#include <dune/gdt/test/hyperbolic/discretizers/base.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class ProblemType,
          class AnsatzSpaceType,
          class MatrixType,
          class VectorType,
          class TestSpaceType = AnsatzSpaceType>
class StationaryContainerBasedDefaultDiscretization;

template <class ProblemImp,
          class FVSpaceImp,
          NumericalFluxes numerical_flux,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method>
class HyperbolicFVDefaultDiscretization;


namespace internal {


template <class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp>
class StationaryContainerBasedDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef StationaryContainerBasedDefaultDiscretization<ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp>
      derived_type;
  typedef ProblemImp ProblemType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef TestSpaceImp TestSpaceType;
  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;
}; // class StationaryContainerBasedDefaultDiscretizationTraits


template <class TestCaseImp,
          class FVSpaceImp,
          NumericalFluxes numerical_flux,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method>
class HyperbolicFVDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef HyperbolicFVDefaultDiscretization<TestCaseImp,
                                            FVSpaceImp,
                                            numerical_flux,
                                            time_stepper_method,
                                            rhs_time_stepper_method>
      derived_type;
  typedef typename TestCaseImp::ProblemType ProblemType;
  typedef FVSpaceImp SpaceType;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef typename Dune::XT::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef std::map<double, DiscreteFunctionType, Dune::GDT::internal::FloatCmpLt> DiscreteSolutionType;
}; // class InStationaryDefaultDiscretizationTraits


} // namespace internal


template <class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp>
class StationaryContainerBasedDefaultDiscretization
    : public ContainerBasedStationaryDiscretizationInterface<internal::
                                                                 StationaryContainerBasedDefaultDiscretizationTraits<ProblemImp,
                                                                                                                     AnsatzSpaceImp,
                                                                                                                     MatrixImp,
                                                                                                                     VectorImp,
                                                                                                                     TestSpaceImp>>
{
  typedef ContainerBasedStationaryDiscretizationInterface<internal::
                                                              StationaryContainerBasedDefaultDiscretizationTraits<ProblemImp,
                                                                                                                  AnsatzSpaceImp,
                                                                                                                  MatrixImp,
                                                                                                                  VectorImp,
                                                                                                                  TestSpaceImp>>
      BaseType;
  typedef StationaryContainerBasedDefaultDiscretization<ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp>
      ThisType;

public:
  using typename BaseType::ProblemType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec,
                                                VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {
  }

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec,
                                                VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {
  }

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {
  }

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {
  }

  StationaryContainerBasedDefaultDiscretization(ThisType&& /*source*/) = default;

  /// \name Required by StationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return problem_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  /// \}
  /// \name Required by ContainerBasedStationaryDiscretizationInterface.
  /// \{

  const MatrixType& system_matrix() const
  {
    return system_matrix_;
  }

  const VectorType& rhs_vector() const
  {
    return rhs_vector_;
  }

  bool has_dirichlet_shift() const
  {
    return has_dirichlet_shift_;
  }

  const VectorType& dirichlet_shift() const
  {
    if (!has_dirichlet_shift_)
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Do not call dirichlet_shift() if has_dirichlet_shift() is false!");
    return dirichlet_shift_;
  }

  /// \}

private:
  const ProblemType& problem_;
  const AnsatzSpaceType ansatz_space_;
  const TestSpaceType test_space_;
  const MatrixType system_matrix_;
  const VectorType rhs_vector_;
  const VectorType dirichlet_shift_;
  const bool has_dirichlet_shift_;
}; // class StationaryContainerBasedDefaultDiscretization


namespace internal {


template <class OperatorType, NumericalFluxes numerical_flux = NumericalFluxes::godunov>
struct AdvectionOperatorCreator
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& /*dx_function*/,
                             const RangeFieldType /*dt*/,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, is_linear, false);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, NumericalFluxes::godunov_with_reconstruction>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& /*dx_function*/,
                             const RangeFieldType /*dt*/,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, is_linear, true);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, NumericalFluxes::laxfriedrichs>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& dx_function,
                             const RangeFieldType dt,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, dx_function, dt, is_linear, false, false);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, NumericalFluxes::laxfriedrichs_with_reconstruction>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& dx_function,
                             const RangeFieldType dt,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, dx_function, dt, is_linear, true, false);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, NumericalFluxes::local_laxfriedrichs>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& dx_function,
                             const RangeFieldType dt,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, dx_function, dt, is_linear, false, true);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, NumericalFluxes::local_laxfriedrichs_with_reconstruction>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux,
                             const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& dx_function,
                             const RangeFieldType dt,
                             const bool is_linear)
  {
    return OperatorType(analytical_flux, boundary_values, dx_function, dt, is_linear, true, true);
  }
};


} // namespace internal


template <class TestCaseImp,
          class FVSpaceImp,
          NumericalFluxes numerical_flux,
          TimeStepperMethods time_stepper_method,
          TimeStepperMethods rhs_time_stepper_method = time_stepper_method>
class HyperbolicFVDefaultDiscretization
    : public FVDiscretizationInterface<internal::HyperbolicFVDefaultDiscretizationTraits<TestCaseImp,
                                                                                         FVSpaceImp,
                                                                                         numerical_flux,
                                                                                         time_stepper_method,
                                                                                         rhs_time_stepper_method>>
{
  typedef FVDiscretizationInterface<internal::HyperbolicFVDefaultDiscretizationTraits<TestCaseImp,
                                                                                      FVSpaceImp,
                                                                                      numerical_flux,
                                                                                      time_stepper_method,
                                                                                      rhs_time_stepper_method>>
      BaseType;
  typedef HyperbolicFVDefaultDiscretization<TestCaseImp,
                                            FVSpaceImp,
                                            numerical_flux,
                                            time_stepper_method,
                                            rhs_time_stepper_method>
      ThisType;

public:
  typedef TestCaseImp TestCaseType;
  using typename BaseType::SpaceType;
  using typename BaseType::ProblemType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::VectorType;
  using typename BaseType::DiscreteFunctionType;

private:
  static const bool linear = ProblemType::linear;
  static const size_t dimDomain = ProblemType::dimDomain;
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
  typedef
      typename TimeStepperFactory<AdvectionOperatorType, DiscreteFunctionType, RangeFieldType, time_stepper_method>::
          TimeStepperType OperatorTimeStepperType;
  typedef typename TimeStepperFactory<RHSOperatorType, DiscreteFunctionType, RangeFieldType, rhs_time_stepper_method>::
      TimeStepperType RHSOperatorTimeStepperType;
  typedef typename Dune::GDT::FractionalTimeStepper<OperatorTimeStepperType, RHSOperatorTimeStepperType>
      TimeStepperType;

public:
  HyperbolicFVDefaultDiscretization(const TestCaseImp& tst_cs, const std::shared_ptr<const SpaceType> fv_space_ptr)
    : test_case_(tst_cs)
    , fv_space_(fv_space_ptr)
  {
  }

  /// \name Required by FVDiscretizationInterface.
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
#if HAVE_EIGEN
    try {
      // get analytical flux, initial and boundary values
      const std::shared_ptr<const AnalyticalFluxType> analytical_flux = problem().flux();
      const std::shared_ptr<const InitialValueType> initial_values = problem().initial_values();
      const std::shared_ptr<const BoundaryValueType> boundary_values = problem().boundary_values();
      const std::shared_ptr<const RHSType> rhs = problem().rhs();

      // create a discrete function for the solution
      DiscreteFunctionType u(*fv_space_, "solution");

      // project initial values
      project(*initial_values, u);

      RangeFieldType t_end = test_case_.t_end();
      const RangeFieldType CFL = problem().CFL();

      // calculate dx and choose initial dt
      Dune::XT::Grid::Dimensions<typename SpaceType::GridLayerType> dimensions(fv_space_->grid_layer());
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
      solution.clear();
      if (problem().has_non_zero_rhs()) {
        // use fractional step method
        RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
        TimeStepperType timestepper(timestepper_op, timestepper_rhs);
        timestepper.solve(t_end, dt, num_save_steps, solution);
      } else {
        timestepper_op.solve(t_end, dt, num_save_steps, solution);
      }

    } catch (Dune::Exception& e) {
      std::cerr << "Dune reported: " << e.what() << std::endl;
      std::abort();
    }
#else // HAVE_EIGEN
    (void)solution; // silence warning during compilation
    static_assert(AlwaysFalse<DiscreteSolutionType>::value, "You are missing eigen!");
#endif // HAVE_EIGEN
  }

  /// \}

  void visualize(const DiscreteSolutionType& solution, const std::string filename) const
  {
    for (size_t ii = 0; ii < solution.size(); ++ii)
      solution[ii].second.visualize(filename, Dune::XT::Common::to_string(ii));
  }

private:
  const TestCaseType& test_case_;
  const std::shared_ptr<const SpaceType> fv_space_;
}; // class HyperbolicFVDefaultDiscretization


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
