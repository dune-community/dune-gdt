// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH

#include <utility>
#include <vector>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/timestepper/adaptive-rungekutta.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/timestepper/fractional-step.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class ProblemType, class AnsatzSpaceType, class MatrixType, class VectorType,
          class TestSpaceType = AnsatzSpaceType>
class StationaryContainerBasedDefaultDiscretization;

template <class ProblemImp, class FVSpaceImp, bool use_lax_friedrichs_flux, bool use_adaptive_timestepper,
          bool use_linear_reconstruction>
class InStationaryDefaultDiscretization;


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


template <class TestCaseImp, class FVSpaceImp, bool use_lax_friedrichs_flux, bool use_adaptive_timestepper,
          bool use_linear_reconstruction>
class InStationaryDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef InStationaryDefaultDiscretization<TestCaseImp, FVSpaceImp, use_lax_friedrichs_flux, use_adaptive_timestepper,
                                            use_linear_reconstruction>
      derived_type;
  typedef typename TestCaseImp::ProblemType ProblemType;
  typedef FVSpaceImp FVSpaceType;
  typedef typename FVSpaceType::RangeFieldType RangeFieldType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef DiscreteFunction<FVSpaceType, VectorType> DiscreteFunctionType;
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

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm, AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp, MatrixType system_mtrx, VectorType rhs_vec,
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

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm, AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx, VectorType rhs_vec, VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {
  }

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm, AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp, MatrixType system_mtrx, VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {
  }

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm, AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx, VectorType rhs_vec)
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
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
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


template <class OperatorType, bool use_lax_friedrichs_flux = false>
struct AdvectionOperatorCreator
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux, const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& /*dx_function*/, const RangeFieldType /*dt*/,
                             const bool is_linear, const bool use_linear_reconstruction)
  {
    return OperatorType(analytical_flux, boundary_values, is_linear, use_linear_reconstruction);
  }
};

template <class OperatorType>
struct AdvectionOperatorCreator<OperatorType, true>
{
  template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, class RangeFieldType>
  static OperatorType create(const AnalyticalFluxType& analytical_flux, const BoundaryValueType& boundary_values,
                             const ConstantFunctionType& dx_function, const RangeFieldType dt, const bool is_linear,
                             const bool /*use_linear_reconstruction*/)
  {
    return OperatorType(analytical_flux, boundary_values, dx_function, dt, is_linear);
  }
};


} // namespace internal


template <class TestCaseImp, class FVSpaceImp, bool use_lax_friedrichs_flux, bool use_adaptive_timestepper,
          bool use_linear_reconstruction>
class InStationaryDefaultDiscretization
    : public NonStationaryDiscretizationInterface<internal::
                                                      InStationaryDefaultDiscretizationTraits<TestCaseImp, FVSpaceImp,
                                                                                              use_lax_friedrichs_flux,
                                                                                              use_adaptive_timestepper,
                                                                                              use_linear_reconstruction>>
{
  typedef NonStationaryDiscretizationInterface<internal::
                                                   InStationaryDefaultDiscretizationTraits<TestCaseImp, FVSpaceImp,
                                                                                           use_lax_friedrichs_flux,
                                                                                           use_adaptive_timestepper,
                                                                                           use_linear_reconstruction>>
      BaseType;
  typedef InStationaryDefaultDiscretization<TestCaseImp, FVSpaceImp, use_lax_friedrichs_flux, use_adaptive_timestepper,
                                            use_linear_reconstruction>
      ThisType;

public:
  typedef TestCaseImp TestCaseType;
  using typename BaseType::ProblemType;
  using typename BaseType::FVSpaceType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::VectorType;
  using typename BaseType::DiscreteFunctionType;

  InStationaryDefaultDiscretization(const TestCaseImp& tst_cs, const std::shared_ptr<const FVSpaceType> fv_space_ptr)
    : test_case_(tst_cs)
    , fv_space_(fv_space_ptr)
  {
  }

  /// \name Required by NonStationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return test_case_.problem();
  }

  const FVSpaceType& fv_space() const
  {
    return *fv_space_;
  }

  using BaseType::solve;

  void solve(DiscreteSolutionType& solution, const bool is_linear) const
  {
#if HAVE_EIGEN
    try {
      // set dimensions
      static const size_t dimDomain = ProblemType::dimDomain;

      // get analytical flux, initial and boundary values
      typedef typename ProblemType::FluxType AnalyticalFluxType;
      typedef typename ProblemType::RHSType RHSType;
      typedef typename ProblemType::InitialValueType InitialValueType;
      typedef typename ProblemType::BoundaryValueType BoundaryValueType;
      typedef typename ProblemType::DomainFieldType DomainFieldType;
      typedef typename ProblemType::RangeFieldType RangeFieldType;
      const std::shared_ptr<const AnalyticalFluxType> analytical_flux = problem().flux();
      const std::shared_ptr<const InitialValueType> initial_values    = problem().initial_values();
      const std::shared_ptr<const BoundaryValueType> boundary_values  = problem().boundary_values();
      const std::shared_ptr<const RHSType> rhs                        = problem().rhs();

      // allocate a discrete function for the concentration and another one to temporary store the update in each step
      typedef DiscreteFunction<FVSpaceType, Dune::Stuff::LA::CommonDenseVector<RangeFieldType>> FVFunctionType;
      FVFunctionType u(*fv_space_, "solution");

      // project initial values
      project(*initial_values, u);

      RangeFieldType t_end = test_case_.t_end();

      const RangeFieldType CFL = problem().CFL();

      // calculate dx and choose t_end and initial dt
      Dune::Stuff::Grid::Dimensions<typename FVSpaceType::GridViewType> dimensions(fv_space_->grid_view());
      RangeFieldType dx = dimensions.entity_width.max();
      if (dimDomain == 2)
        dx /= std::sqrt(2);
      RangeFieldType dt = CFL * dx;

      // define operator types
      typedef typename Dune::Stuff::Functions::
          Constant<typename FVSpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
              ConstantFunctionType;
      typedef
          typename std::conditional<use_lax_friedrichs_flux,
                                    typename Dune::GDT::AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                                                                       BoundaryValueType,
                                                                                       ConstantFunctionType>,
                                    typename Dune::GDT::AdvectionGodunovOperator<AnalyticalFluxType,
                                                                                 BoundaryValueType>>::type OperatorType;
      typedef typename Dune::GDT::AdvectionRHSOperator<RHSType> RHSOperatorType;

      // create right hand side operator
      RHSOperatorType rhs_operator(*rhs);

      const ConstantFunctionType dx_function(dx);

      // create advection operator
      OperatorType advection_operator =
          internal::AdvectionOperatorCreator<OperatorType, use_lax_friedrichs_flux>::create(
              *analytical_flux, *boundary_values, dx_function, dt, is_linear, use_linear_reconstruction);

      typedef
          typename std::conditional<use_adaptive_timestepper,
                                    typename Dune::GDT::AdaptiveRungeKuttaTimeStepper<OperatorType, FVFunctionType>,
                                    typename Dune::GDT::ExplicitRungeKuttaTimeStepper<OperatorType, FVFunctionType>>::
              type OperatorTimeStepperType;
      OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

      // do the time steps
      const size_t num_save_steps = 100;
      solution.clear();
      if (problem().has_non_zero_rhs()) {
        // use fractional step method
        typedef typename std::
            conditional<use_adaptive_timestepper,
                        typename Dune::GDT::AdaptiveRungeKuttaTimeStepper<RHSOperatorType, FVFunctionType>,
                        typename Dune::GDT::ExplicitRungeKuttaTimeStepper<RHSOperatorType, FVFunctionType>>::type
                RHSOperatorTimeStepperType;
        RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
        typedef typename Dune::GDT::FractionalTimeStepper<OperatorTimeStepperType, RHSOperatorTimeStepperType>
            TimeStepperType;
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
    static_assert(AlwaysFalse<DiscreteSolutionType>::value, "You are missing eigen!");
#endif // HAVE_EIGEN
  }

  /// \}

private:
  const TestCaseType& test_case_;
  const std::shared_ptr<const FVSpaceType> fv_space_;
}; // class InStationaryDefaultDiscretization


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
