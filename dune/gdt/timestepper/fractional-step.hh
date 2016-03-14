// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH
#define DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/la/container.hh>

#include "interface.hh"


namespace Dune {
namespace GDT {
namespace TimeStepper {

/**
 * Takes two time steppers and performs a simple fractional step scheme (see e.g. LeVeque, "Finite Volume Methods for
 * Hyperbolic Problems", 2002, Section 17.1. In each time step t, the first stepper is applied up to time t+dt and the
 * result is taken as input for the second time stepper, which is also applied up to time t+dt. Initial values and time
 * are taken from the first time stepper.
 */
template< class FirstStepperImp, class SecondStepperImp >
class FractionalStep
    : public TimeStepperInterface< typename FirstStepperImp::DiscreteFunctionType, typename FirstStepperImp::TimeFieldType >
{
  typedef TimeStepperInterface< typename FirstStepperImp::DiscreteFunctionType, typename FirstStepperImp::TimeFieldType > BaseType;
public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;

  typedef FirstStepperImp FirstStepperType;
  typedef SecondStepperImp SecondStepperType;

  using BaseType::current_solution;
  using BaseType::current_time;
  using BaseType::solution;

  FractionalStep(FirstStepperType& first_stepper, SecondStepperType& second_stepper)
    : BaseType(first_stepper.current_time(), first_stepper.current_solution())
    , first_stepper_(first_stepper)
    , second_stepper_(second_stepper)
  {
    first_stepper_.set_current_solution_pointer(current_solution());
    first_stepper_.set_solution_pointer(solution());
    second_stepper_.set_current_solution_pointer(current_solution());
    second_stepper_.set_solution_pointer(solution());
  } // constructor

  TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt) override final
  {
    auto& t = current_time();
    const TimeFieldType actual_dt = std::min(dt, max_dt);
    const auto dt_1 = first_stepper_.solve(t+actual_dt, dt, -1, false);
    const auto dt_2 = second_stepper_.solve(t+actual_dt, dt, -1, false);
    t += actual_dt;
    return std::min(dt_1, dt_2);
  } // ... step(...)

private:
  FirstStepperType& first_stepper_;
  SecondStepperType& second_stepper_;
}; // class FractionalStep


} // namespace TimeStepper
} // namespace Stuff
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH
