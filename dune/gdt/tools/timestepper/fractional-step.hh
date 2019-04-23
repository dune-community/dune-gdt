// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH
#define DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/la/container.hh>

#include "interface.hh"


namespace Dune {
namespace GDT {


/**
 * Takes two time steppers and performs a simple fractional step scheme (see e.g. LeVeque, "Finite Volume Methods for
 * Hyperbolic Problems", 2002, Section 17.1). In each time step t, the first stepper is applied up to time t+dt and the
 * result is taken as input for the second time stepper, which is also applied up to time t+dt. Initial values and time
 * are taken from the first time stepper.
 */
template <class FirstStepperImp, class SecondStepperImp>
class FractionalTimeStepper : public TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType>
{
  using BaseType = TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType>;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;

  using FirstStepperType = FirstStepperImp;
  using SecondStepperType = SecondStepperImp;

  using BaseType::current_solution;
  using BaseType::current_time;
  using BaseType::solution;

  FractionalTimeStepper(FirstStepperType& first_stepper, SecondStepperType& second_stepper)
    : BaseType(first_stepper.current_time(), first_stepper.current_solution())
    , first_stepper_(first_stepper)
    , second_stepper_(second_stepper)
  {
    first_stepper_.set_current_solution_pointer(current_solution());
    first_stepper_.set_solution_pointer(solution());
    second_stepper_.set_current_solution_pointer(current_solution());
    second_stepper_.set_solution_pointer(solution());
  } // constructor

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    auto& t = current_time();
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    const auto dt_1 = first_stepper_.solve(t + actual_dt, dt, -1, 0, false);
    const auto dt_2 = second_stepper_.solve(t + actual_dt, dt_1, -1, 0, false);
    t += actual_dt;
    return dt_2;
  } // ... step(...)

  virtual RangeFieldType step_first(const RangeFieldType dt, const RangeFieldType max_dt) override
  {
    auto& t = current_time();
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    return first_stepper_.solve(t + actual_dt, dt, -1, 0, false);
  } // ... step(...)

  virtual RangeFieldType step_second(const RangeFieldType dt, const RangeFieldType actual_dt) override
  {
    auto& t = current_time();
    const auto dt2 = second_stepper_.solve(t + actual_dt, dt, -1, 0, false);
    t += actual_dt;
    return dt2;
  } // ... step(...)

private:
  FirstStepperType& first_stepper_;
  SecondStepperType& second_stepper_;
}; // class FractionalTimeStepper

/**
 * Takes two time steppers and performs a strang splitting fractional step scheme (see e.g. LeVeque, "Finite Volume
 * Methods for
 * Hyperbolic Problems", 2002, Section 17.1). In each time step t, the first stepper is applied up to time t+dt/2, then
 * the
 * second time stepper is applied up to time t+dt, then the first stepper is applied again up to time t+dt. Initial
 * values and time
 * are taken from the first time stepper.
 */
template <class FirstStepperImp, class SecondStepperImp>
class StrangSplittingTimeStepper : public TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType>
{
  using BaseType = TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType>;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;

  using FirstStepperType = FirstStepperImp;
  using SecondStepperType = SecondStepperImp;

  using BaseType::current_solution;
  using BaseType::current_time;
  using BaseType::solution;

  StrangSplittingTimeStepper(FirstStepperType& first_stepper, SecondStepperType& second_stepper)
    : BaseType(first_stepper.current_time(), first_stepper.current_solution())
    , first_stepper_(first_stepper)
    , second_stepper_(second_stepper)
  {
    first_stepper_.set_current_solution_pointer(current_solution());
    first_stepper_.set_solution_pointer(solution());
    second_stepper_.set_current_solution_pointer(current_solution());
    second_stepper_.set_solution_pointer(solution());
  } // constructor

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    auto& t = current_time();
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    first_stepper_.solve(t + actual_dt / 2, actual_dt / 2, static_cast<size_t>(-1), 0, false);
    second_stepper_.solve(t + actual_dt, actual_dt, static_cast<size_t>(-1), 0, false);
    first_stepper_.solve(t + actual_dt, actual_dt / 2, static_cast<size_t>(-1), 0, false);
    t += actual_dt;
    return dt;
  } // ... step(...)

private:
  FirstStepperType& first_stepper_;
  SecondStepperType& second_stepper_;
}; // class StrangSplittingTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_FRACTIONAL_STEP_HH
