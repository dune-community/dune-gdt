// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_FACTORY_HH
#define DUNE_GDT_TIMESTEPPER_FACTORY_HH

#include <dune/gdt/timestepper/interface.hh>
#include <dune/gdt/timestepper/adaptive-rungekutta.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/timestepper/implicit-rungekutta.hh>
#include <dune/gdt/timestepper/matrix-exponential.hh>
#include <dune/gdt/timestepper/fractional-step.hh>


namespace Dune {
namespace GDT {


template <class OperatorImp, class DiscreteFunctionImp, TimeStepperMethods method>
struct TimeStepperFactory
{
  typedef typename DiscreteFunctionImp::RangeFieldType RangeFieldType;

  typedef typename std::
      conditional<method == TimeStepperMethods::bogacki_shampine || method == TimeStepperMethods::dormand_prince
                      || method == TimeStepperMethods::adaptive_rungekutta_other,
                  typename Dune::GDT::AdaptiveRungeKuttaTimeStepper<OperatorImp, DiscreteFunctionImp, method>,
                  typename std::
                      conditional<method == TimeStepperMethods::implicit_euler
                                      || method == TimeStepperMethods::implicit_midpoint
                                      || method == TimeStepperMethods::trapezoidal_rule,
                                  typename Dune::GDT::
                                      DiagonallyImplicitRungeKuttaTimeStepper<OperatorImp, DiscreteFunctionImp, method>,
                                  typename std::
                                      conditional<method == TimeStepperMethods::matrix_exponential,
                                                  typename Dune::GDT::MatrixExponentialTimeStepper<OperatorImp,
                                                                                                   DiscreteFunctionImp>,
                                                  typename Dune::GDT::ExplicitRungeKuttaTimeStepper<OperatorImp,
                                                                                                    DiscreteFunctionImp,
                                                                                                    method>>::type>::
                          type>::type TimeStepperType;


  static std::unique_ptr<TimeStepperType> create(const OperatorImp& op,
                                                 const DiscreteFunctionImp& initial_values,
                                                 const RangeFieldType r = 1.0,
                                                 const RangeFieldType t_0 = 0.0,
                                                 const RangeFieldType tol = 1e-4)
  {
    return XT::Common::make_unique<TimeStepperType>(op, initial_values, r, t_0, tol);
  }
}; // class TimeStepperFactory


template <class FirstTimeStepperType, class SecondTimeStepperType, TimeStepperSplittingMethods method>
struct TimeStepperSplittingFactory
{
  typedef typename std::conditional<method == TimeStepperSplittingMethods::fractional_step,
                                    FractionalTimeStepper<FirstTimeStepperType, SecondTimeStepperType>,
                                    StrangSplittingTimeStepper<FirstTimeStepperType, SecondTimeStepperType>>::type
      TimeStepperType;


  static std::unique_ptr<TimeStepperType> create(const FirstTimeStepperType& first_stepper,
                                                 const SecondTimeStepperType& second_stepper)
  {
    return XT::Common::make_unique<TimeStepperType>(first_stepper, second_stepper);
  }
}; // class TimeStepperFactory


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_FACTORY_HH
