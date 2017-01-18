// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
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


namespace Dune {
namespace GDT {


template <class OperatorImp, class DiscreteFunctionImp, class TimeFieldImp, TimeStepperMethods method>
struct TimeStepperFactory
{
  typedef typename std::
      conditional<method == TimeStepperMethods::bogacki_shampine || method == TimeStepperMethods::dormand_prince
                      || method == TimeStepperMethods::adaptive_rungekutta_other,
                  typename Dune::GDT::
                      AdaptiveRungeKuttaTimeStepper<OperatorImp, DiscreteFunctionImp, TimeFieldImp, method>,
                  typename Dune::GDT::
                      ExplicitRungeKuttaTimeStepper<OperatorImp, DiscreteFunctionImp, TimeFieldImp, method>>::type
          TimeStepperType;

  static TimeStepperType create(const OperatorImp& op,
                                const DiscreteFunctionImp& initial_values,
                                const typename DiscreteFunctionImp::RangeFieldType r = 1.0,
                                const TimeFieldImp t_0 = 0.0,
                                const typename DiscreteFunctionImp::RangeFieldType tol = 1e-4)
  {
    return TimeStepperType(op, initial_values, r, t_0, tol);
  }
}; // class TimeStepperFactory


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_FACTORY_HH
