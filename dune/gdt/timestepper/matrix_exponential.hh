// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
#define DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH

#include "interface.hh"

namespace Dune {
namespace GDT {


/** \brief Time stepper solving linear equation d_t u = Au + b by matrix exponential
 */
template <class OperatorImp, class DiscreteFunctionImp, class TimeFieldImp = double>
class MatrixExponentialTimeStepper : public TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp> BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;
  using typename BaseType::DataHandleType;

  typedef OperatorImp OperatorType;

  using BaseType::current_solution;
  using BaseType::current_time;

  MatrixExponentialTimeStepper(const OperatorType& op,
                               const DiscreteFunctionType& initial_values,
                               const TimeFieldImp t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , op_(op)
  {
  }

  virtual TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt) override final
  {
    const TimeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();

    op_.apply_matrix_exponential(u_n, t, actual_dt);

    // augment time
    t += actual_dt;

    return dt;
  } // ... step(...)

private:
  const OperatorType& op_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
