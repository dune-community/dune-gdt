// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_TIMESTEPPER_ENUMS_HH
#define DUNE_GDT_TIMESTEPPER_ENUMS_HH

namespace Dune {
namespace GDT {


enum class TimeStepperMethods
{
  euler_heun,
  bogacki_shampine,
  dormand_prince,
  adaptive_rungekutta_other,
  explicit_euler,
  explicit_rungekutta_second_order_ssp,
  explicit_rungekutta_third_order_ssp,
  explicit_rungekutta_classic_fourth_order,
  explicit_rungekutta_other,
  implicit_euler,
  implicit_midpoint,
  trapezoidal_rule,
  diagonally_implicit_other,
  matrix_exponential
};

enum class TimeStepperSplittingMethods
{
  fractional_step,
  strang
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_ENUMS_HH
