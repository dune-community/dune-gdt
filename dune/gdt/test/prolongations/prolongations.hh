// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_TEST_PROLONGATION_HH
#define DUNE_GDT_TEST_PROLONGATION_HH

#include <dune/gdt/prolongations.hh>
#include <dune/gdt/test/prolongations/base.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct ProlongationTest : public internal::ProlongationOperatorsBase<SpaceType, SpaceType>
{
  typedef internal::ProlongationOperatorsBase<SpaceType, SpaceType> BaseType;

  using BaseType::prepare;

  void produces_correct_results(const double& tolerance)
  {
    prepare(tolerance);

    prolong(fine_space_.grid_layer(), coarse_discrete_function_, fine_discrete_function_);
    auto fine_l2_error =
        make_l2_operator(fine_space_.grid_layer(), 2)->induced_norm(function_ - fine_discrete_function_);
    EXPECT_LE(fine_l2_error, tolerance);
    fine_discrete_function_.vector() *= 0.0;

    prolong(coarse_discrete_function_, fine_discrete_function_);
    fine_l2_error = make_l2_operator(fine_space_.grid_layer(), 2)->induced_norm(function_ - fine_discrete_function_);
    EXPECT_LE(fine_l2_error, tolerance);
  } // ... produces_correct_results(...)

  using BaseType::function_;
  using BaseType::fine_space_;
  using BaseType::coarse_discrete_function_;
  using BaseType::fine_discrete_function_;
}; // struct ProlongationTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATION_HH
