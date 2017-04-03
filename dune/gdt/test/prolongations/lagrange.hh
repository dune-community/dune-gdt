// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_LAGRANGE_HH
#define DUNE_GDT_TEST_PROLONGATIONS_LAGRANGE_HH

#include <dune/common/unused.hh>

#include <dune/gdt/prolongations/lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct LagrangeProlongationLocalizableOperatorTest
    : public LocalizableProlongationOperatorBase<SpaceType, SpaceType, LagrangeProlongationLocalizableOperator>
{
  typedef LocalizableProlongationOperatorBase<SpaceType, SpaceType, LagrangeProlongationLocalizableOperator> BaseType;
  using typename BaseType::ProlongationOperatorType;

  void constructible_by_ctor(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_layer = this->fine_space_.grid_layer();
    const auto& source = this->coarse_discrete_function_;
    auto& range = this->fine_discrete_function_;

    DUNE_UNUSED ProlongationOperatorType op(grid_layer, source, range);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_layer = this->fine_space_.grid_layer();
    const auto& source = this->coarse_discrete_function_;
    auto& range = this->fine_discrete_function_;

    auto w_gv DUNE_UNUSED = make_lagrange_prolongation_localizable_operator(grid_layer, source, range);
    auto wo_gv DUNE_UNUSED = make_lagrange_prolongation_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
};


template <class SpaceType>
struct LagrangeProlongationOperatorTest
    : public ProlongationOperatorBase<SpaceType, SpaceType, LagrangeProlongationOperator>
{
  typedef ProlongationOperatorBase<SpaceType, SpaceType, LagrangeProlongationOperator> BaseType;
  using typename BaseType::ProlongationOperatorType;

  void constructible_by_ctor(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_layer = this->fine_space_.grid_layer();

    DUNE_UNUSED ProlongationOperatorType op(grid_layer);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_layer = this->fine_space_.grid_layer();

    auto op DUNE_UNUSED = make_lagrange_prolongation_operator(grid_layer);
  } // ... constructible_by_factory(...)

  void free_function_callable(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_layer = this->fine_space_.grid_layer();
    const auto& source = this->coarse_discrete_function_;
    auto& range = this->fine_discrete_function_;

    prolong_lagrange(grid_layer, source, range);
    prolong_lagrange(source, range);
  } // ... free_function_callable(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_LAGRANGE_HH
