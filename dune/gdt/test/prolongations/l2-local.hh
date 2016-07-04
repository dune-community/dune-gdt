// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH
#define DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH

#include <dune/common/unused.hh>

#include <dune/gdt/prolongations/l2-local.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct L2LocalProlongationLocalizableOperatorTest
    : public LocalizableProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationLocalizableOperator>
{
  typedef LocalizableProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationLocalizableOperator> BaseType;
  using typename BaseType::ProlongationOperatorType;

  void constructible_by_ctor(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_view     = this->fine_space_.grid_view();
    const auto& source = this->coarse_discrete_function_;
    auto& range        = this->fine_discrete_function_;

    ProlongationOperatorType DUNE_UNUSED(w_over_integrate)(0, grid_view, source, range);
    ProlongationOperatorType DUNE_UNUSED(wo_over_integrate)(grid_view, source, range);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_view     = this->fine_space_.grid_view();
    const auto& source = this->coarse_discrete_function_;
    auto& range        = this->fine_discrete_function_;

    auto DUNE_UNUSED(w_gv_w_oi)   = make_local_l2_prolongation_localizable_operator(grid_view, source, range, 1);
    auto DUNE_UNUSED(w_gv_wo_oi)  = make_local_l2_prolongation_localizable_operator(grid_view, source, range);
    auto DUNE_UNUSED(wo_gv_w_oi)  = make_local_l2_prolongation_localizable_operator(source, range, 1);
    auto DUNE_UNUSED(wo_gv_wo_oi) = make_local_l2_prolongation_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
};


template <class SpaceType>
struct L2LocalProlongationOperatorTest
    : public ProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationOperator>
{
  typedef ProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationOperator> BaseType;
  using typename BaseType::ProlongationOperatorType;

  void constructible_by_ctor(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_view     = this->fine_space_.grid_view();
    const auto& source = this->coarse_discrete_function_;
    auto& range        = this->fine_discrete_function_;

    ProlongationOperatorType DUNE_UNUSED(w_over_integrate)(0, grid_view);
    ProlongationOperatorType DUNE_UNUSED(wo_over_integrate)(grid_view);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_view = this->fine_space_.grid_view();

    auto DUNE_UNUSED(w_over_integrate)  = make_local_l2_prolongation_operator(grid_view, 1);
    auto DUNE_UNUSED(wo_over_integrate) = make_local_l2_prolongation_operator(grid_view);
  } // ... constructible_by_factory(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH
