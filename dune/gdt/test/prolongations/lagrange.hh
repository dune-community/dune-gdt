// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

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

    auto grid_view     = this->fine_space_.grid_view();
    const auto& source = this->coarse_discrete_function_;
    auto& range        = this->fine_discrete_function_;

    ProlongationOperatorType DUNE_UNUSED(op)(grid_view, source, range);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory(const double tolerance = 1e-15)
  {
    this->prepare(tolerance);

    auto grid_view     = this->fine_space_.grid_view();
    const auto& source = this->coarse_discrete_function_;
    auto& range        = this->fine_discrete_function_;

    auto DUNE_UNUSED(w_gv) = make_lagrange_prolongation_localizable_operator(grid_view, source, range);
    auto DUNE_UNUSED(wo_gv) = make_lagrange_prolongation_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_LAGRANGE_HH
