// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_PROJECTIONS_L2_LOCAL_HH
#define DUNE_GDT_TEST_PROJECTIONS_L2_LOCAL_HH

#include <dune/gdt/projections/l2-local.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct L2LocalProjectionLocalizableOperatorTest
    : public LocalizableProjectionOperatorBase<SpaceType,
                                               L2LocalProjectionLocalizableOperator<
                                                   typename SpaceType::GridViewType,
                                                   typename internal::OperatorBaseTraits<SpaceType>::FunctionType,
                                                   typename internal::OperatorBaseTraits<SpaceType>::
                                                       DiscreteFunctionType>>
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    auto w_grid_view_w_over_integrate DUNE_UNUSED =
        make_local_l2_projection_localizable_operator(grid_view, source, range, 1);
    auto w_grid_view_wo_over_integrate DUNE_UNUSED =
        make_local_l2_projection_localizable_operator(grid_view, source, range);
    auto wo_grid_view_w_over_integrate DUNE_UNUSED = make_local_l2_projection_localizable_operator(source, range, 1);
    auto wo_grid_view_wo_over_integrate DUNE_UNUSED = make_local_l2_projection_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
};


template <class SpaceType>
struct L2LocalProjectionOperatorTest
    : public ProjectionOperatorBase<SpaceType, L2LocalProjectionOperator<typename SpaceType::GridViewType, double>>
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();

    auto w_over_integrate DUNE_UNUSED = make_local_l2_projection_operator(grid_view, 1);
    auto wo_over_integrate DUNE_UNUSED = make_local_l2_projection_operator(grid_view);
  } // ... constructible_by_factory(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_L2_LOCAL_HH
