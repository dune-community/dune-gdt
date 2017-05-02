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

#ifndef DUNE_GDT_TEST_PROJECTIONS_L2_HH
#define DUNE_GDT_TEST_PROJECTIONS_L2_HH

#include <dune/gdt/projections/l2.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {

template <class SpaceType>
struct L2ProjectionLocalizableOperatorTest
    : public LocalizableProjectionOperatorBase<SpaceType,
                                               L2ProjectionLocalizableOperator<
                                                   typename SpaceType::GridLayerType,
                                                   typename internal::OperatorBaseTraits<SpaceType>::FunctionType,
                                                   typename internal::OperatorBaseTraits<SpaceType>::
                                                       DiscreteFunctionType>>
{
  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    auto w_grid_layer_w_over_integrate DUNE_UNUSED =
        make_l2_projection_localizable_operator(grid_layer, source, range, 1);
    auto w_grid_layer_wo_over_integrate DUNE_UNUSED =
        make_l2_projection_localizable_operator(grid_layer, source, range);
    auto wo_grid_layer_w_over_integrate DUNE_UNUSED = make_l2_projection_localizable_operator(source, range, 1);
    auto wo_grid_layer_wo_over_integrate DUNE_UNUSED = make_l2_projection_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
};

template <class SpaceType>
struct L2ProjectionOperatorTest
    : public ProjectionOperatorBase<SpaceType, L2ProjectionOperator<typename SpaceType::GridLayerType, double>>
{
  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();
    auto op_w_over_integrate DUNE_UNUSED = make_l2_projection_operator(grid_layer, 1);
    auto op_wo_over_integrate DUNE_UNUSED = make_l2_projection_operator(grid_layer);
  } // ... constructible_by_factory(...)

  void free_function_callable()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    project_l2(grid_layer, source, range);
    project_l2(source, range);
  } // ... free_function_callable(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_L2_HH
