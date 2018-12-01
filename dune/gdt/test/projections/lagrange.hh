// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH
#define DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH

#include <dune/gdt/projections/lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {

template <class SpaceType>
struct LagrangeProjectionLocalizableOperatorTest
  : public LocalizableProjectionOperatorBase<
        SpaceType,
        LagrangeProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                              typename internal::OperatorBaseTraits<SpaceType>::FunctionType,
                                              typename internal::OperatorBaseTraits<SpaceType>::DiscreteFunctionType>>
{
  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    auto& range = this->discrete_function_;

    auto with_grid_layer DUNE_UNUSED = make_lagrange_projection_localizable_operator(grid_layer, source, range);
    auto wo_grid_layer DUNE_UNUSED = make_lagrange_projection_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
}; // struct LagrangeProjectionLocalizableOperatorTest

template <class SpaceType>
struct LagrangeProjectionOperatorTest
  : public ProjectionOperatorBase<SpaceType, LagrangeProjectionOperator<typename SpaceType::GridLayerType, double>>
{
  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();
    auto op DUNE_UNUSED = make_lagrange_projection_operator(grid_layer);
  } // ... constructible_by_factory(...)

  void free_function_callable()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    auto& range = this->discrete_function_;

    Dune::GDT::project_lagrange(grid_layer, source, range);
    Dune::GDT::project_lagrange(source, range);
  } // ... free_function_callable(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH
