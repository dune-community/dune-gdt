// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH
#define DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH

#include <dune/gdt/projections/lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


struct LagrangeProjectionLocalizableOperatorTest
    : public LocalizableProjectionOperatorBase<SPACETYPE,
                                               LagrangeProjectionLocalizableOperator<
                                                   typename SPACETYPE::GridViewType,
                                                   typename internal::OperatorBaseTraits<SPACETYPE>::FunctionType,
                                                   typename internal::OperatorBaseTraits<SPACETYPE>::
                                                       DiscreteFunctionType>>
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    auto& range           = this->discrete_function_;

    auto with_grid_view DUNE_UNUSED = make_lagrange_projection_localizable_operator(grid_view, source, range);
    auto wo_grid_view DUNE_UNUSED   = make_lagrange_projection_localizable_operator(source, range);
  } // ... constructible_by_factory(...)
}; // struct LagrangeProjectionLocalizableOperatorTest


struct LagrangeProjectionOperatorTest
    : public ProjectionOperatorBase<SPACETYPE, LagrangeProjectionOperator<typename SPACETYPE::GridViewType, double>>
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    auto op DUNE_UNUSED   = make_lagrange_projection_operator(grid_view);
  } // ... constructible_by_factory(...)

  void free_function_callable()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    auto& range           = this->discrete_function_;

    Dune::GDT::project_lagrange(grid_view, source, range);
    Dune::GDT::project_lagrange(source, range);
  } // ... free_function_callable(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_LAGRANGE_HH
