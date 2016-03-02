// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PROJECTIONS_L2_GLOBAL_HH
#define DUNE_GDT_TEST_PROJECTIONS_L2_GLOBAL_HH

#include <dune/gdt/projections/l2-global.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template< class SpaceType >
struct L2GlobalProjectionLocalizableOperatorTest
  : public LocalizableProjectionOperatorBase< SpaceType, L2GlobalProjectionLocalizableOperator<
        typename SpaceType::GridViewType,
        typename internal::OperatorBaseTraits< SpaceType >::FunctionType,
        typename internal::OperatorBaseTraits< SpaceType >::DiscreteFunctionType > >
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    auto DUNE_UNUSED(w_grid_view_w_over_integrate)
        = make_global_l2_projection_localizable_operator(grid_view, source, range, 1);
    auto DUNE_UNUSED(w_grid_view_wo_over_integrate)
        = make_global_l2_projection_localizable_operator(grid_view, source, range);
    auto DUNE_UNUSED(wo_grid_view_w_over_integrate)
        = make_global_l2_projection_localizable_operator(           source, range, 1);
    auto DUNE_UNUSED(wo_grid_view_wo_over_integrate)
        = make_global_l2_projection_localizable_operator(           source, range);
  } // ... constructible_by_factory(...)
};


template< class SpaceType >
struct L2GlobalProjectionOperatorTest
  : public ProjectionOperatorBase< SpaceType
                                 , L2GlobalProjectionOperator< typename SpaceType::GridViewType, double > >
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();

    auto DUNE_UNUSED(w_over_integrate)  = make_global_l2_projection_operator(grid_view, 1);
    auto DUNE_UNUSED(wo_over_integrate) = make_global_l2_projection_operator(grid_view);
  } // ... constructible_by_factory(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_L2_GLOBAL_HH
