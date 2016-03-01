#ifndef DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH

#include <dune/gdt/operators/projections/l2.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template< class SpaceType >
struct L2ProjectionLocalizableOperatorTest
  : public LocalizableProjectionOperatorBase< SpaceType, L2ProjectionLocalizableOperator<
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
        = make_l2_projection_localizable_operator(grid_view, source, range, 1);
    auto DUNE_UNUSED(w_grid_view_wo_over_integrate)
        = make_l2_projection_localizable_operator(grid_view, source, range);
    auto DUNE_UNUSED(wo_grid_view_w_over_integrate)
        = make_l2_projection_localizable_operator(           source, range, 1);
    auto DUNE_UNUSED(wo_grid_view_wo_over_integrate)
        = make_l2_projection_localizable_operator(           source, range);
  } // ... constructible_by_factory(...)
};


template< class SpaceType >
struct L2ProjectionOperatorTest
  : public ProjectionOperatorBase< SpaceType
                                 , L2ProjectionOperator< typename SpaceType::GridViewType, double > >
{
  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    auto DUNE_UNUSED(op_w_over_integrate)  = make_l2_projection_operator(grid_view, 1);
    auto DUNE_UNUSED(op_wo_over_integrate) = make_l2_projection_operator(grid_view);
  } // ... constructible_by_factory(...)

  void free_function_callable()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    Dune::GDT::project_l2(grid_view, source, range);
    Dune::GDT::project_l2(source, range);
  } // ... free_function_callable(...)
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH
