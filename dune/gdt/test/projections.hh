#ifndef DUNE_GDT_TEST_PROJECTIONS_HH
#define DUNE_GDT_TEST_PROJECTIONS_HH

#include <dune/gdt/projections.hh>
#include <dune/gdt/test/projections/base.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct ProjectionTest : public internal::ProjectionOperatorBase<SpaceType>
{
  void produces_correct_results(const double& tolerance = 1e-15)
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->function_;
    auto& range           = this->discrete_function_;

    project(grid_view, source, range);
    project(source, range);

    this->measure_error(tolerance);
  } // ... produces_correct_results(...)
}; // struct ProjectionTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_HH
