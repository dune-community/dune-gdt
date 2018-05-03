// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   Rene Milk       (2016 - 2018)

#ifndef DUNE_GDT_TEST_PROJECTIONS_HH
#define DUNE_GDT_TEST_PROJECTIONS_HH

#include <dune/gdt/projections.hh>
#include <dune/gdt/test/projections/base.hh>
#include <dune/gdt/test/spaces/cg.hh>

namespace Dune {
namespace GDT {
namespace Test {

template <class SPACETYPE>
struct ProjectionTest : public internal::ProjectionOperatorBase<SPACETYPE>
{
  void produces_correct_results(const double& tolerance = 1e-15)
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->function_;
    auto& range = this->discrete_function_;

    project(grid_layer, source, range);
    project(source, range);

    this->measure_error(tolerance);
  } // ... produces_correct_results(...)
}; // struct ProjectionTest

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_HH
