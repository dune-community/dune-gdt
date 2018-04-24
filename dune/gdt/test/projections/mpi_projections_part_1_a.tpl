// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)
//
// reserved.

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/projections.hh>

#include <dune/gdt/test/projections/lagrange.hh>
#include <dune/gdt/test/projections/l2-global.hh>


#include <dune/gdt/spaces/cg.hh>

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt/default.hh>

#include <dune/gdt/test/projections/projections.hh>
#include <dune/gdt/test/projections/l2.hh>

using namespace Dune::GDT::Test;

// clang-format off
{% for SpaceType,Name in config.spaces_with_names %}


typedef L2ProjectionOperatorTest<{{SpaceType}}> L2ProjectionOperatorTest_{{Name}};
TEST_F(L2ProjectionOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2ProjectionOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2ProjectionOperatorTest_{{Name}}, free_function_callable)
{
  this->free_function_callable();
}
TEST_F(L2ProjectionOperatorTest_{{Name}}, produces_correct_results)
{
  // RT : 0.0925927
  {% if 'FvSpace' in SpaceType %}
    const double tolerance = 0.096226;
  {% elif 'RaviartThomasSpace' in SpaceType %}
    const double tolerance = 0.0925927;
  {% else %}
    using Grid = Dune::XT::Grid::extract_grid_t<typename L2ProjectionOperatorTest_{{Name}}::GridLayerType>;
    const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  {% endif %}
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}
{% endfor %}
// clang-format on
