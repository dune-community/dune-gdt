// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk (2017 - 2018)
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

typedef ProjectionTest<{{SpaceType}}> ProjectionTest_{{Name}};
TEST_F(ProjectionTest_{{Name}}, produces_correct_results)
{
  {% if 'FvSpace' in SpaceType %}
    const double tolerance = 0.096226;
  {% elif 'RaviartThomasSpace' in SpaceType %}
    const double tolerance = 0.0925927;
  {% else %}
    const auto tolerance = this->default_tolerance;
  {% endif %}
  this->produces_correct_results(tolerance);
}
{% endfor %}
// clang-format on
