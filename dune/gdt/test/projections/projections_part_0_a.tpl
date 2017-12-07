// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights
// reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/projections.hh>

#include <dune/gdt/test/projections/lagrange.hh>
#include <dune/gdt/test/projections/l2-global.hh>

#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-functions-wrapper.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

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
  {% elif 'DunePdelabRtSpaceWrapper' in SpaceType %}
    const double tolerance = 0.0925927;
  {% else %}
    const auto tolerance = this->default_tolerance;
  {% endif %}
  this->produces_correct_results(tolerance);
}
{% endfor %}
// clang-format on
