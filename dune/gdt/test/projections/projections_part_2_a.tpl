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

typedef L2ProjectionLocalizableOperatorTest<{{SpaceType}}> L2ProjectionLocalizableOperatorTest_{{Name}};
TEST_F(L2ProjectionLocalizableOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2ProjectionLocalizableOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2ProjectionLocalizableOperatorTest_{{Name}}, produces_correct_results)
{
  // RT : 0.096226
  {% if 'FvSpace' in SpaceType %}
    const double tolerance = 0.096226;
  {% elif 'DunePdelabRtSpaceWrapper' in SpaceType %}
    const double tolerance = 0.0925927;
  {% else %}
    typedef Dune::XT::Grid::extract_grid_t<L2ProjectionLocalizableOperatorTest_{{Name}}::GridLayerType> Grid;
    const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  {% endif %}
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}
{% endfor %}
// clang-format on
