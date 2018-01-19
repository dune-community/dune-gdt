// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2018)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/prolongations/l2-local.hh>
#include <dune/gdt/test/prolongations/l2.hh>
#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-functions-wrapper.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

using namespace Dune::GDT::Test;

// clang-format off
{% for SpaceType,Name in config.spaces_with_names %}

typedef L2LocalProlongationOperatorTest<{{SpaceType}}>
  L2LocalProlongationOperatorTest_{{Name}};
typedef L2LocalProlongationLocalizableOperatorTest<{{SpaceType}}>
  L2LocalProlongationLocalizableOperatorTest_{{Name}};

{% if 'FvSpace' in SpaceType %}
  const double {{Name}}_tolerance = 1.45e-1;
{% elif 'DunePdelabRtSpaceWrapper' in SpaceType %}
    const auto {{Name}}_tolerance = pdelab_rt_tolerance<L2LocalProlongationOperatorTest_{{Name}}>();
{% elif 'FemCg' in SpaceType %}
    const auto {{Name}}_tolerance = fem_cg_tolerance<L2LocalProlongationOperatorTest_{{Name}}>();
{% else %}
  const auto {{Name}}_tolerance = Dune::XT::Grid::is_alugrid<Dune::XT::Grid::extract_grid_t<{{SpaceType}}::GridLayerType>>::value
      ? L2LocalProlongationOperatorTest_{{Name}}::alugrid_tolerance
      : L2LocalProlongationOperatorTest_{{Name}}::default_tolerance;
{% endif %}

TEST_F(L2LocalProlongationOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2LocalProlongationOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory({{Name}}_tolerance);
}
TEST_F(L2LocalProlongationOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results({{Name}}_tolerance);
}

TEST_F(L2LocalProlongationLocalizableOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2LocalProlongationLocalizableOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory({{Name}}_tolerance);
}
TEST_F(L2LocalProlongationLocalizableOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results({{Name}}_tolerance);
  this->produces_correct_results({{Name}}_tolerance);
}


{% endfor %}
// clang-format on
