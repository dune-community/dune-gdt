// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/prolongations/prolongations.hh>
#include <dune/gdt/test/prolongations/l2.hh>
#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

using namespace Dune::GDT::Test;

// clang-format off
{% for SpaceType,Name in config.spaces_with_names %}

typedef ProlongationTest<{{SpaceType}}>
  ProlongationTest_{{Name}};
typedef L2ProlongationOperatorTest<{{SpaceType}}>
  L2ProlongationOperatorTest_{{Name}};
typedef L2ProlongationLocalizableOperatorTest<{{SpaceType}}>
  L2ProlongationLocalizableOperatorTest_{{Name}};

{% if 'FvSpace' in SpaceType %}
  const double {{Name}}_tolerance = 1.45e-1;
{% elif 'DunePdelabRtSpaceWrapper' in SpaceType %}
    const auto {{Name}}_tolerance = pdelab_rt_tolerance<L2ProlongationOperatorTest_{{Name}}>();
{% elif 'FemCg' in SpaceType %}
    const auto {{Name}}_tolerance = fem_cg_tolerance<L2ProlongationOperatorTest_{{Name}}>();
{% else %}
  const auto {{Name}}_tolerance = Dune::XT::Grid::is_alugrid<Dune::XT::Grid::extract_grid_t<typename {{SpaceType}}::GridLayerType>>::value ? L2ProlongationOperatorTest_{{Name}}::alugrid_tolerance : L2ProlongationOperatorTest_{{Name}}::default_tolerance;
{% endif %}

TEST_F(ProlongationTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results({{Name}}_tolerance);
  this->produces_correct_results({{Name}}_tolerance);
}

TEST_F(L2ProlongationOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2ProlongationOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2ProlongationOperatorTest_{{Name}}, free_function_callable)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2ProlongationOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results({{Name}}_tolerance);
  this->produces_correct_results({{Name}}_tolerance);
}

TEST_F(L2ProlongationLocalizableOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor({{Name}}_tolerance);
}
TEST_F(L2ProlongationLocalizableOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory({{Name}}_tolerance);
}
TEST_F(L2ProlongationLocalizableOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results({{Name}}_tolerance);
  this->produces_correct_results({{Name}}_tolerance);
}

{% endfor %}
// clang-format on
