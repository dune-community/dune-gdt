// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include <dune/gdt/test/spaces/cg.hh>
#include <dune/gdt/test/spaces/dg.hh>
#include <dune/gdt/test/spaces/fv.hh>
#include <dune/gdt/test/spaces/rt.hh>

// clang-format off
{% for SpaceType,Name in config.spaces_with_names %}

{% if 'DunePdelabRtSpaceWrapper' in SpaceType %}
  typedef RT_Space<{{SpaceType}}> TestType_{{Name}};
{% else %}
  typedef SpaceBase<{{SpaceType}}>
    TestType_{{Name}};
{% endif %}

TEST_F(TestType_{{Name}}, fulfills_interface)
{
  this->fulfills_interface();
}
TEST_F(TestType_{{Name}}, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}
TEST_F(TestType_{{Name}}, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}
TEST_F(TestType_{{Name}}, check_for_correct_copy)
{
  this->check_for_correct_copy();
}

{% if 'DunePdelabRtSpaceWrapper' in SpaceType %}
  TEST_F(TestType_{{Name}}, matches_raviart_thomas_signature)
  {
    this->matches_raviart_thomas_signature();
  }

  {% if 'simplex' in Name and '::conforming' in Name %}
    TEST_F(RT_Space_{{Name}}, fulfills_raviart_thomas_2d_simplicial_interface)
    {
      this->fulfills_raviart_thomas_2d_simplicial_interface();
    }
  {% endif %}

{% endif %}

{% if 'CgSpaceWrapper' in SpaceType %}
  typedef P1Q1_CG_Space<{{SpaceType}}> P1Q1_CG_Space_{{Name}};
  TEST_F(P1Q1_CG_Space_{{Name}}, fulfills_continuous_interface)
  {
    this->fulfills_continuous_interface();
  }
  TEST_F(P1Q1_CG_Space_{{Name}}, maps_correctly)
  {
    this->maps_correctly();
  }
{% endif %}

{% if 'DgSpaceWrapper' in SpaceType %}
  typedef P1Q1_DG_Space<{{SpaceType}}> P1Q1_DG_Space_{{Name}};
  TEST_F(P1Q1_DG_Space_{{Name}}, fulfills_continuous_interface)
  {
    this->maps_correctly();
  }
{% endif %}


{% endfor %}
// clang-format on
