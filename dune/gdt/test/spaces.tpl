// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt/default.hh>

#include <dune/gdt/test/grids.hh>

// clang-format off
{% for SpaceType,Name in config.spaces_with_names %}

typedef SpaceBase<{{SpaceType}}> TestType_{{Name}};

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

{% endfor %}
// clang-format on
