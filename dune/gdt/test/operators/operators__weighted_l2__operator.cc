// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first

#include "weighted-l2.hh"
#include <dune/gdt/test/spaces/fv/default.hh>

using namespace Dune::GDT::Test;

typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(WeightedL2OperatorTest, ConstantSpaces);
TYPED_TEST(WeightedL2OperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(WeightedL2OperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(WeightedL2OperatorTest, apply_is_callable)
{
  this->apply_is_callable();
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_constant_arguments)
{
  this->correct_for_constant_arguments(this->dimDomain == 3 ? 2.14e-14 : 3.5531e-15);
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_linear_arguments)
{
  this->correct_for_linear_arguments(this->dimDomain == 3 ? 7.11e-15 : 3.553e-15);
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_quadratic_arguments)
{
  this->correct_for_quadratic_arguments(3.553e-15);
}
