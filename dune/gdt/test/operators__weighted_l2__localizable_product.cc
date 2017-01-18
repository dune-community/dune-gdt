// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first

#include "operators/weighted-l2.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(WeightedL2LocalizableProductTest, ConstantSpaces);
TYPED_TEST(WeightedL2LocalizableProductTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(WeightedL2LocalizableProductTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(WeightedL2LocalizableProductTest, is_localizable_product)
{
  this->is_localizable_product();
}
TYPED_TEST(WeightedL2LocalizableProductTest, correct_for_constant_arguments)
{
  this->correct_for_constant_arguments(this->dimDomain == 3 ? 2.14e-14 : 1e-15);
}
TYPED_TEST(WeightedL2LocalizableProductTest, correct_for_linear_arguments)
{
  this->correct_for_linear_arguments(this->dimDomain == 3 ? 5.33e-15 : 3.57e-15);
}
TYPED_TEST(WeightedL2LocalizableProductTest, correct_for_quadratic_arguments)
{
  this->correct_for_quadratic_arguments(3.57e-15);
}
