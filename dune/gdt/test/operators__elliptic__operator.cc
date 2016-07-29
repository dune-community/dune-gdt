// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include <dune/stuff/test/main.hxx> // <- this one has to come first

#include "operators/elliptic.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;

typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(EllipticOperatorTest, ConstantSpaces);
TYPED_TEST(EllipticOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(EllipticOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(EllipticOperatorTest, apply_is_callable)
{
  this->apply_is_callable();
}
TYPED_TEST(EllipticOperatorTest, apply2_correct_for_constant_arguments)
{
  this->correct_for_constant_arguments(this->dimDomain == 3 ? 4.27e-14 : 1e-15);
}
TYPED_TEST(EllipticOperatorTest, apply2_correct_for_linear_arguments)
{
  this->correct_for_linear_arguments(this->dimDomain == 1 ? 3.56e-15 : (this->dimDomain == 2 ? 7.11e-15 : 1.43e-14));
}
TYPED_TEST(EllipticOperatorTest, apply2_correct_for_quadratic_arguments)
{
#ifndef NDEBUG
  const double tolerance = 1e-15;
#else
  const double tolerance = 1.77636e-15;
#endif
  this->correct_for_quadratic_arguments(this->dimDomain == 1 ? tolerance : 3.56e-15);
}
