// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first

#include "elliptic.hh"
#include <dune/gdt/test/spaces/fv/default.hh>

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(EllipticLocalizableProductTest, ConstantSpaces);
TYPED_TEST(EllipticLocalizableProductTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(EllipticLocalizableProductTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(EllipticLocalizableProductTest, is_localizable_product)
{
  this->is_localizable_product();
}
TYPED_TEST(EllipticLocalizableProductTest, correct_for_constant_arguments)
{
  this->correct_for_constant_arguments(this->dimDomain == 3 ? 4.27e-14 : 1e-15);
}
TYPED_TEST(EllipticLocalizableProductTest, correct_for_linear_arguments)
{
  this->correct_for_linear_arguments(this->dimDomain == 1 ? 3.56e-15 : (this->dimDomain == 2 ? 7.11e-15 : 1.43e-14));
}
TYPED_TEST(EllipticLocalizableProductTest, correct_for_quadratic_arguments)
{
#ifndef NDEBUG
  const double tolerance = 1e-15;
#else
  const double tolerance = 1.78e-15;
#endif
  this->correct_for_quadratic_arguments(this->dimDomain == 1 ? tolerance : 3.56e-15);
}
