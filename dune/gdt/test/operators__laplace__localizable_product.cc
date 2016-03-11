// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first

#include "operators/laplace.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(LaplaceLocalizableProductTest, ConstantSpaces);
TYPED_TEST(LaplaceLocalizableProductTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LaplaceLocalizableProductTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LaplaceLocalizableProductTest, is_localizable_product)
{
  this->is_localizable_product();
}
TYPED_TEST(LaplaceLocalizableProductTest, correct_for_constant_arguments)
{
  this->correct_for_constant_arguments();
}
TYPED_TEST(LaplaceLocalizableProductTest, correct_for_linear_arguments)
{
  this->correct_for_linear_arguments();
}
TYPED_TEST(LaplaceLocalizableProductTest, correct_for_quadratic_arguments)
{
  this->correct_for_quadratic_arguments();
}
