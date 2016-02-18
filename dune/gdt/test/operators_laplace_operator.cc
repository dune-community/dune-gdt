// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"

#include <dune/gdt/tests/operators/laplace.hh>

using namespace Dune::GDT::Tests;

typedef testing::Types< SPACE_FV_SGRID(1, 1)
                      , SPACE_FV_SGRID(2, 1)
                      , SPACE_FV_SGRID(3, 1)
                      > ConstantSpaces;

TYPED_TEST_CASE(LaplaceOperatorTest, ConstantSpaces);
TYPED_TEST(LaplaceOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor();
}
TYPED_TEST(LaplaceOperatorTest, constructible_by_factory) {
  this->constructible_by_factory();
}
TYPED_TEST(LaplaceOperatorTest, apply_is_callable) {
  this->apply_is_callable();
}
TYPED_TEST(LaplaceOperatorTest, apply2_correct_for_constant_arguments) {
  this->correct_for_constant_arguments();
}
TYPED_TEST(LaplaceOperatorTest, apply2_correct_for_linear_arguments) {
  this->correct_for_linear_arguments();
}
TYPED_TEST(LaplaceOperatorTest, apply2_correct_for_quadratic_arguments) {
  this->correct_for_quadratic_arguments();
}
