// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first

#include "operators/weighted-l2.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;

typedef testing::Types< SPACE_FV_SGRID(1, 1)
                      , SPACE_FV_SGRID(2, 1)
                      , SPACE_FV_SGRID(3, 1)
                      > ConstantSpaces;

TYPED_TEST_CASE(WeightedL2OperatorTest, ConstantSpaces);
TYPED_TEST(WeightedL2OperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor();
}
TYPED_TEST(WeightedL2OperatorTest, constructible_by_factory) {
  this->constructible_by_factory();
}
TYPED_TEST(WeightedL2OperatorTest, apply_is_callable) {
  this->apply_is_callable();
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_constant_arguments) {
  this->correct_for_constant_arguments();
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_linear_arguments) {
  this->correct_for_linear_arguments();
}
TYPED_TEST(WeightedL2OperatorTest, apply2_correct_for_quadratic_arguments) {
  this->correct_for_quadratic_arguments();
}
