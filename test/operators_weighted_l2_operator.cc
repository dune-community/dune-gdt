// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"
#include "spaces_cg_pdelab.hh"

#include <dune/gdt/tests/operators/weighted-l2.hh>

using namespace Dune::GDT::Tests;

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
