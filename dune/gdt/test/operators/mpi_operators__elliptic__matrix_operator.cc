// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first

#include <dune/gdt/test/spaces/dg.hh>

#include "elliptic.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_DG_YASPGRID(1, 1, 3), SPACE_DG_YASPGRID(2, 1, 3), SPACE_DG_YASPGRID(3, 1, 3)> CubicSpaces;
TYPED_TEST_CASE(EllipticMatrixOperatorTest, CubicSpaces);


TYPED_TEST(EllipticMatrixOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(EllipticMatrixOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(EllipticMatrixOperatorTest, is_matrix_operator)
{
  this->is_matrix_operator();
}

TYPED_TEST(EllipticMatrixOperatorTest, correct_for_constant_arguments)
{
  const double relax_factor = this->space_.grid_layer().grid().comm().size() > 1 ? 2 : 1;
  this->correct_for_quadratic_arguments(6.90e-13 * relax_factor);
}

TYPED_TEST(EllipticMatrixOperatorTest, correct_for_linear_arguments)
{
  const double relax_factor = this->space_.grid_layer().grid().comm().size() > 1 ? 2 : 1;
  this->correct_for_linear_arguments(EllipticDefaultTolerances::linear * relax_factor);
}

TYPED_TEST(EllipticMatrixOperatorTest, correct_for_quadratic_arguments)
{
  const double relax_factor = this->space_.grid_layer().grid().comm().size() > 1 ? 2 : 1;
  this->correct_for_quadratic_arguments(EllipticDefaultTolerances::quadratic * relax_factor);
}
