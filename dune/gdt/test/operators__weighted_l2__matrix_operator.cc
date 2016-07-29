// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include <dune/stuff/test/main.hxx> // <- this one has to come first

#include "operators/weighted-l2.hh"
#include "spaces/fv/default.hh"
#include "spaces/dg/fem.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;


#if HAVE_DUNE_FEM

typedef testing::Types<SPACE_DG_FEM_SGRID(1, 1, 2), SPACE_DG_FEM_SGRID(2, 1, 2), SPACE_DG_FEM_SGRID(3, 1, 2)>
    QuadraticSpaces;
TYPED_TEST_CASE(WeightedL2MatrixOperatorTest, QuadraticSpaces);

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types<SPACE_CG_PDELAB_SGRID(1, 1, 1), SPACE_CG_PDELAB_SGRID(2, 1, 1), SPACE_CG_PDELAB_SGRID(3, 1, 1)>
    LinearSpaces;
TYPED_TEST_CASE(WeightedL2MatrixOperatorTest, LinearSpaces);

#else // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB

typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1)> ConstantSpaces;
TYPED_TEST_CASE(WeightedL2MatrixOperatorTest, ConstantSpaces);

#endif // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB


TYPED_TEST(WeightedL2MatrixOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(WeightedL2MatrixOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(WeightedL2MatrixOperatorTest, is_matrix_operator)
{
  this->is_matrix_operator();
}
TYPED_TEST(WeightedL2MatrixOperatorTest, correct_for_constant_arguments)
{
#ifndef NDEBUG
  const double tolerance = 2.85e-14;
#else
  const double tolerance = 4.27e-14;
#endif
  this->correct_for_constant_arguments(this->dimDomain == 1 ? 1.43e-14 : (this->dimDomain == 2 ? 7.11e-15 : tolerance));
}

#if HAVE_DUNE_FEM || HAVE_DUNE_PDELAB
TYPED_TEST(WeightedL2MatrixOperatorTest, correct_for_linear_arguments)
{
#ifndef NDEBUG
  const double tolerance = 1e-15;
#else
  const double tolerance = 8.89e-15;
#endif
  this->correct_for_linear_arguments(this->dimDomain == 3 ? 2.49e-14 : tolerance);
}
#else
TEST(DISABLED_WeightedL2MatrixOperatorTest, correct_for_linear_arguments)
{
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif

#if HAVE_DUNE_FEM
TYPED_TEST(WeightedL2MatrixOperatorTest, correct_for_quadratic_arguments)
{
#ifndef NDEBUG
  const double tolerance = 1e-15;
#else
  const double tolerance = 1.78e-15;
#endif
  this->correct_for_quadratic_arguments(this->dimDomain == 3 ? 5.33e-15 : tolerance);
}
#else
TEST(DISABLED_WeightedL2MatrixOperatorTest, correct_for_quadratic_arguments)
{
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif
