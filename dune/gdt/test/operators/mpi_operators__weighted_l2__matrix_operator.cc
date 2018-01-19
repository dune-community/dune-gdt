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

#include "weighted-l2.hh"
#include <dune/gdt/test/spaces/fv/default.hh>
#include <dune/gdt/test/spaces/dg/fem.hh>
#include <dune/gdt/test/spaces/cg/pdelab.hh>

using namespace Dune::GDT::Test;


#if HAVE_DUNE_FEM

typedef testing::Types<SPACE_DG_FEM_YASPGRID(1, 1, 2), SPACE_DG_FEM_YASPGRID(2, 1, 2), SPACE_DG_FEM_YASPGRID(3, 1, 2)>
    QuadraticSpaces;
TYPED_TEST_CASE(WeightedL2MatrixOperatorTest, QuadraticSpaces);

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::
    Types<SPACE_CG_PDELAB_YASPGRID(1, 1, 1), SPACE_CG_PDELAB_YASPGRID(2, 1, 1), SPACE_CG_PDELAB_YASPGRID(3, 1, 1)>
        LinearSpaces;
TYPED_TEST_CASE(WeightedL2MatrixOperatorTest, LinearSpaces);

#else // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB

typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)> ConstantSpaces;
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
  const double tolerance = 8.53e-14;
#endif
  this->correct_for_constant_arguments(this->dimDomain == 1 ? 2.2e-14 : (this->dimDomain == 2 ? 2.85e-14 : tolerance));
}

#if HAVE_DUNE_FEM || HAVE_DUNE_PDELAB
TYPED_TEST(WeightedL2MatrixOperatorTest, correct_for_linear_arguments)
{
#ifndef NDEBUG
  const double tolerance = 1e-15;
#else
  const double tolerance = 1.25e-14;
#endif
  this->correct_for_linear_arguments(this->dimDomain == 3 ? 2.67e-14 : tolerance);
}
#else
TEST(DISABLED_WeightedL2MatrixOperatorTest, correct_for_linear_arguments)
{
  std::cerr << Dune::XT::Common::colorStringRed("Missing dependencies!") << std::endl;
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
  this->correct_for_quadratic_arguments(this->dimDomain == 3 ? 1.43e-14 : tolerance);
}
#else
TEST(DISABLED_WeightedL2MatrixOperatorTest, correct_for_quadratic_arguments)
{
  std::cerr << Dune::XT::Common::colorStringRed("Missing dependencies!") << std::endl;
}
#endif
