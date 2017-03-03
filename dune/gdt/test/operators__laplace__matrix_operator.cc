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

#include "operators/laplace.hh"
#include "spaces/fv/default.hh"
#include "spaces/dg/fem.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;


#if HAVE_DUNE_FEM

typedef testing::Types<SPACE_DG_FEM_YASPGRID(1, 1, 3), SPACE_DG_FEM_YASPGRID(2, 1, 3), SPACE_DG_FEM_YASPGRID(3, 1, 3)>
    CubicSpaces;
TYPED_TEST_CASE(LaplaceMatrixOperatorTest, CubicSpaces);

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types<SPACE_CG_PDELAB_YASPGRID(1, 1, 1),
                       SPACE_CG_PDELAB_YASPGRID(2, 1, 1),
                       SPACE_CG_PDELAB_YASPGRID(3, 1, 1)>
    LinearSpaces;
TYPED_TEST_CASE(LaplaceMatrixOperatorTest, LinearSpaces);

#else // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB

typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)> ConstantSpaces;
TYPED_TEST_CASE(LaplaceMatrixOperatorTest, ConstantSpaces);

#endif // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB


TYPED_TEST(LaplaceMatrixOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LaplaceMatrixOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LaplaceMatrixOperatorTest, is_matrix_operator)
{
  this->is_matrix_operator();
}

#if HAVE_DUNE_FEM || HAVE_DUNE_PDELAB
TYPED_TEST(LaplaceMatrixOperatorTest, correct_for_constant_arguments)
{
  this->correct_for_constant_arguments();
}
#else
TEST(DISABLED_LaplaceMatrixOperatorTest, correct_for_constant_arguments)
{
  std::cerr << Dune::XT::Common::colorStringRed("Missing dependencies!") << std::endl;
}
#endif

#if HAVE_DUNE_FEM
TYPED_TEST(LaplaceMatrixOperatorTest, correct_for_linear_arguments)
{
  this->correct_for_linear_arguments();
}
TYPED_TEST(LaplaceMatrixOperatorTest, correct_for_quadratic_arguments)
{
  this->correct_for_quadratic_arguments();
}
#else // HAVE_DUNE_FEM
TEST(DISABLED_LaplaceMatrixOperatorTest, correct_for_linear_arguments)
{
  std::cerr << Dune::XT::Common::colorStringRed("Missing dependencies!") << std::endl;
}
TEST(DISABLED_LaplaceMatrixOperatorTest, correct_for_quadratic_arguments)
{
  std::cerr << Dune::XT::Common::colorStringRed("Missing dependencies!") << std::endl;
}
#endif // HAVE_DUNE_FEM
