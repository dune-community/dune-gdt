// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"
#include "spaces_cg_pdelab.hh"

#include <dune/gdt/tests/operators/elliptic.hh>

using namespace Dune::GDT::Tests;


#if HAVE_DUNE_FEM

typedef testing::Types< SPACE_DG_FEM_SGRID(1, 1, 3)
                      , SPACE_DG_FEM_SGRID(2, 1, 3)
                      , SPACE_DG_FEM_SGRID(3, 1, 3)
                      > CubicSpaces;
TYPED_TEST_CASE(EllipticMatrixOperatorTest, CubicSpaces);

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types< SPACE_CG_PDELAB_SGRID(1, 1, 1)
                      , SPACE_CG_PDELAB_SGRID(2, 1, 1)
                      , SPACE_CG_PDELAB_SGRID(3, 1, 1)
                      > LinearSpaces;
TYPED_TEST_CASE(EllipticMatrixOperatorTest, LinearSpaces);

#else // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB

typedef testing::Types< SPACE_FV_SGRID(1, 1)
                      , SPACE_FV_SGRID(2, 1)
                      , SPACE_FV_SGRID(3, 1)
                      > ConstantSpaces;
TYPED_TEST_CASE(EllipticMatrixOperatorTest, ConstantSpaces);

#endif // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB


TYPED_TEST(EllipticMatrixOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor();
}
TYPED_TEST(EllipticMatrixOperatorTest, constructible_by_factory) {
  this->constructible_by_factory();
}
TYPED_TEST(EllipticMatrixOperatorTest, is_matrix_operator) {
  this->is_matrix_operator();
}

#if HAVE_DUNE_FEM || HAVE_DUNE_PDELAB
TYPED_TEST(EllipticMatrixOperatorTest, correct_for_constant_arguments) {
  this->correct_for_constant_arguments();
}
#else
TEST(DISABLED_EllipticMatrixOperatorTest, correct_for_constant_arguments) {
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif

#if HAVE_DUNE_FEM
TYPED_TEST(EllipticMatrixOperatorTest, correct_for_linear_arguments) {
  this->correct_for_linear_arguments();
}
TYPED_TEST(EllipticMatrixOperatorTest, correct_for_quadratic_arguments) {
  this->correct_for_quadratic_arguments();
}
#else // HAVE_DUNE_FEM
TEST(DISABLED_EllipticMatrixOperatorTest, correct_for_linear_arguments) {
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
TEST(DISABLED_EllipticMatrixOperatorTest, correct_for_quadratic_arguments) {
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif // HAVE_DUNE_FEM
