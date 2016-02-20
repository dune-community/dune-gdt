// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"
#include "spaces_cg_pdelab.hh"

#include <dune/gdt/tests/operators/l2.hh>

using namespace Dune::GDT::Tests;


#if HAVE_DUNE_FEM

typedef testing::Types< SPACE_DG_FEM_SGRID(1, 1, 2)
                      , SPACE_DG_FEM_SGRID(2, 1, 2)
                      , SPACE_DG_FEM_SGRID(3, 1, 2)
                      > QuadraticSpaces;
TYPED_TEST_CASE(L2MatrixOperatorTest, QuadraticSpaces);

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types< SPACE_CG_PDELAB_SGRID(1, 1, 1)
                      , SPACE_CG_PDELAB_SGRID(2, 1, 1)
                      , SPACE_CG_PDELAB_SGRID(3, 1, 1)
                      > LinearSpaces;
TYPED_TEST_CASE(L2MatrixOperatorTest, LinearSpaces);

#else // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB

typedef testing::Types< SPACE_FV_SGRID(1, 1)
                      , SPACE_FV_SGRID(2, 1)
                      , SPACE_FV_SGRID(3, 1)
                      > ConstantSpaces;
TYPED_TEST_CASE(L2MatrixOperatorTest, ConstantSpaces);

#endif // HAVE_DUNE_FEM || HAVE_DUNE_PDELAB


TYPED_TEST(L2MatrixOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor();
}
TYPED_TEST(L2MatrixOperatorTest, constructible_by_factory) {
  this->constructible_by_factory();
}
TYPED_TEST(L2MatrixOperatorTest, is_matrix_operator) {
  this->is_matrix_operator();
}
TYPED_TEST(L2MatrixOperatorTest, correct_for_constant_arguments) {
  this->correct_for_constant_arguments();
}

#if HAVE_DUNE_FEM || HAVE_DUNE_PDELAB
TYPED_TEST(L2MatrixOperatorTest, correct_for_linear_arguments) {
  this->correct_for_linear_arguments();
}
#else
TEST(DISABLED_L2MatrixOperatorTest, correct_for_linear_arguments) {
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif

#if HAVE_DUNE_FEM
TYPED_TEST(L2MatrixOperatorTest, correct_for_quadratic_arguments) {
  this->correct_for_quadratic_arguments();
}
#else
TEST(DISABLED_L2MatrixOperatorTest, correct_for_quadratic_arguments) {
  std::cerr << DSC::colorStringRed("Missing dependencies!") << std::endl;
}
#endif
