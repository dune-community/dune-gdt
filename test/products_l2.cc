// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"
#include "products_l2.hh"


typedef testing::Types< SPACE_FV_YASPGRID(1, 1)
                      , SPACE_FV_YASPGRID(2, 1)
                      , SPACE_FV_YASPGRID(3, 1)
                      > ConstantSpaces;

TYPED_TEST_CASE(L2LocalizableProduct, ConstantSpaces);
TYPED_TEST(L2LocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2LocalizableProduct, ConstantSpaces);
TYPED_TEST(L2LocalizableProduct, constant_arguments) {
  this->constant_arguments();
}

TYPED_TEST_CASE(L2LocalizableProduct, ConstantSpaces);
TYPED_TEST(L2LocalizableProduct, linear_arguments) {
  this->linear_arguments();
}

TYPED_TEST_CASE(L2LocalizableProduct, ConstantSpaces);
TYPED_TEST(L2LocalizableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}


#if HAVE_DUNE_FEM

typedef testing::Types<
                        SPACE_DG_FEM_YASPGRID(1, 1, 2)
                      , SPACE_DG_FEM_YASPGRID(2, 1, 2)
                      , SPACE_DG_FEM_YASPGRID(3, 1, 2)
                      > QuadraticSpaces;

TYPED_TEST_CASE(L2AssemblableProduct, QuadraticSpaces);
TYPED_TEST(L2AssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(L2AssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(L2AssemblableProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(L2AssemblableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}

#else // HAVE_DUNE_FEM

TYPED_TEST_CASE(L2AssemblableProduct, ConstantSpaces);
TYPED_TEST(L2AssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(L2AssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TEST(DISABLED_L2AssemblableProduct, linear_arguments)    {}
TEST(DISABLED_L2AssemblableProduct, quadratic_arguments) {}

#endif // HAVE_DUNE_FEM


TYPED_TEST_CASE(L2Product, ConstantSpaces);
TYPED_TEST(L2Product, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(L2Product, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(L2Product, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(L2Product, quadratic_arguments) {
  this->quadratic_arguments();
}

