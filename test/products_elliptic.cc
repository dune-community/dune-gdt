// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"
#include "spaces_cg_pdelab.hh"

#include "products_elliptic.hh"


typedef testing::Types< SPACE_FV_YASPGRID(1, 1)
                      , SPACE_FV_YASPGRID(2, 1)
                      , SPACE_FV_YASPGRID(3, 1)
                      > ConstantSpaces;

TYPED_TEST_CASE(EllipticLocalizableProduct, ConstantSpaces);
TYPED_TEST(EllipticLocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(EllipticLocalizableProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(EllipticLocalizableProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(EllipticLocalizableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}

TYPED_TEST_CASE(SimplifiedEllipticLocalizableProduct, ConstantSpaces);
TYPED_TEST(SimplifiedEllipticLocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SimplifiedEllipticLocalizableProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(SimplifiedEllipticLocalizableProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(SimplifiedEllipticLocalizableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}


#if HAVE_DUNE_FEM

typedef testing::Types<
                        SPACE_DG_FEM_SGRID(1, 1, 3)
                      , SPACE_DG_FEM_SGRID(2, 1, 3)
                      , SPACE_DG_FEM_SGRID(3, 1, 3)
                      > CubicSpaces;

TYPED_TEST_CASE(EllipticAssemblableProduct, CubicSpaces);
TYPED_TEST(EllipticAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(EllipticAssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(EllipticAssemblableProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(EllipticAssemblableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}

TYPED_TEST_CASE(SimplifiedEllipticAssemblableProduct, CubicSpaces);
TYPED_TEST(SimplifiedEllipticAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SimplifiedEllipticAssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(SimplifiedEllipticAssemblableProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(SimplifiedEllipticAssemblableProduct, quadratic_arguments) {
  this->quadratic_arguments();
}

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types<
                        SPACE_CG_PDELAB_YASPGRID(1, 1, 1)
                      , SPACE_CG_PDELAB_YASPGRID(2, 1, 1)
                      , SPACE_CG_PDELAB_YASPGRID(3, 1, 1)
                      > LinearSpaces;

TYPED_TEST_CASE(EllipticAssemblableProduct, LinearSpaces);
TYPED_TEST(EllipticAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(EllipticAssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TEST(DISABLED_EllipticAssemblableProduct, linear_arguments)    {}
TEST(DISABLED_EllipticAssemblableProduct, quadratic_arguments) {}

TYPED_TEST_CASE(SimplifiedEllipticAssemblableProduct, LinearSpaces);
TYPED_TEST(SimplifiedEllipticAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SimplifiedEllipticAssemblableProduct, constant_arguments) {
  this->constant_arguments();
}
TEST(DISABLED_SimplifiedEllipticAssemblableProduct, linear_arguments)    {}
TEST(DISABLED_SimplifiedEllipticAssemblableProduct, quadratic_arguments) {}

#else // HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

TEST(DISABLED_EllipticAssemblableProduct, fulfills_interface)  {}
TEST(DISABLED_EllipticAssemblableProduct, constant_arguments)  {}
TEST(DISABLED_EllipticAssemblableProduct, linear_arguments)    {}
TEST(DISABLED_EllipticAssemblableProduct, quadratic_arguments) {}

TEST(DISABLED_SimplifiedEllipticAssemblableProduct, fulfills_interface)  {}
TEST(DISABLED_SimplifiedEllipticAssemblableProduct, constant_arguments)  {}
TEST(DISABLED_SimplifiedEllipticAssemblableProduct, linear_arguments)    {}
TEST(DISABLED_SimplifiedEllipticAssemblableProduct, quadratic_arguments) {}

#endif // HAVE_DUNE_PDELAB // HAVE_DUNE_FEM


TYPED_TEST_CASE(EllipticProduct, ConstantSpaces);
TYPED_TEST(EllipticProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(EllipticProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(EllipticProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(EllipticProduct, quadratic_arguments) {
  this->quadratic_arguments();
}

TYPED_TEST_CASE(SimplifiedEllipticProduct, ConstantSpaces);
TYPED_TEST(SimplifiedEllipticProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SimplifiedEllipticProduct, constant_arguments) {
  this->constant_arguments();
}
TYPED_TEST(SimplifiedEllipticProduct, linear_arguments) {
  this->linear_arguments();
}
TYPED_TEST(SimplifiedEllipticProduct, quadratic_arguments) {
  this->quadratic_arguments();
}
