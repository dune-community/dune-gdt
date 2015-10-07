// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"

#include "products_swipdgpenalty.hh"

typedef testing::Types< SPACE_FV_YASPGRID(1, 1)
//                      , SPACE_FV_YASPGRID(2, 1)
//                      , SPACE_FV_YASPGRID(3, 1)
                      > ConstantSpaces;

TYPED_TEST_CASE(SwipdgPenaltyLocalizableProduct, ConstantSpaces);
TYPED_TEST(SwipdgPenaltyLocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SwipdgPenaltyLocalizableProduct, continuous_arguments) {
  this->continuous_arguments();
}
TYPED_TEST(SwipdgPenaltyLocalizableProduct, discontinuous_arguments) {
  this->discontinuous_arguments();
}
TEST(DISABLED_SwipdgPenaltyLocalizableProduct, dim_2) {}
TEST(DISABLED_SwipdgPenaltyLocalizableProduct, dim_3) {}


#if HAVE_DUNE_FEM

typedef testing::Types<
                        SPACE_DG_FEM_YASPGRID(1, 1, 1)
//                      , SPACE_DG_FEM_YASPGRID(2, 1, 1)
//                      , SPACE_DG_FEM_YASPGRID(3, 1, 1)
                      > LinearSpaces;

TYPED_TEST_CASE(SwipdgPenaltyAssemblableProduct, LinearSpaces);
TYPED_TEST(SwipdgPenaltyAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SwipdgPenaltyAssemblableProduct, continuous_arguments) {
  this->continuous_arguments();
}
TYPED_TEST(SwipdgPenaltyAssemblableProduct, discontinuous_arguments) {
  this->discontinuous_arguments();
}

#else // HAVE_DUNE_FEM

TYPED_TEST_CASE(SwipdgPenaltyAssemblableProduct, ConstantSpaces);
TYPED_TEST(SwipdgPenaltyAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SwipdgPenaltyAssemblableProduct, discontinuous_arguments) {
  this->discontinuous_arguments();
}
TEST(DISABLED_SwipdgPenaltyAssemblableProduct, continuous_arguments) {}

#endif // HAVE_DUNE_FEM

TEST(DISABLED_SwipdgPenaltyAssemblableProduct, dim_2) {}
TEST(DISABLED_SwipdgPenaltyAssemblableProduct, dim_3) {}


TYPED_TEST_CASE(SwipdgPenaltyProduct, ConstantSpaces);
TYPED_TEST(SwipdgPenaltyProduct, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(SwipdgPenaltyProduct, continuous_arguments) {
  this->continuous_arguments();
}
TYPED_TEST(SwipdgPenaltyProduct, discontinuous_arguments) {
  this->discontinuous_arguments();
}
TEST(DISABLED_SwipdgPenaltyProduct, dim_2) {}
TEST(DISABLED_SwipdgPenaltyProduct, dim_3) {}
