// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "spaces_dg.hh"
#include "spaces_dg_fem.hh"

#if HAVE_DUNE_FEM


typedef testing::Types<
                        SPACES_DG_FEM(1)
                      , SPACES_DG_FEM(2)
# if HAVE_ALUGRID
                      , SPACES_DG_FEM_ALUGRID(1)
                      , SPACES_DG_FEM_ALUGRID(2)
# endif
                      > DG_Spaces_Fem;

TYPED_TEST_CASE(DG_Space, DG_Spaces_Fem);
TYPED_TEST(DG_Space, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(DG_Space, copy_and_move_ctor) {
  this->copy_and_move_ctor();
}
TYPED_TEST(DG_Space, mapper_fulfills_interface) {
  this->mapper_fulfills_interface();
}
TYPED_TEST(DG_Space, basefunctionset_fulfills_interface) {
  this->basefunctionset_fulfills_interface();
}


typedef testing::Types<
                        SPACES_DG_FEM(1)
# if HAVE_ALUGRID
                      , SPACES_DG_FEM_ALUGRID(1)
# endif
                      > P1Q1_DG_Spaces_Fem;

TYPED_TEST_CASE(P1Q1_DG_Space, P1Q1_DG_Spaces_Fem);
TYPED_TEST(P1Q1_DG_Space, maps_correctly) {
  this->maps_correctly();
}

#else // HAVE_DUNE_FEM

TEST(DISABLED_SpaceBase, fulfills_interface)                 {}
TEST(DISABLED_SpaceBase, copy_and_move_ctor)                 {}
TEST(DISABLED_SpaceBase, mapper_fulfills_interface)          {}
TEST(DISABLED_SpaceBase, basefunctionset_fulfills_interface) {}
TEST(DISABLED_P1Q1_DG_Space, maps_correctly)                 {}

#endif // HAVE_DUNE_FEM
