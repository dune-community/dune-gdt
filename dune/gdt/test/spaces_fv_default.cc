// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include<dune/stuff/test/main.hxx>

#include "spaces/fv.hh"
#include "spaces/fv/default.hh"


typedef testing::Types< SPACES_FV
#if HAVE_ALUGRID
                      , SPACES_FV_ALUGRID
#endif
                      > FV_Spaces;

TYPED_TEST_CASE(FV_Space, FV_Spaces);
TYPED_TEST(FV_Space, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(FV_Space, mapper_fulfills_interface) {
  this->mapper_fulfills_interface();
}
TYPED_TEST(FV_Space, basefunctionset_fulfills_interface) {
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(FV_Space, check_for_correct_copy) {
  this->check_for_correct_copy();
}
