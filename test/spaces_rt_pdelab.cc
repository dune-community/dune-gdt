// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_GDT_TEST_SPACES_RT_CHECK 1

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include "spaces_rt.hh"
#include "spaces_rt_pdelab.hh"


typedef testing::Types< SPACES_RT_PDELAB
#if HAVE_ALUGRID
                      , SPACES_RT_PDELAB_ALUGRID
#endif
                      > RT_Spaces;

TYPED_TEST_CASE(RT_Space, RT_Spaces);
TYPED_TEST(RT_Space, fulfills_interface) {
  this->fulfills_interface();
}
TYPED_TEST(RT_Space, mapper_fulfills_interface) {
  this->mapper_fulfills_interface();
}
TYPED_TEST(RT_Space, basefunctionset_fulfills_interface) {
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(RT_Space, check_for_correct_copy) {
  this->check_for_correct_copy();
}
