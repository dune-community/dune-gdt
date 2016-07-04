// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#include <dune/stuff/test/main.hxx>

#include "spaces/cg.hh"
#include "spaces/cg/pdelab.hh"

#if HAVE_DUNE_PDELAB

typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    CG_Spaces_Pdelab;

TYPED_TEST_CASE(CG_Space, CG_Spaces_Pdelab);
TYPED_TEST(CG_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST(CG_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}
TYPED_TEST(CG_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(CG_Space, check_for_correct_copy)
{
  this->check_for_correct_copy();
}


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    P1Q1_CG_Spaces_Pdelab;

TYPED_TEST_CASE(P1Q1_CG_Space, P1Q1_CG_Spaces_Pdelab);
TYPED_TEST(P1Q1_CG_Space, fulfills_continuous_interface)
{
  this->fulfills_continuous_interface();
}
TYPED_TEST(P1Q1_CG_Space, maps_correctly)
{
  this->maps_correctly();
}

#else // HAVE_DUNE_PDELAB

TEST(DISABLED_CG_Space, fulfills_interface)
{
}
TEST(DISABLED_CG_Space, mapper_fulfills_interface)
{
}
TEST(DISABLED_CG_Space, basefunctionset_fulfills_interface)
{
}
TEST(DISABLED_CG_Space, check_for_correct_copy)
{
}
TEST(DISABLED_P1Q1_CG_Space, fulfills_continuous_interface)
{
}
TEST(DISABLED_P1Q1_CG_Space, maps_correctly)
{
}

#endif // HAVE_DUNE_PDELAB
