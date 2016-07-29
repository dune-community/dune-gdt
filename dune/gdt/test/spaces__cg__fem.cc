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
#include "spaces/cg/fem.hh"

#if HAVE_DUNE_FEM

typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       > CG_Spaces_Fem;

TYPED_TEST_CASE(LeafCG_Space, CG_Spaces_Fem);
TYPED_TEST(LeafCG_Space, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(LeafCG_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}
TYPED_TEST(LeafCG_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(LeafCG_Space, check_for_correct_copy)
{
  this->check_for_correct_copy();
}

typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       > P1Q1_CG_Spaces_Fem;

TYPED_TEST_CASE(LeafP1Q1_CG_Space, P1Q1_CG_Spaces_Fem);
TYPED_TEST(LeafP1Q1_CG_Space, fulfills_continuous_interface)
{
  this->fulfills_continuous_interface();
}
TYPED_TEST(LeafP1Q1_CG_Space, maps_correctly)
{
  this->maps_correctly();
}

#else // HAVE_DUNE_FEM

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

#endif // HAVE_DUNE_FEM
