// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)

#define DUNE_GDT_TEST_SPACES_RT_CHECK 1

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include "spaces/rt.hh"
#include "spaces/rt/pdelab.hh"

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_RT_PDELAB
#if HAVE_ALUGRID
                       ,
                       SPACES_RT_PDELAB_ALUGRID
#endif
                       >
    RT_Spaces;

TYPED_TEST_CASE(RT_Space, RT_Spaces);
TYPED_TEST(RT_Space, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(RT_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}
TYPED_TEST(RT_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(RT_Space, check_for_correct_copy)
{
  this->check_for_correct_copy();
}
TYPED_TEST(RT_Space, matches_raviart_thomas_signature)
{
  this->matches_raviart_thomas_signature();
}


#if HAVE_ALUGRID

typedef testing::Types<SPACE_RT_PDELAB_ALUCONFORMGRID(2)> RT_2d_Spaces;

TYPED_TEST_CASE(RT_2d_simplicial_Space, RT_2d_Spaces);
TYPED_TEST(RT_2d_simplicial_Space, fulfills_raviart_thomas_2d_simplicial_interface)
{
  this->fulfills_raviart_thomas_2d_simplicial_interface();
}

#else // HAVE_ALUGRID

TEST(DISABLED_RT_2d_simplicial_Space, fulfills_raviart_thomas_2d_simplicial_interface)
{
}

#endif // HAVE_ALUGRID
#else // HAVE_DUNE_PDELAB


TEST(DISABLED_RT_Space, fulfills_interface)
{
}
TEST(DISABLED_RT_Space, mapper_fulfills_interface)
{
}
TEST(DISABLED_RT_Space, basefunctionset_fulfills_interface)
{
}
TEST(DISABLED_RT_Space, check_for_correct_copy)
{
}
TEST(DISABLED_RT_Space, matches_raviart_thomas_signature)
{
}
TEST(DISABLED_RT_2d_simplicial_Space, fulfills_raviart_thomas_2d_simplicial_interface)
{
}


#endif // HAVE_DUNE_PDELAB
