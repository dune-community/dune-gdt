// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "spaces_dg_common.hh"

typedef testing::Types<
#if HAVE_DUNE_FEM
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > P1Q1_Spaces;

TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, fulfills_interface)
{
  this->fulfills_interface();
}


TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}


TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Discontinuous_Lagrange, P1Q1_Spaces);
TYPED_TEST(P1Q1_Discontinuous_Lagrange, maps_correctly)
{
  this->maps_correctly();
}

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
