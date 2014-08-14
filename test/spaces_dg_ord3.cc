// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "spaces_dg_common.hh"

typedef testing::Types<
#if HAVE_DUNE_FEM
                        Q3_DISCONTINUOUS_LAGRANGE_SPACES_FEM
# if HAVE_ALUGRID
                      , P3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
                      , Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
# endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
                      , Q3_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
# if HAVE_ALUGRID
                      , Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
# endif
#endif // HAVE_DUNE_PDELAB
                      > P3Q3_Spaces;


TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}


#include <dune/stuff/test/test_main.hh>
