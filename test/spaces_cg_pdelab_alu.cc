// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "spaces_cg.hh"


typedef testing::Types<
#if HAVE_DUNE_PDELAB
#if HAVE_ALUGRID
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB, Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB

    > P1Q1_Continuous_Lagrange_Spaces;


TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, fulfills_continuous_interface)
{
  this->fulfills_continuous_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, maps_correctly)
{
  this->maps_correctly();
}


#include <dune/stuff/test/test_main.hh>
