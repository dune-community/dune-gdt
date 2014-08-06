#include "spaces_cg.hh"

#if HAVE_DUNE_FEM

typedef testing::Types<P1_CONTINUOUS_LAGRANGE_SPACES_FEM, Q1_CONTINUOUS_LAGRANGE_SPACES_FEM>
    P1Q1_Continuous_Lagrange_Spaces;


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

#endif // HAVE_DUNE_FEM

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
