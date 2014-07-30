#include "operators_products_prolong.hh"

typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES> ProlongationOperatorSpaceTypes;

typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES> LagrangeProlongationOperatorSpaceTypes;

#undef L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES


TYPED_TEST_CASE(LagrangeProlongationOperator, LagrangeProlongationOperatorSpaceTypes);
TYPED_TEST(LagrangeProlongationOperator, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(ProlongationOperator, ProlongationOperatorSpaceTypes);
TYPED_TEST(ProlongationOperator, produces_correct_results)
{
  this->produces_correct_results();
}

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
