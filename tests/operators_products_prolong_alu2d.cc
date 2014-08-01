#include "operators_products_prolong.hh"

#if HAVE_ALUGRID

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO                                                         \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */                              \
  std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,              \
            Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,             \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>          \
                                                                                                                       \
      , std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,        \
                  Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,       \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>          \
                                                                                                                       \
      /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */                          \
      , std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,        \
                  Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,       \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>          \
                                                                                                                       \
      , std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,        \
                  Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,       \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>

#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID                                                                   \
  /* all combinations which have Spaces::DiscontinuousLagrange::FemBased as FineSpaceType */                           \
  std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,              \
            Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,          \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>,      \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLevelGridPartType, 1, double, 1>>       \
                                                                                                                       \
      , std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,        \
                  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,    \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>,      \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLevelGridPartType, 1, double, 1>>

typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO> ProlongationOperatorSpaceTypes;


typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO> LagrangeProlongationOperatorSpaceTypes;


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

#endif

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
