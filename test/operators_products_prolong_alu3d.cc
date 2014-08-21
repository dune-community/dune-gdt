// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products_prolong.hh"

#if HAVE_ALUGRID

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_THREE                                                       \
  std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,              \
            Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,             \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluCube3dLevelGridPartType, 1, double, 1>,             \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluCube3dLevelGridPartType, 1, double, 1>>,            \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,         \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>


/* those below do not work in 3d any more! */
#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID                                                                   \
  std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,              \
            Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,          \
      std::pair<Dune::GDT::Spaces::ContinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,          \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>,      \
      std::pair<Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>,       \
                Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLevelGridPartType, 1, double, 1>>

#endif // HAVE_ALUGRID


typedef testing::Types<
#if HAVE_ALUGRID
    L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
    > L2ProlongationOperatorSpaceTypes;


typedef testing::Types<
#if HAVE_ALUGRID
    LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_THREE, L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif // HAVE_ALUGRID
    > ProlongationOperatorSpaceTypes;


typedef testing::Types<
#if HAVE_ALUGRID
    LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_THREE
#endif
    > LagrangeProlongationOperatorSpaceTypes;

#undef L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES


TYPED_TEST_CASE(L2ProlongationOperator, L2ProlongationOperatorSpaceTypes);
TYPED_TEST(L2ProlongationOperator, produces_correct_results)
{
  this->produces_correct_results();
}


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


#include <dune/stuff/test/test_main.cxx>
