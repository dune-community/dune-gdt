// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#if 0
#include "operators_products_prolong.hh"

#if HAVE_ALUGRID

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO \
  /* all combinations which have Spaces::CG::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  \
  /* all combinations which have Spaces::CG::FemBased as FineSpaceType */ \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > >

#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID \
  /* all combinations which have Spaces::DG::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::CG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DG::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > >

#endif // HAVE_ALUGRID

typedef testing::Types<
#if HAVE_ALUGRID
                        LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO
#endif
                      > ProlongationOperatorSpaceTypes;


typedef testing::Types<
#if HAVE_ALUGRID
                        LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO
#endif
                      > LagrangeProlongationOperatorSpaceTypes;


#undef L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES

TYPED_TEST_CASE(LagrangeProlongationOperator, LagrangeProlongationOperatorSpaceTypes);
TYPED_TEST(LagrangeProlongationOperator, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(ProlongationOperator, ProlongationOperatorSpaceTypes);
TYPED_TEST(ProlongationOperator, produces_correct_results) {
  this->produces_correct_results();
}
#endif // 0

TEST(DISABLED_LagrangeProlongationOperator, produces_correct_results) {}
TEST(DISABLED_ProlongationOperator, produces_correct_results) {}
