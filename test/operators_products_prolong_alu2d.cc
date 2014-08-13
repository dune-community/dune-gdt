// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products_prolong.hh"

#if HAVE_ALUGRID

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID_TWO \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > >

#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID \
  /* all combinations which have Spaces::DiscontinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > >

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


#include <dune/stuff/test/test_main.hh>
