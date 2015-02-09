// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "spaces_dg_fem.hh"
#include "operators_projections_l2.hh"

#if HAVE_DUNE_FEM


typedef testing::Types< SPACES_DG_FEM(1)
                      , SPACES_DG_FEM(2)
#if HAVE_ALUGRID
                      , SPACES_DG_FEM_ALUGRID(1)
                      , SPACES_DG_FEM_ALUGRID(2)
#endif // HAVE_ALUGRID
                      > SpaceTypes;

TYPED_TEST_CASE(L2ProjectionOperator, SpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}
TYPED_TEST(L2ProjectionOperator, free_project_l2_function_works) {
 this->free_project_l2_function_works();
}

TYPED_TEST_CASE(ProjectionOperator, SpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}
TYPED_TEST(ProjectionOperator, free_project_function_works) {
 this->free_project_function_works();
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_L2ProjectionOperator, produces_correct_results) {}
TEST(DISABLED_L2ProjectionOperator, free_project_l2_function_works) {}
TEST(DISABLED_ProjectionOperator, produces_correct_results) {}
TEST(DISABLED_ProjectionOperator, free_project_function_works) {}


#endif // HAVE_DUNE_FEM
