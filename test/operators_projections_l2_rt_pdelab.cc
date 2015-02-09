// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/stuff/test/main.hxx>

#include "spaces_rt_pdelab.hh"
#include "operators_projections_l2.hh"

#if HAVE_DUNE_PDELAB


typedef testing::Types< SPACES_RT_PDELAB
#if HAVE_ALUGRID
                      , SPACES_RT_PDELAB_ALUGRID
#endif // HAVE_ALUGRID
                      > SpaceTypes;

TYPED_TEST_CASE(L2ProjectionOperator, SpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
 this->produces_correct_results(0.0925927);
}
TYPED_TEST(L2ProjectionOperator, free_project_l2_function_works) {
 this->free_project_l2_function_works(0.0925927);
}

TYPED_TEST_CASE(ProjectionOperator, SpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results) {
 this->produces_correct_results(0.0925927);
}
TYPED_TEST(ProjectionOperator, free_project_function_works) {
 this->free_project_function_works(0.0925927);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProjectionOperator, produces_correct_results) {}
TEST(DISABLED_L2ProjectionOperator, free_project_l2_function_works) {}
TEST(DISABLED_ProjectionOperator, produces_correct_results) {}
TEST(DISABLED_ProjectionOperator, free_project_function_works) {}


#endif // HAVE_DUNE_PDELAB
