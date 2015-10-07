// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "operators_projections_l2.hh"


typedef testing::Types< SPACE_FV_YASPGRID(1, 1)
                      , SPACE_FV_YASPGRID(2, 1)
                      , SPACE_FV_YASPGRID(3, 1)
#if HAVE_ALUGRID
                      , SPACE_FV_ALUCONFORMGRID(2, 1)
                      , SPACE_FV_ALUCONFORMGRID(3, 1)
                      , SPACE_FV_ALUCUBEGRID(2, 1)
                      , SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_ALUGRID
                      > SpaceTypes;

TYPED_TEST_CASE(L2ProjectionOperator, SpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
 this->produces_correct_results(0.096226);
}
TYPED_TEST(L2ProjectionOperator, free_project_l2_function_works) {
 this->free_project_l2_function_works(0.096226);
}

TYPED_TEST_CASE(ProjectionOperator, SpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results) {
 this->produces_correct_results(0.096226);
}
TYPED_TEST(ProjectionOperator, free_project_function_works) {
 this->free_project_function_works(0.096226);
}
