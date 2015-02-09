// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "spaces_cg_pdelab.hh"
#include "operators_projections_lagrange.hh"

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif // HAVE_ALUGRID
                       > SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionOperator, SpaceTypes);
TYPED_TEST(LagrangeProjectionOperator, produces_correct_results)
{
  this->produces_correct_results();
}
TYPED_TEST(LagrangeProjectionOperator, free_project_function_works)
{
  this->free_project_function_works();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_LagrangeProjectionOperator, produces_correct_results)
{
}
TEST(DISABLED_LagrangeProjectionOperator, free_project_function_works)
{
}


#endif // HAVE_DUNE_PDELAB
