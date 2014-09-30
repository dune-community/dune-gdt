// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "spaces_cg_pdelab.hh"
#include "operators_projections_l2.hh"

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif // HAVE_ALUGRID
                       > SpaceTypes;

TYPED_TEST_CASE(L2ProjectionOperator, SpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(ProjectionOperator, SpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results)
{
  this->produces_correct_results();
}
TYPED_TEST(ProjectionOperator, apply_projection_works)
{
  this->apply_projection_works();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProjectionOperator, produces_correct_results)
{
}
TEST(DISABLED_ProjectionOperator, apply_projection_works)
{
}


#endif // HAVE_DUNE_PDELAB
