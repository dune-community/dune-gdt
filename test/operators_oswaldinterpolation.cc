// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#include "spaces_dg_fem.hh"
#include "spaces_dg_pdelab.hh"
#include "operators_oswaldinterpolation.hh"

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_DG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_DG_FEM_ALUGRID(1)
#endif
                       > SpaceTypes;

TYPED_TEST_CASE(Oswald_Interpolation_Operator, SpaceTypes);
TYPED_TEST(Oswald_Interpolation_Operator, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_Oswald_Interpolation_Operator, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
