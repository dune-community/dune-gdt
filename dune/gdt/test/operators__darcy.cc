// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include "spaces/cg/fem.hh"
#include "spaces/rt/pdelab.hh"

#include "operators/darcy.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM && HAVE_DUNE_PDELAB && HAVE_ALUGRID && !defined(__GNUC__)

typedef testing::Types<
                        /*std::pair< SPACE_CG_FEM_ALUCONFORMGRID(2, 1, 1), SPACE_CG_FEM_ALUCONFORMGRID(2, 2, 1) > // <- TODO: enable once #40 is resolved
                      ,*/ std::pair< SPACE_CG_FEM_ALUCONFORMGRID(2, 1, 1), SPACE_RT_PDELAB_ALUCONFORMGRID(2) >
                      > SpaceTypes;

TYPED_TEST_CASE(DarcyOperatorTest, SpaceTypes);
TYPED_TEST(DarcyOperatorTest, produces_correct_results) {
  this->produces_correct_results();
}


#else // HAVE_DUNE_FEM && HAVE_DUNE_PDELAB && HAVE_ALUGRID && !defined(__GNUC__)


TEST(DISABLED_DarcyOperatorTest, produces_correct_results) {}


#endif // HAVE_DUNE_FEM && HAVE_DUNE_PDELAB && HAVE_ALUGRID && !defined(__GNUC__)
