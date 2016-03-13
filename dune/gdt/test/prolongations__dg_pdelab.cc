// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "prolongations.hh"
#include "spaces/dg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types< SPACES_DG_PDELAB_LEVEL(1)
#if HAVE_ALUGRID
                      , SPACES_DG_PDELAB_ALUGRID_LEVEL(1)
#endif
                      > SpaceTypes;

TYPED_TEST_CASE(ProlongationTest, SpaceTypes);
TYPED_TEST(ProlongationTest, produces_correct_results) {
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_ProlongationTest, produces_correct_results) {}


#endif // HAVE_DUNE_PDELAB
