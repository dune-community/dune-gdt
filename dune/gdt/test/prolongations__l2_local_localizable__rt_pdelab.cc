// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "prolongations/l2.hh"
#include "spaces/rt/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types< SPACES_RT_PDELAB_LEVEL
#if HAVE_ALUGRID && !defined(__GNUC__)
                      , SPACES_RT_PDELAB_ALUGRID_LEVEL
#endif
                      > SpaceTypes;

TYPED_TEST_CASE(L2LocalProlongationLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProlongationLocalizableOperatorTest, produces_correct_results) {
  this->produces_correct_results(0.204125);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2LocalProlongationLocalizableOperatorTest, produces_correct_results) {}


#endif // HAVE_DUNE_PDELAB
