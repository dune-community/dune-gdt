// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "prolongations/l2-local.hh"
#include "spaces/dg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types< SPACES_DG_PDELAB_LEVEL(1)
#if HAVE_ALUGRID && !defined(__GNUC__)
                      , SPACES_DG_PDELAB_ALUGRID_LEVEL(1)
#endif
                      > SpaceTypes;

TYPED_TEST_CASE(L2LocalProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProlongationOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor();
}
TYPED_TEST(L2LocalProlongationOperatorTest, constructible_by_factory) {
  this->constructible_by_factory();
}
TYPED_TEST(L2LocalProlongationOperatorTest, produces_correct_results) {
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2LocalProlongationOperatorTest, constructible_by_ctor) {}
TEST(DISABLED_L2LocalProlongationOperatorTest, constructible_by_factory) {}
TEST(DISABLED_L2LocalProlongationOperatorTest, produces_correct_results) {}


#endif // HAVE_DUNE_PDELAB
