// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1

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

TYPED_TEST_CASE(L2ProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor(this->dimDomain == 3 ? 2.05e-1 : 1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_factory) {
  this->constructible_by_factory(this->dimDomain == 3 ? 2.05e-1 : 1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, free_function_callable) {
  this->free_function_callable(this->dimDomain == 3 ? 2.05e-1 : 1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, produces_correct_results) {
  this->produces_correct_results(this->dimDomain == 3 ? 2.05e-1 : 1.45e-1);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProlongationOperatorTest, constructible_by_ctor) {}
TEST(DISABLED_L2ProlongationOperatorTest, constructible_by_factory) {}
TEST(DISABLED_L2ProlongationOperatorTest, free_function_callable) {}
TEST(DISABLED_L2ProlongationOperatorTest, produces_correct_results) {}


#endif // HAVE_DUNE_PDELAB
