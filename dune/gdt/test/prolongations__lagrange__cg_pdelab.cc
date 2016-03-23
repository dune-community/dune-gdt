// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- This one has to come first!

#include "prolongations/lagrange.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB_LEVEL(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID_LEVEL(1)
#endif
                       > SpaceTypes;

TYPED_TEST_CASE(LagrangeProlongationOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProlongationOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProlongationOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProlongationOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(LagrangeProlongationOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_LagrangeProlongationOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProlongationOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProlongationOperatorTest, free_function_callable)
{
}
TEST(DISABLED_LagrangeProlongationOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
