// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/stuff/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(LagrangeProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, free_function_callable)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
