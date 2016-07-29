// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/stuff/test/main.hxx>

#include "projections/l2-global.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProjectionOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
