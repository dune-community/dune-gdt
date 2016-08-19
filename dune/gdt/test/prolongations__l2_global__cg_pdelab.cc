// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include "prolongations/l2-global.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB_LEVEL(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID_LEVEL(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProlongationOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProlongationOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProlongationOperatorTest, produces_correct_results)
{
  this->produces_correct_results(pdelab_cg_tolerance(*this));
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2GlobalProlongationOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProlongationOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProlongationOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
