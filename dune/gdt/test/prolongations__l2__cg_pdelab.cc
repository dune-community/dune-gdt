// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#include <dune/xt/common/test/main.hxx>

#include "prolongations/l2.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB_LEVEL(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID_LEVEL(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2ProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2ProlongationOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(L2ProlongationOperatorTest, produces_correct_results)
{
  this->produces_correct_results(pdelab_cg_tolerance(*this));
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProlongationOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2ProlongationOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2ProlongationOperatorTest, free_function_callable)
{
}
TEST(DISABLED_L2ProlongationOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
