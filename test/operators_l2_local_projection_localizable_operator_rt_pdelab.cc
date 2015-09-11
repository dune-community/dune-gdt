// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING 1

#include <dune/stuff/test/main.hxx>

#include <dune/gdt/tests/operators/projections/l2.hh>

#include "spaces_rt_pdelab.hh"

using namespace Dune::GDT::Tests;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_RT_PDELAB
#if HAVE_ALUGRID
                       ,
                       SPACES_RT_PDELAB_ALUGRID
#endif
                       > SpaceTypes;

TYPED_TEST_CASE(L2LocalProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2LocalProjectionLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results(0.0925927);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2LocalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2LocalProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
