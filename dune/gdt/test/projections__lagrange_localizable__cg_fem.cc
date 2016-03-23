// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       > SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
