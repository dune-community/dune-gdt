// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include <dune/stuff/test/main.hxx>

#include "projections/l2.hh"
#include "spaces/dg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_DG_PDELAB(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_DG_PDELAB_ALUGRID(1)
#endif
                       > SpaceTypes;

TYPED_TEST_CASE(L2ProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2ProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2ProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2ProjectionLocalizableOperatorTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::Stuff::Grid::is_alugrid<Grid>::value ? LocalizableProjectionOperator_alugrid_tolerance
                                                                    : LocalizableProjectionOperator_default_tolerance;
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2ProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2ProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
