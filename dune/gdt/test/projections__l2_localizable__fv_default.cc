// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016 - 2017)

#include <dune/xt/common/test/main.hxx>

#include "projections/l2.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_YASPGRID(1, 1),
                       SPACE_FV_YASPGRID(2, 1),
                       SPACE_FV_YASPGRID(3, 1)
#if HAVE_ALUGRID
                           ,
                       SPACE_FV_ALUCONFORMGRID(2, 1),
                       SPACE_FV_ALUCONFORMGRID(3, 1),
                       SPACE_FV_ALUCUBEGRID(2, 1),
                       SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_ALUGRID
                       >
    SpaceTypes;

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
  this->produces_correct_results(0.096226);
}
