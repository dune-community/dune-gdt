// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#include <dune/xt/common/test/main.hxx>

#include "prolongations/l2-local.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_YASPGRID_LEVEL(1, 1),
                       SPACE_FV_YASPGRID_LEVEL(2, 1),
                       SPACE_FV_YASPGRID_LEVEL(3, 1)
#if HAVE_ALUGRID
                           ,
                       SPACE_FV_ALUCONFORMGRID_LEVEL(2, 1),
                       SPACE_FV_ALUCONFORMGRID_LEVEL(3, 1),
                       SPACE_FV_ALUCUBEGRID_LEVEL(2, 1),
                       SPACE_FV_ALUCUBEGRID_LEVEL(3, 1)
#endif // HAVE_ALUGRID
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2LocalProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProlongationOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor(1.45e-1);
}
TYPED_TEST(L2LocalProlongationOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory(1.45e-1);
}
TYPED_TEST(L2LocalProlongationOperatorTest, produces_correct_results)
{
  this->produces_correct_results(1.45e-1);
}
