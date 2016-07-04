// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)

#include <dune/stuff/test/main.hxx>

#include "projections/l2.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1), SPACE_FV_YASPGRID(1, 1),
                       SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)
#if HAVE_DUNE_ALUGRID
                                                    ,
                       SPACE_FV_ALUCONFORMGRID(2, 1), SPACE_FV_ALUCONFORMGRID(3, 1), SPACE_FV_ALUCUBEGRID(2, 1),
                       SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_DUNE_ALUGRID
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2ProjectionOperatorTest, SpaceTypes);
TYPED_TEST(L2ProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2ProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2ProjectionOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(L2ProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results(0.096226);
}
