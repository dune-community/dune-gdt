// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)

#include <dune/stuff/test/main.hxx>

#include "projections.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1), SPACE_FV_YASPGRID(1, 1),
                       SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)
#if HAVE_ALUGRID
                                                    ,
                       SPACE_FV_ALUCONFORMGRID(2, 1), SPACE_FV_ALUCONFORMGRID(3, 1), SPACE_FV_ALUCUBEGRID(2, 1),
                       SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_ALUGRID
                       >
    SpaceTypes;

TYPED_TEST_CASE(ProjectionTest, SpaceTypes);
TYPED_TEST(ProjectionTest, produces_correct_results)
{
  this->produces_correct_results(9.63e-2);
}
