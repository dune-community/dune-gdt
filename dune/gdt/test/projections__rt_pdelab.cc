// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include "projections.hh"
#include "spaces/rt/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_RT_PDELAB
#if HAVE_ALUGRID
                       ,
                       SPACES_RT_PDELAB_ALUGRID
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(ProjectionTest, SpaceTypes);
TYPED_TEST(ProjectionTest, produces_correct_results)
{
  this->produces_correct_results(9.26e-2);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_ProjectionTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
