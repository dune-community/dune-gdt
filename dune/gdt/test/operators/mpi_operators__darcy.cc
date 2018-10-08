// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include <dune/gdt/test/spaces/cg.hh>
#include <dune/gdt/test/spaces/rt.hh>

#include "darcy.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_ALUGRID

typedef testing::Types<
    /*std::pair< SPACE_CG_ALUCONFORMGRID(2, 1, 1), SPACE_CG_ALUCONFORMGRID(2, 2, 1) > // <- TODO: enable once #40 is resolved
                      ,*/ std::pair<SPACE_CG_ALUCONFORMGRID(2, 1, 1), SPACE_RT_ALUCONFORMGRID(2)>>
    SpaceTypes;

TYPED_TEST_CASE(DarcyOperatorTest, SpaceTypes);
TYPED_TEST(DarcyOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_ALUGRID


TEST(DISABLED_DarcyOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_ALUGRID
