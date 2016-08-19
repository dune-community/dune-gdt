// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/stuff/grid/information.hh>

#include "prolongations.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM_LEVEL(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID_LEVEL(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(ProlongationTest, SpaceTypes);
TYPED_TEST(ProlongationTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::Stuff::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_ProlongationTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
