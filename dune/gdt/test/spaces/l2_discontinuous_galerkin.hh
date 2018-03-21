// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH
#define DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH

#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceTest
    : public SpaceTestBase<DiscontinuousLagrangeSpace<typename XT::Grid::GridProvider<G>::LeafGridViewType, p, r>, p>
{
  XT::Grid::GridProvider<G> grid_provider;

  DiscontinuousLagrangeSpaceTest(XT::Grid::GridProvider<G>&& grid_provider_)
    : grid_provider(grid_provider_)
  {
    this->grid_view = std::make_shared<typename XT::Grid::GridProvider<G>::LeafGridViewType>(grid_provider.leaf_view());
  }

  void gives_correct_identification()
  {
    ASSERT_NE(this->space, nullptr);
    EXPECT_EQ(SpaceType::discontinuous_lagrange, this->space->type());
    EXPECT_EQ(p, this->space->min_polorder());
    EXPECT_EQ(p, this->space->max_polorder());
    for (int diff_order : {0, 1, 2, 3, 4, 5})
      EXPECT_FALSE(this->space->continuous(diff_order));
    EXPECT_FALSE(this->space->continuous_normal_components());
    EXPECT_TRUE(this->space->is_lagrangian());
  }
}; // struct DiscontinuousLagrangeSpace


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, p>
{
  DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, p>(make_simplicial_grid<G>())
  {
  }
};


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnCubicLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, p>
{
  DiscontinuousLagrangeSpaceOnCubicLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, p>(make_cubic_grid<G>())
  {
  }
};


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnPrismLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, p>
{
  DiscontinuousLagrangeSpaceOnPrismLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, p>(make_prism_grid<G>())
  {
  }
};


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnMixedLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, p>
{
  DiscontinuousLagrangeSpaceOnMixedLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, p>(make_mixed_grid<G>())
  {
  }
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH
