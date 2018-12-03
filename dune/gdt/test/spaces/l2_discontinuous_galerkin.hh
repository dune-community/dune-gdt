// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH
#define DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH

#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class G, size_t r, class R, int p>
struct DiscontinuousLagrangeSpaceTest
  : public SpaceTestBase<DiscontinuousLagrangeSpace<typename XT::Grid::GridProvider<G>::LeafGridViewType, r, R>, p>
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


template <class G, size_t r, class R, int p>
struct DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, R, p>
{
  DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, R, p>(make_simplicial_grid<G>())
  {}
};


template <class G, size_t r, class R, int p>
struct DiscontinuousLagrangeSpaceOnCubicLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, R, p>
{
  DiscontinuousLagrangeSpaceOnCubicLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, R, p>(make_cubic_grid<G>())
  {}
};


template <class G, size_t r, class R, int p>
struct DiscontinuousLagrangeSpaceOnPrismLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, R, p>
{
  DiscontinuousLagrangeSpaceOnPrismLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, R, p>(make_prism_grid<G>())
  {}
};


template <class G, size_t r, class R, int p>
struct DiscontinuousLagrangeSpaceOnMixedLeafViewTest : public DiscontinuousLagrangeSpaceTest<G, r, R, p>
{
  DiscontinuousLagrangeSpaceOnMixedLeafViewTest()
    : DiscontinuousLagrangeSpaceTest<G, r, R, p>(make_mixed_grid<G>())
  {}
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH
