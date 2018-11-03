// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_SPACES_H1_CONTINUOUS_LAGRANGE_HH
#define DUNE_GDT_TEST_SPACES_H1_CONTINUOUS_LAGRANGE_HH

#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class G, class R, int p>
struct ContinuousLagrangeSpaceTest
    : public SpaceTestBase<ContinuousLagrangeSpace<typename XT::Grid::GridProvider<G>::LeafGridViewType, R>, p>
{
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;

  XT::Grid::GridProvider<G> grid_provider;

  ContinuousLagrangeSpaceTest(XT::Grid::GridProvider<G>&& grid_provider_)
    : grid_provider(grid_provider_)
  {
    this->grid_view = std::make_shared<typename XT::Grid::GridProvider<G>::LeafGridViewType>(grid_provider.leaf_view());
  }

  void gives_correct_identification()
  {
    ASSERT_NE(this->space, nullptr);
    EXPECT_EQ(SpaceType::continuous_lagrange, this->space->type());
    EXPECT_EQ(p, this->space->min_polorder());
    EXPECT_EQ(p, this->space->max_polorder());
    EXPECT_TRUE(this->space->continuous(0));
    for (int diff_order : {1, 2, 3, 4, 5})
      EXPECT_FALSE(this->space->continuous(diff_order));
    EXPECT_TRUE(this->space->continuous_normal_components());
    EXPECT_TRUE(this->space->is_lagrangian());
  }

  void mapper_maps_correctly()
  {
    ASSERT_NE(this->space, nullptr);
    // collect all global ids that are associated with a global lagrange point
    std::map<FieldVector<D, d>, std::set<size_t>, XT::Common::FieldVectorLess>
        global_lagrange_point_to_global_indices_map;
    for (auto&& element : elements(*this->grid_view)) {
      const auto global_indices = this->space->mapper().global_indices(element);
      EXPECT_LE(this->space->mapper().local_size(element), global_indices.size());
      const auto lagrange_points = this->space->finite_element(element.geometry().type()).lagrange_points();
      EXPECT_EQ(lagrange_points.size(), this->space->mapper().local_size(element));
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto global_lagrange_point = element.geometry().global(lagrange_points[ii]);
        const auto global_index = this->space->mapper().global_index(element, ii);
        EXPECT_EQ(global_indices[ii], global_index);
        global_lagrange_point_to_global_indices_map[global_lagrange_point].insert(global_index);
      }
    }
    // check that all global lagrange points have indeed one and only one global DoF id ...
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : global_lagrange_point_to_global_indices_map) {
      const auto global_DoF_indices_per_point = entry.second;
      EXPECT_EQ(global_DoF_indices_per_point.size(), 1);
      global_DoF_indices.insert(*(global_DoF_indices_per_point.begin()));
    }
    EXPECT_EQ(global_lagrange_point_to_global_indices_map.size(), global_DoF_indices.size());
    // ... and that the numbering is consecutive
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      EXPECT_EQ(global_DoF_id, count);
      ++count;
    }
    EXPECT_EQ(global_DoF_indices.size(), this->space->mapper().size());
  } // ... mapper_maps_correctly(...)
}; // struct ContinuousLagrangeSpace


template <class G, class R, int p>
struct ContinuousLagrangeSpaceOnSimplicialLeafViewTest : public ContinuousLagrangeSpaceTest<G, R, p>
{
  ContinuousLagrangeSpaceOnSimplicialLeafViewTest()
    : ContinuousLagrangeSpaceTest<G, R, p>(make_simplicial_grid<G>())
  {
  }
};


template <class G, class R, int p>
struct ContinuousLagrangeSpaceOnCubicLeafViewTest : public ContinuousLagrangeSpaceTest<G, R, p>
{
  ContinuousLagrangeSpaceOnCubicLeafViewTest()
    : ContinuousLagrangeSpaceTest<G, R, p>(make_cubic_grid<G>())
  {
  }
};


template <class G, class R, int p>
struct ContinuousLagrangeSpaceOnPrismLeafViewTest : public ContinuousLagrangeSpaceTest<G, R, p>
{
  ContinuousLagrangeSpaceOnPrismLeafViewTest()
    : ContinuousLagrangeSpaceTest<G, R, p>(make_prism_grid<G>())
  {
  }
};


template <class G, class R, int p>
struct ContinuousLagrangeSpaceOnMixedLeafViewTest : public ContinuousLagrangeSpaceTest<G, R, p>
{
  ContinuousLagrangeSpaceOnMixedLeafViewTest()
    : ContinuousLagrangeSpaceTest<G, R, p>(make_mixed_grid<G>())
  {
  }
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_SPACES_H1_CONTINUOUS_LAGRANGE_HH
