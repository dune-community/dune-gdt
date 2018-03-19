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

#include <memory>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/spaces/l2/discontinuous-galerkin.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class GridViewType, size_t r, int p>
struct DiscontinuousLagrangeSpaceTest : public SpaceTestBase<DiscontinuousLagrangeSpace<GridViewType, p, r>, p>
{
  using BaseType = SpaceTestBase<DiscontinuousLagrangeSpace<GridViewType, p, r>, p>;

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
struct DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest
    : public DiscontinuousLagrangeSpaceTest<typename XT::Grid::GridProvider<G>::LeafGridViewType, r, p>
{
  using GridProviderType = XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  GridProviderType grid_provider;

  DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest() // (i) negative coordinates and not the same as the reference
      : grid_provider(XT::Grid::make_cube_grid<G>(-1.5, -1, 3).grid_ptr()) //                                element,
  { //                                                   (ii) at least 3 elements to have fully inner ones,
    grid_provider.global_refine(1); //                  (iii) refine at least once to obtain all kinds of orientations
    this->grid_view = std::make_shared<LeafGridViewType>(grid_provider.leaf_view());
  }
};


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnCubicLeafViewTest
    : public DiscontinuousLagrangeSpaceTest<typename XT::Grid::GridProvider<G>::LeafGridViewType, r, p>
{
  using GridProviderType = XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;

  DiscontinuousLagrangeSpaceOnCubicLeafViewTest()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    FieldVector<D, d> lower_left(-1.5); //        (i) negative coordinates and not the same as the reference element
    FieldVector<D, d> upper_right(-1.);
    std::array<unsigned int, d> num_elements; // (ii) at least 3 elements to have fully inner ones
    std::fill(num_elements.begin(), num_elements.end(), 3);
    grid_provider = std::make_shared<GridProviderType>(
        StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
    this->grid_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
  }
}; // struct DiscontinuousLagrangeSpaceOnCubicLeafViewTest


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnPrismLeafViewTest
    : public DiscontinuousLagrangeSpaceTest<typename XT::Grid::GridProvider<G>::LeafGridViewType, r, p>
{
  using GridProviderType = XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;

  DiscontinuousLagrangeSpaceOnPrismLeafViewTest()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    if (d == 3) {
      GridFactory<G> factory;
      for (auto&& vertex : {XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                            XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                            XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                            XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                            XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                            XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.})}) {
        factory.insertVertex(vertex);
      }
      factory.insertElement(GeometryType(GeometryType::prism, 3), {0, 1, 2, 3, 4, 5});
      grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
      grid_provider->global_refine(1);
      this->grid_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
    } else {
      // cannot use ASSERT_... in a ctor
      EXPECT_TRUE(false) << "Does not make sense in " << d << "d!\n"
                         << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
      grid_provider = nullptr;
      this->grid_view = nullptr;
    }
  } // DiscontinuousLagrangeSpaceOnPrismLeafViewTest(...)
}; // struct DiscontinuousLagrangeSpaceOnPrismLeafViewTest


template <class G, size_t r, int p>
struct DiscontinuousLagrangeSpaceOnMixedLeafViewTest
    : public DiscontinuousLagrangeSpaceTest<typename XT::Grid::GridProvider<G>::LeafGridViewType, r, p>
{
  using GridProviderType = XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  DiscontinuousLagrangeSpaceOnMixedLeafViewTest()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    switch (d) {
      case 1: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Does not make sense in 1d (all cubes are simplices)!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        this->grid_view = nullptr;
        break;
      }
      case 2: {
        GridFactory<G> factory;
        for (auto&& vertex : {XT::Common::FieldVector<D, d>({-1., -1.5}),
                              XT::Common::FieldVector<D, d>({-1., -1.25}),
                              XT::Common::FieldVector<D, d>({-1., -1.}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.5}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.25}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.}),
                              XT::Common::FieldVector<D, d>({-1.75, -1.25})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(GeometryType(GeometryType::cube, 2), {3, 0, 4, 1});
        factory.insertElement(GeometryType(GeometryType::cube, 2), {4, 1, 5, 2});
        factory.insertElement(GeometryType(GeometryType::simplex, 2), {4, 6, 3});
        factory.insertElement(GeometryType(GeometryType::simplex, 2), {4, 5, 6});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        this->grid_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      case 3: {
        GridFactory<G> factory;
        for (auto&& vertex : {XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                              XT::Common::FieldVector<D, d>({-1., -1.25, -1.}),
                              XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.}),
                              XT::Common::FieldVector<D, d>({-1.5, -1., -1.}),
                              XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                              XT::Common::FieldVector<D, d>({-1., -1.25, -1.5}),
                              XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                              XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.5}),
                              XT::Common::FieldVector<D, d>({-1.5, -1., -1.5}),
                              XT::Common::FieldVector<D, d>({-1.75, -1.25, -1.})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(GeometryType(GeometryType::cube, 3), {3, 0, 4, 1, 9, 6, 10, 7});
        factory.insertElement(GeometryType(GeometryType::cube, 3), {4, 1, 5, 2, 10, 7, 11, 8});
        factory.insertElement(GeometryType(GeometryType::simplex, 3), {4, 12, 3, 10});
        factory.insertElement(GeometryType(GeometryType::simplex, 3), {4, 5, 12, 10});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        this->grid_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      default: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Not implemented yet for dimension " << d << "!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        this->grid_view = nullptr;
      }
    }
  } // DiscontinuousLagrangeSpaceOnMixedLeafViewTest(...)
}; // struct DiscontinuousLagrangeSpaceOnMixedLeafViewTest


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_SPACES_L2_DISCONTINUOUS_GALERKIN_HH
