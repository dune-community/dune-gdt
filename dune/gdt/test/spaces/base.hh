// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014, 2016 - 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016 - 2018)

#ifndef DUNE_GDT_TEST_SPACES_BASE_HH
#define DUNE_GDT_TEST_SPACES_BASE_HH

#include <memory>
#include <type_traits>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


/// \todo Implement lagrange_basis_jacobians_seem_to_be_correct in the vector-valued case!
template <class Space, int p>
struct SpaceTestBase : public ::testing::Test
{
  static_assert(is_space<Space>::value, "");
  using SpaceType = Space;

  using GridViewType = typename SpaceType::GridViewType;
  using D = typename SpaceType::D;
  static const constexpr size_t d = SpaceType::d;
  using R = typename SpaceType::R;
  static const constexpr size_t r = SpaceType::r;
  static const constexpr size_t rC = SpaceType::rC;

  using GlobalBasisType = typename SpaceType::GlobalBasisType;
  using MapperType = typename SpaceType::MapperType;
  using LocalFiniteElementFamilyType = typename SpaceType::LocalFiniteElementFamilyType;

  std::shared_ptr<GridViewType> grid_view;
  std::shared_ptr<SpaceType> space;

  ~SpaceTestBase() = default;

  void SetUp() override final
  {
    ASSERT_NE(grid_view, nullptr) << "Any derived test has to create the grid_view on construction!";
    space = std::shared_ptr<SpaceType>(new SpaceType(*grid_view, p));
  }

  void TearDown() override final
  {
    space.reset();
  }

  void basis_of_lagrange_space_exists_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& element : elements(*grid_view))
      EXPECT_EQ(r * numLagrangePoints(element.geometry().type().id(), d, p), space->basis().localize(element)->size());
  }

  void basis_of_lagrange_space_exists_on_each_element_with_correct_order()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& element : elements(*grid_view))
      EXPECT_EQ(p, space->basis().localize(element)->order());
  }

  void lagrange_points_exist_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& geometry_type : grid_view->indexSet().types(0))
      EXPECT_EQ(numLagrangePoints(geometry_type.id(), d, p),
                space->finite_elements().get(geometry_type, p).lagrange_points().size());
  }

  void basis_is_lagrange_basis(const double& tolerance = 1e-15)
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& element : elements(*grid_view)) {
      const auto basis = space->basis().localize(element);
      const auto lagrange_points = space->finite_elements().get(element.geometry().type(), p).lagrange_points();
      EXPECT_EQ(lagrange_points.size(), basis->size() / r);
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto values = basis->evaluate_set(lagrange_points[ii]);
        for (size_t rr = 0; rr < r; ++rr) {
          for (size_t jj = 0; jj < lagrange_points.size(); ++jj) {
            EXPECT_TRUE(XT::Common::FloatCmp::eq(
                values[rr * lagrange_points.size() + jj][rr], ii == jj ? 1. : 0., tolerance, tolerance))
                << "ii = " << ii << "\nrr = " << rr << "\njj = " << jj
                << "\nlagrange_points[ii] = " << lagrange_points[ii]
                << "\nbasis->evaluate_set(lagrange_points[ii])[jj] = " << values[jj];
          }
        }
      }
    }
  } // ... basis_is_lagrange_basis(...)

  // I am too lazy to implement this in the vector-valued case.
  template <size_t r_ = r, typename = typename std::enable_if<r_ == r && r_ == 1, void>::type>
  void basis_jacobians_of_lagrange_space_seem_to_be_correct()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& element : elements(*grid_view)) {
      const auto& reference_element = ReferenceElements<D, d>::general(element.geometry().type());
      const auto basis = space->basis().localize(element);
      const double h = 1e-6;
      for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.geometry().type(), basis->order())) {
        const auto& xx = quadrature_point.position();
        const auto& J_inv_T = element.geometry().jacobianInverseTransposed(xx);
        const auto jacobians = basis->jacobians_of_set(xx);
        EXPECT_EQ(basis->size(), jacobians.size());
        const auto values_xx = basis->evaluate_set(xx);
        EXPECT_EQ(basis->size(), values_xx.size());
        auto approximate_jacobians = jacobians;
        // compute approximate partial derivatives
        for (size_t dd = 0; dd < d; ++dd) {
          // try to find a suitable x + h
          auto xx_plus_h = xx;
          xx_plus_h[dd] += h;
          if (!reference_element.checkInside(xx_plus_h)) {
            xx_plus_h[dd] -= 2. * h;
          }
          ASSERT_TRUE(reference_element.checkInside(xx_plus_h))
              << "xx_plus_h = " << xx_plus_h << " is not inside the reference element!";
          const auto values_xx_plus_h = basis->evaluate_set(xx_plus_h);
          EXPECT_EQ(basis->size(), values_xx_plus_h.size());
          for (size_t ii = 0; ii < basis->size(); ++ii) {
            approximate_jacobians[ii][0][dd] = (values_xx_plus_h[ii] - values_xx[ii]) / (xx_plus_h[dd] - xx[dd]);
            if (xx_plus_h[dd] - xx[dd] < 0)
              approximate_jacobians[ii][0][dd] *= -1.;
          }
        }
        // transform
        auto tmp_jac = approximate_jacobians[0][0];
        for (size_t ii = 0; ii < basis->size(); ++ii) {
          J_inv_T.mv(approximate_jacobians[ii][0], tmp_jac);
          approximate_jacobians[ii][0] = tmp_jac;
        }
        // check
        double tolerance = 1e-4;
        for (size_t ii = 0; ii < basis->size(); ++ii)
          EXPECT_TRUE(XT::Common::FloatCmp::eq(jacobians[ii][0], approximate_jacobians[ii][0], tolerance, tolerance))
              << "ii = " << ii << "\njacobians[ii][0] = " << jacobians[ii][0] << "\n"
              << "approximate_jacobians[ii][0] = " << approximate_jacobians[ii][0] << "\n"
              << "absolue L_infty error: " << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() << "\n"
              << "relative L_infty error: "
              << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() / jacobians[ii][0].infinity_norm();
      }
    }
  } // ... basis_jacobians_of_lagrange_space_seem_to_be_correct(...)

  void mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(space->is_lagrangian()) << "Do not call this test otherwise!";
    for (auto&& element : elements(*grid_view))
      EXPECT_EQ(r * numLagrangePoints(element.geometry().type().id(), d, p), space->mapper().local_size(element));
  }

  void mapper_reports_correct_max_num_DoFs()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    size_t max_num_dofs = 0;
    for (auto&& element : elements(*grid_view))
      max_num_dofs = std::max(max_num_dofs, space->mapper().local_size(element));
    EXPECT_LE(max_num_dofs, space->mapper().max_local_size());
  }

  void mapper_of_discontinuous_space_maps_correctly()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    // we want to check that the numbering is consecutive and that each global index exists only once
    std::set<size_t> global_indices;
    // we test both call variants
    std::set<size_t> map_to_global;
    for (auto&& element : elements(*grid_view)) {
      for (const auto& global_index : space->mapper().global_indices(element))
        global_indices.insert(global_index);
      for (size_t ii = 0; ii < space->mapper().local_size(element); ++ii)
        map_to_global.insert(space->mapper().global_index(element, ii));
    }
    // check for consecutive numbering
    EXPECT_EQ(0, *global_indices.begin());
    EXPECT_EQ(global_indices.size() - 1, *global_indices.rbegin());
    EXPECT_EQ(0, *map_to_global.begin());
    EXPECT_EQ(map_to_global.size() - 1, *map_to_global.rbegin());
    // check that the mapper is of the same opinion
    EXPECT_EQ(space->mapper().size(), global_indices.size());
    EXPECT_EQ(space->mapper().size(), map_to_global.size());
    // check that each global index is unique
    for (const auto& global_index : global_indices)
      EXPECT_EQ(1, global_indices.count(global_index));
    for (const auto& global_index : map_to_global)
      EXPECT_EQ(1, map_to_global.count(global_index));
  } // ... mapper_of_discontinuous_space_maps_correctly(...)

  void local_interpolation_seems_to_be_correct()
  {
    ASSERT_NE(grid_view, nullptr);
    ASSERT_NE(space, nullptr);
    for (const auto& geometry_type : grid_view->indexSet().types(0)) {
      const auto& finite_element = space->finite_elements().get(geometry_type, p);
      const auto& shape_functions = finite_element.basis();
      ASSERT_EQ(finite_element.size(), shape_functions.size());
      ASSERT_EQ(finite_element.size(), finite_element.interpolation().size());
      for (size_t ii = 0; ii < shape_functions.size(); ++ii) {
        const auto dofs = finite_element.interpolation().interpolate(
            [&](const auto& x) { return shape_functions.evaluate(x)[ii]; }, shape_functions.order());
        ASSERT_GE(dofs.size(), shape_functions.size());
        for (size_t jj = 0; jj < shape_functions.size(); ++jj)
          EXPECT_TRUE(XT::Common::FloatCmp::eq(ii == jj ? 1. : 0., dofs[jj]))
              << "\nii == jj ? 1. : 0. = " << (ii == jj ? 1. : 0.) << "\ndofs[jj] = " << dofs[jj];
      }
    }
  } // ... local_interpolation_seems_to_be_correct(...)
}; // struct SpaceTestBase


template <class G>
XT::Grid::GridProvider<G> make_simplicial_grid()
{ //                                                  (i) Negative coordinates and not the same as the reference element
  auto grid = XT::Grid::make_cube_grid<G>(-1.5, -1, 3); // (ii) at least 3 elements to have fully inner ones,
  grid.global_refine(1); //                              (iii) refine at least once to obtain all kinds of orientations.
  return grid;
}


using SimplicialGrids = ::testing::Types<ONED_1D,
                                         YASP_1D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_2D_SIMPLEX_CONFORMING,
                                         ALU_2D_SIMPLEX_NONCONFORMING
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_2D
#endif
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_3D_SIMPLEX_CONFORMING,
                                         ALU_3D_SIMPLEX_NONCONFORMING
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_3D
#endif
                                         >;


template <class G>
XT::Grid::GridProvider<G> make_cubic_grid()
{
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;
  FieldVector<D, d> lower_left(-1.5); //        (i) Negative coordinates and not the same as the reference element
  FieldVector<D, d> upper_right(-1.);
  std::array<unsigned int, d> num_elements; // (ii) at least 3 elements to have fully inner ones.
  std::fill(num_elements.begin(), num_elements.end(), 3);
  XT::Grid::GridProvider<G> grid(StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
  return grid;
} // ... make_cubic_grid(...)


using CubicGrids = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_2D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                    ,
                                    UG_2D
#endif
                                    ,
                                    YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_3D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                    ,
                                    UG_3D
#endif
                                    ,
                                    YASP_4D_EQUIDISTANT_OFFSET>;


template <class G,
          typename this_disables_dimensions_without_prisms =
              typename std::enable_if<XT::Grid::is_grid<G>::value && G::dimension == 3, void>::type>
XT::Grid::GridProvider<G> make_prism_grid()
{
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;
  GridFactory<G> factory;
  for (auto&& vertex : {XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                        XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                        XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                        XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                        XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                        XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.})}) {
    factory.insertVertex(vertex);
  }
  factory.insertElement(GeometryTypes::prism, {0, 1, 2, 3, 4, 5});
  XT::Grid::GridProvider<G> grid(factory.createGrid());
  grid.global_refine(1);
  return grid;
} // ... make_prism_grid(...)


using PrismGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D
#endif
    >;


template <
    class G,
    typename this_disables_dimensions_whithout_mixed_elements =
        typename std::enable_if<XT::Grid::is_grid<G>::value && (G::dimension == 2 || G::dimension == 3), void>::type>
XT::Grid::GridProvider<G> make_mixed_grid()
{
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;
  if (d == 2) {
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
    factory.insertElement(GeometryTypes::cube(2), {3, 0, 4, 1});
    factory.insertElement(GeometryTypes::cube(2), {4, 1, 5, 2});
    factory.insertElement(GeometryTypes::simplex(2), {4, 6, 3});
    factory.insertElement(GeometryTypes::simplex(2), {4, 5, 6});
    XT::Grid::GridProvider<G> grid(factory.createGrid());
    grid.global_refine(1);
    return grid;
  } else if (d == 3) {
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
    factory.insertElement(GeometryTypes::cube(3), {3, 0, 4, 1, 9, 6, 10, 7});
    factory.insertElement(GeometryTypes::cube(3), {4, 1, 5, 2, 10, 7, 11, 8});
    factory.insertElement(GeometryTypes::simplex(3), {4, 12, 3, 10});
    factory.insertElement(GeometryTypes::simplex(3), {4, 5, 12, 10});
    XT::Grid::GridProvider<G> grid(factory.createGrid());
    grid.global_refine(1);
    return grid;
  } else
    DUNE_THROW(InvalidStateException, "");
} // ... make_mixed_grid(...)


using MixedGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D,
    UG_3D
#endif
    >;


using MixedGridsWithConformingIntersections = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D
#endif
    >;


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_SPACES_BASE_HH
