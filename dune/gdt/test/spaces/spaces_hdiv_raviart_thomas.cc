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

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <algorithm>
#include <memory>
#include <tuple>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/refinement.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>


template <class GridViewType, int p>
struct RtSpace : public ::testing::Test
{
  static_assert(p == 0, "The space cannot handle higher orders (yet)!");
  using SpaceType = Dune::GDT::RaviartThomasSpace<GridViewType, double>;
  using D = typename SpaceType::D;
  static const constexpr size_t d = SpaceType::d;

  virtual std::shared_ptr<GridViewType> grid_view() = 0;

  std::shared_ptr<SpaceType> space;

  ~RtSpace() = default;

  void SetUp() override final
  {
    ASSERT_NE(grid_view(), nullptr);
    space = std::shared_ptr<SpaceType>(new SpaceType(*grid_view(), p));
  }

  void TearDown() override final
  {
    space.reset();
  }

  void gives_correct_identification()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    EXPECT_EQ(Dune::GDT::SpaceType::raviart_thomas, space->type());
    EXPECT_EQ(p, space->min_polorder());
    EXPECT_EQ(p, space->max_polorder());
    for (int diff_order : {0, 1, 2, 3, 4, 5})
      EXPECT_FALSE(space->continuous(diff_order));
    EXPECT_TRUE(space->continuous_normal_components());
    EXPECT_FALSE(space->is_lagrangian());
  }

  void basis_exists_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.type());
      EXPECT_EQ(reference_element.size(1), space->basis().localize(element)->size());
    }
  }

  void basis_exists_on_each_element_with_correct_order()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view()))
      EXPECT_EQ(1, space->basis().localize(element)->order());
  }

  void mapper_reports_correct_num_DoFs_on_each_element()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.type());
      EXPECT_EQ(reference_element.size(1), space->mapper().local_size(element));
    }
  }

  void mapper_reports_correct_max_num_DoFs()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    size_t max_num_dofs = 0;
    for (auto&& element : elements(*grid_view()))
      max_num_dofs = std::max(max_num_dofs, space->mapper().local_size(element));
    EXPECT_LE(max_num_dofs, space->mapper().max_local_size());
  }

  // this only works for conforming intersections!
  void mapper_maps_correctly()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    Dune::GDT::FiniteVolumeMapper<GridViewType> entity_indices(*grid_view());
    // determine global numbering of intersections
    std::map<size_t, std::map<size_t, size_t>> element_and_local_intersection_index_to_global_intersection_index;
    size_t tmp_counter = 0;
    for (auto&& element : elements(*grid_view())) {
      const auto global_element_index = entity_indices.global_index(element, 0);
      for (auto&& intersection : intersections(*grid_view(), element)) {
        const auto element_local_intersection_index = intersection.indexInInside();
        if (!intersection.neighbor()
            || (entity_indices.global_index(element, 0) < entity_indices.global_index(intersection.outside(), 0))) {
          // this is the place to handle this intersection
          const auto global_intersection_index = tmp_counter;
          ++tmp_counter;
          // store the global intersection index for this element
          element_and_local_intersection_index_to_global_intersection_index[global_element_index]
                                                                           [element_local_intersection_index] =
                                                                               global_intersection_index;
          if (intersection.neighbor()) {
            // as well as for the neighbor
            const auto global_neighbor_index = entity_indices.global_index(intersection.outside(), 0);
            const auto neighbor_local_intersection_index = intersection.indexInOutside();
            element_and_local_intersection_index_to_global_intersection_index[global_neighbor_index]
                                                                             [neighbor_local_intersection_index] =
                                                                                 global_intersection_index;
          }
        }
      }
    }
    // collect all global ids that are associated with an intersection
    std::map<size_t, std::set<size_t>> global_intersection_index_to_global_indices_map;
    for (auto&& element : elements(*grid_view())) {
      const auto global_DoF_indices = space->mapper().global_indices(element);
      EXPECT_LE(space->mapper().local_size(element), global_DoF_indices.size());
      const auto intersection_to_local_DoF_indices_map =
          space->finite_elements().get(element.type(), p).coefficients().local_key_indices(1);
      EXPECT_EQ(intersection_to_local_DoF_indices_map.size(), space->mapper().local_size(element));
      for (auto&& intersection : intersections(*grid_view(), element)) {
        const auto local_intersection_index = intersection.indexInInside();
        const auto local_DoF_indices = intersection_to_local_DoF_indices_map[local_intersection_index];
        EXPECT_EQ(1, local_DoF_indices.size());
        const auto local_DoF_index = *local_DoF_indices.begin();
        const auto global_DoF_index = space->mapper().global_index(element, local_DoF_index);
        EXPECT_EQ(global_DoF_indices[local_DoF_index], global_DoF_index);
        const auto global_intersection_index = element_and_local_intersection_index_to_global_intersection_index
                                                   .at(entity_indices.global_index(element, 0))
                                                   .at(local_intersection_index);
        global_intersection_index_to_global_indices_map[global_intersection_index].insert(global_DoF_index);
      }
    }
    // check that all intersections have indeed one and only one global DoF id ...
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : global_intersection_index_to_global_indices_map) {
      const auto global_DoF_indices_per_intersection = entry.second;
      EXPECT_EQ(global_DoF_indices_per_intersection.size(), 1);
      global_DoF_indices.insert(*(global_DoF_indices_per_intersection.begin()));
    }
    EXPECT_EQ(global_intersection_index_to_global_indices_map.size(), global_DoF_indices.size());
    // ... and that the numbering is consecutive
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      EXPECT_EQ(global_DoF_id, count);
      ++count;
    }
    EXPECT_EQ(global_DoF_indices.size(), space->mapper().size());
  } // ... mapper_maps_correctly(...)

  void local_DoF_indices_exist_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.type());
      const auto intersection_to_local_DoF_indices_map =
          space->finite_elements().get(element.type(), p).coefficients().local_key_indices(1);
      EXPECT_EQ(reference_element.size(1), intersection_to_local_DoF_indices_map.size());
      for (const auto& DoF_indices_per_intersection : intersection_to_local_DoF_indices_map)
        EXPECT_EQ(1, DoF_indices_per_intersection.size());
    }
  }

  void basis_is_rt_basis()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    Dune::GDT::FiniteVolumeMapper<GridViewType> entity_indices(*grid_view());
    for (auto&& element : elements(*grid_view())) {
      const auto basis = space->basis().localize(element);
      const auto intersection_to_local_DoF_indices_map =
          space->finite_elements().get(element.type(), p).coefficients().local_key_indices(1);
      for (auto&& intersection : intersections(*grid_view(), element)) {
        const auto xx_in_element_coordinates = intersection.geometry().center();
        const auto xx_in_reference_element_coordinates = element.geometry().local(xx_in_element_coordinates);
        const auto xx_in_reference_intersection_coordinates =
            intersection.geometryInInside().local(xx_in_reference_element_coordinates);
        const auto normal = intersection.integrationOuterNormal(xx_in_reference_intersection_coordinates);
        const auto basis_values = basis->evaluate_set(xx_in_reference_element_coordinates);
        const auto intersection_index = intersection.indexInInside();
        for (const auto& DoF_index : intersection_to_local_DoF_indices_map[intersection_index]) {
          double switch_ = 1;
          if (intersection.neighbor()
              && entity_indices.global_index(element, 0) < entity_indices.global_index(intersection.outside(), 0))
            switch_ *= -1.;
          for (size_t ii = 0; ii < basis->size(); ++ii)
            EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(
                (ii == DoF_index ? 1. : 0.) * switch_, basis_values[ii] * normal, 1e-14, 1e-14))
                << "ii = " << ii << "\nDoF_index = " << DoF_index << "\nii == DoF_index ? 1. : 0. "
                << (ii == DoF_index ? 1. : 0.) << "\nswitch_ = " << switch_
                << "\nbasis_values[ii] * normal = " << basis_values[ii] * normal;
        }
      }
    }
  } // ... basis_is_rt_basis(...)

  void basis_jacobians_seem_to_be_correct()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(false) << "continue here";
    //    for (auto&& element : elements(*grid_view())) {
    //      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.type());
    //      const auto basis = space->base_function_set(element);
    //      const double h = 1e-6;
    //      for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.type(),
    //      basis->order())) {
    //        const auto& xx = quadrature_point.position();
    //        const auto& J_inv_T = element.geometry().jacobianInverseTransposed(xx);
    //        const auto jacobians = basis->jacobians_of_set(xx);
    //        EXPECT_EQ(basis->size(), jacobians.size());
    //        const auto values_xx = basis->evaluate_set(xx);
    //        EXPECT_EQ(basis->size(), values_xx.size());
    //        auto approximate_jacobians = jacobians;
    //        // compute approximate partial derivatives
    //        for (size_t dd = 0; dd < d; ++dd) {
    //          // try to find a suitable x + h
    //          auto xx_plus_h = xx;
    //          xx_plus_h[dd] += h;
    //          if (!reference_element.checkInside(xx_plus_h)) {
    //            xx_plus_h[dd] -= 2. * h;
    //          }
    //          ASSERT_TRUE(reference_element.checkInside(xx_plus_h)) << "xx_plus_h = " << xx_plus_h
    //                                                                << " is not inside the reference element!";
    //          const auto values_xx_plus_h = basis->evaluate_set(xx_plus_h);
    //          EXPECT_EQ(basis->size(), values_xx_plus_h.size());
    //          for (size_t ii = 0; ii < basis->size(); ++ii) {
    //            approximate_jacobians[ii][0][dd] = (values_xx_plus_h[ii] - values_xx[ii]) / (xx_plus_h[dd] - xx[dd]);
    //            if (xx_plus_h[dd] - xx[dd] < 0)
    //              approximate_jacobians[ii][0][dd] *= -1.;
    //          }
    //        }
    //        // transform
    //        auto tmp_jac = approximate_jacobians[0][0];
    //        for (size_t ii = 0; ii < basis->size(); ++ii) {
    //          J_inv_T.mv(approximate_jacobians[ii][0], tmp_jac);
    //          approximate_jacobians[ii][0] = tmp_jac;
    //        }
    //        // check
    //        for (size_t ii = 0; ii < basis->size(); ++ii)
    //          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(jacobians[ii][0], approximate_jacobians[ii][0], 1e-4, 1e-4))
    //              << "ii = " << ii << "\njacobians[ii][0] = " << jacobians[ii][0] << "\n"
    //              << "approximate_jacobians[ii][0] = " << approximate_jacobians[ii][0] << "\n"
    //              << "absolue L_infty error: " << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() <<
    //              "\n"
    //              << "relative L_infty error: "
    //              << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() /
    //              jacobians[ii][0].infinity_norm();
    //      }
    //    }
  } // ... basis_jacobians_seem_to_be_correct(...)

  void local_interpolation_seems_to_be_correct()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (const auto& geometry_type : grid_view()->indexSet().types(0)) {
      const auto& finite_element = space->finite_elements().get(geometry_type, p);
      const auto& shape_functions = finite_element.basis();
      ASSERT_EQ(finite_element.size(), shape_functions.size());
      ASSERT_EQ(finite_element.size(), finite_element.interpolation().size());
      for (size_t ii = 0; ii < shape_functions.size(); ++ii) {
        const auto dofs = finite_element.interpolation().interpolate(
            [&](const auto& x) { return shape_functions.evaluate(x)[ii]; }, shape_functions.order());
        ASSERT_GE(dofs.size(), shape_functions.size());
        for (size_t jj = 0; jj < shape_functions.size(); ++jj)
          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(ii == jj ? 1. : 0., dofs[jj]))
              << "\nii == jj ? 1. : 0. = " << (ii == jj ? 1. : 0.) << "\ndofs[jj] = " << dofs[jj];
      }
    }
  } // ... local_interpolation_seems_to_be_correct(...)
}; // struct RtSpace


template <class G, int p>
struct RtSpaceOnSimplicialLeafView : public RtSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  GridProviderType grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  RtSpaceOnSimplicialLeafView() //        (i) negative coordinates and not the same as the reference
    : grid_provider(Dune::XT::Grid::make_cube_grid<G>(-1.5, -1, 3).grid_ptr()) //                          element,
  { //                                                   (ii) at least 3 elements to have fully inner ones,
    grid_provider.global_refine(1); //                  (iii) refine at least once to obtain all kinds of orientations
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider.leaf_view());
  }

  ~RtSpaceOnSimplicialLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct RtSpaceOnSimplicialLeafView


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
using Order0SimplicialRtSpace = RtSpaceOnSimplicialLeafView<G, 0>;
TYPED_TEST_CASE(Order0SimplicialRtSpace, SimplicialGrids);
TYPED_TEST(Order0SimplicialRtSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0SimplicialRtSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0SimplicialRtSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0SimplicialRtSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0SimplicialRtSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0SimplicialRtSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0SimplicialRtSpace, local_DoF_indices_exist_on_each_element_with_correct_size)
{
  this->local_DoF_indices_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0SimplicialRtSpace, basis_is_rt_basis)
{
  this->basis_is_rt_basis();
}
// TYPED_TEST(Order0SimplicialRtSpace, basis_jacobians_seem_to_be_correct)
//{
//  this->basis_jacobians_seem_to_be_correct();
//}
TYPED_TEST(Order0SimplicialRtSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G, int p>
struct RtSpaceOnCubicLeafView : public RtSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;
  using BaseType = RtSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>;
  using BaseType::d;
  using typename BaseType::D;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  RtSpaceOnCubicLeafView()
  {
    Dune::FieldVector<D, d> lower_left(-1.5); //  (i) negative coordinates and not the same as the reference element
    Dune::FieldVector<D, d> upper_right(-1.);
    std::array<unsigned int, d> num_elements; // (ii) at least 3 elements to have fully inner ones
    std::fill(num_elements.begin(), num_elements.end(), 3);
    grid_provider = std::make_shared<GridProviderType>(
        Dune::StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
  }

  ~RtSpaceOnCubicLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct RtSpaceOnCubicLeafView


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
                                    >;

template <class G>
using Order0CubicRtSpace = RtSpaceOnCubicLeafView<G, 0>;
TYPED_TEST_CASE(Order0CubicRtSpace, CubicGrids);
TYPED_TEST(Order0CubicRtSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0CubicRtSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0CubicRtSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0CubicRtSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0CubicRtSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0CubicRtSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0CubicRtSpace, local_DoF_indices_exist_on_each_element_with_correct_size)
{
  this->local_DoF_indices_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0CubicRtSpace, basis_is_rt_basis)
{
  this->basis_is_rt_basis();
}
// TYPED_TEST(Order0CubicRtSpace, basis_jacobians_seem_to_be_correct)
//{
//  this->basis_jacobians_seem_to_be_correct();
//}
TYPED_TEST(Order0CubicRtSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G, int p>
struct RtSpaceOnMixedLeafView : public RtSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;
  using BaseType = RtSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>;
  using BaseType::d;
  using typename BaseType::D;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  RtSpaceOnMixedLeafView()
  {
    switch (d) {
      case 1: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Does not make sense in 1d (all cubes are simplices)!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        leaf_view = nullptr;
        break;
      }
      case 2: {
        Dune::GridFactory<G> factory;
        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(Dune::GeometryTypes::cube(2), {3, 0, 4, 1});
        factory.insertElement(Dune::GeometryTypes::cube(2), {4, 1, 5, 2});
        factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 6, 3});
        factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 5, 6});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      case 3: {
        Dune::GridFactory<G> factory;
        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25, -1.})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(Dune::GeometryTypes::cube(3), {3, 0, 4, 1, 9, 6, 10, 7});
        factory.insertElement(Dune::GeometryTypes::cube(3), {4, 1, 5, 2, 10, 7, 11, 8});
        factory.insertElement(Dune::GeometryTypes::simplex(3), {4, 12, 3, 10});
        factory.insertElement(Dune::GeometryTypes::simplex(3), {4, 5, 12, 10});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      default: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Not implemented yet for dimension " << d << "!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        leaf_view = nullptr;
      }
    }
  } // RtSpaceOnMixedLeafView(...)

  ~RtSpaceOnMixedLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct RtSpaceOnMixedLeafView


using MixedGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D
//,
//    UG_3D // Intersections are non-conforming between simplices and cubes in 3d, which the mapper cannot handle yet!
#endif
    >;


template <class G>
using Order0MixedRtSpace = RtSpaceOnMixedLeafView<G, 0>;
TYPED_TEST_CASE(Order0MixedRtSpace, MixedGrids);
TYPED_TEST(Order0MixedRtSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0MixedRtSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0MixedRtSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0MixedRtSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0MixedRtSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0MixedRtSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0MixedRtSpace, local_DoF_indices_exist_on_each_element_with_correct_size)
{
  this->local_DoF_indices_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0MixedRtSpace, basis_is_rt_basis)
{
  this->basis_is_rt_basis();
}
// TYPED_TEST(Order0MixedRtSpace, basis_jacobians_seem_to_be_correct)
//{
//  this->basis_jacobians_seem_to_be_correct();
//}
TYPED_TEST(Order0MixedRtSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}
