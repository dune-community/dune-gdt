// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <algorithm>
#include <memory>
#include <tuple>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/refinement.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/spaces/rt/default.hh>


template <class GridLayerType, int p>
struct RtSpace : public ::testing::Test
{
  static_assert(p == 0, "The space cannot handle higher orders (yet)!");
  using SpaceType = Dune::GDT::RaviartThomasSpace<GridLayerType, p>;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;

  virtual std::shared_ptr<GridLayerType> grid_layer() = 0;

  std::shared_ptr<SpaceType> space;

  ~RtSpace() = default;

  void SetUp() override final
  {
    ASSERT_NE(grid_layer(), nullptr);
    space = std::shared_ptr<SpaceType>(new SpaceType(*grid_layer()));
  }

  void TearDown() override final
  {
    space.reset();
  }

  void basis_exists_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_layer())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.geometry().type());
      EXPECT_EQ(reference_element.size(1), space->base_function_set(element).size());
    }
  }

  void basis_exists_on_each_element_with_correct_order()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_layer()))
      EXPECT_EQ(1, space->base_function_set(element).order());
  }

  void mapper_reports_correct_num_DoFs_on_each_element()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_layer())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.geometry().type());
      EXPECT_EQ(reference_element.size(1), space->mapper().numDofs(element));
    }
  }

  void mapper_reports_correct_max_num_DoFs()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    size_t max_num_dofs = 0;
    for (auto&& element : elements(*grid_layer()))
      max_num_dofs = std::max(max_num_dofs, space->mapper().numDofs(element));
    EXPECT_LE(max_num_dofs, space->mapper().maxNumDofs());
  }

  void mapper_maps_correctly()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    // collect all global ids that are associated with a global lagrange point
    std::map<Dune::FieldVector<D, d>, std::set<size_t>, Dune::XT::Common::VectorLess>
        global_intersection_centers_to_global_indices_map;
    for (auto&& element : elements(*grid_layer())) {
      const auto global_indices = space->mapper().globalIndices(element);
      EXPECT_LE(space->mapper().numDofs(element), global_indices.size());
      const auto intersection_to_DoF_map = space->local_DoF_indices(element);
      EXPECT_EQ(intersection_to_DoF_map.size(), space->mapper().numDofs(element));
      for (auto&& intersection : intersections(*grid_layer(), element)) {
        const auto intersection_index = intersection.indexInInside();
        const auto local_DoF_index = intersection_to_DoF_map[intersection_index];
        const auto global_index = space->mapper().mapToGlobal(element, local_DoF_index);
        EXPECT_EQ(global_indices[local_DoF_index], global_index);
        global_intersection_centers_to_global_indices_map[intersection.geometry().center()].insert(global_index);
      }
    }
    // check that all intersections have indeed one and only one global DoF id ...
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : global_intersection_centers_to_global_indices_map) {
      const auto global_DoF_indices_per_intersection = entry.second;
      EXPECT_EQ(global_DoF_indices_per_intersection.size(), 1);
      global_DoF_indices.insert(*(global_DoF_indices_per_intersection.begin()));
    }
    EXPECT_EQ(global_intersection_centers_to_global_indices_map.size(), global_DoF_indices.size());
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
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_layer())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.geometry().type());
      EXPECT_EQ(reference_element.size(1), space->local_DoF_indices(element).size());
    }
  }

  void basis_is_rt_basis()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    Dune::GDT::ZeroOrderScalarDiscontinuousMapper<GridLayerType> entity_indices(*grid_layer());
    for (auto&& element : elements(*grid_layer())) {
      const auto basis = space->base_function_set(element);
      const auto intersection_to_DoF_index_map = space->local_DoF_indices(element);
      for (auto&& intersection : intersections(*grid_layer(), element)) {
        const auto xx_in_element_coordinates = intersection.geometry().center();
        const auto xx_in_reference_element_coordinates = element.geometry().local(xx_in_element_coordinates);
        const auto xx_in_reference_intersection_coordinates =
            intersection.geometryInInside().local(xx_in_reference_element_coordinates);
        const auto normal = intersection.integrationOuterNormal(xx_in_reference_intersection_coordinates);
        const auto basis_values = basis.evaluate(xx_in_reference_element_coordinates);
        const auto intersection_index = intersection.indexInInside();
        const auto DoF_index = intersection_to_DoF_index_map[intersection_index];
        double switch_ = 1;
        if (intersection.neighbor()
            && entity_indices.mapToGlobal(element, 0) < entity_indices.mapToGlobal(intersection.outside(), 0))
          switch_ *= -1.;
        for (size_t ii = 0; ii < basis.size(); ++ii)
          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(
              (ii == DoF_index ? 1. : 0.) * switch_, basis_values[ii] * normal, 1e-14, 1e-14))
              << "ii = " << ii << "\nDoF_index = " << DoF_index << "\nii == DoF_index ? 1. : 0. "
              << (ii == DoF_index ? 1. : 0.) << "\nbasis_values[ii] * normal = " << basis_values[ii] * normal;
      }
    }
  } // ... basis_is_rt_basis(...)

  void basis_jacobians_seem_to_be_correct()
  {
    ASSERT_NE(grid_layer(), nullptr);
    ASSERT_NE(space, nullptr);
    ASSERT_TRUE(false) << "continue here";
    //    for (auto&& element : elements(*grid_layer())) {
    //      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.geometry().type());
    //      const auto basis = space->base_function_set(element);
    //      const double h = 1e-6;
    //      for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(),
    //      basis.order())) {
    //        const auto& xx = quadrature_point.position();
    //        const auto& J_inv_T = element.geometry().jacobianInverseTransposed(xx);
    //        const auto jacobians = basis.jacobian(xx);
    //        EXPECT_EQ(basis.size(), jacobians.size());
    //        const auto values_xx = basis.evaluate(xx);
    //        EXPECT_EQ(basis.size(), values_xx.size());
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
    //          const auto values_xx_plus_h = basis.evaluate(xx_plus_h);
    //          EXPECT_EQ(basis.size(), values_xx_plus_h.size());
    //          for (size_t ii = 0; ii < basis.size(); ++ii) {
    //            approximate_jacobians[ii][0][dd] = (values_xx_plus_h[ii] - values_xx[ii]) / (xx_plus_h[dd] - xx[dd]);
    //            if (xx_plus_h[dd] - xx[dd] < 0)
    //              approximate_jacobians[ii][0][dd] *= -1.;
    //          }
    //        }
    //        // transform
    //        auto tmp_jac = approximate_jacobians[0][0];
    //        for (size_t ii = 0; ii < basis.size(); ++ii) {
    //          J_inv_T.mv(approximate_jacobians[ii][0], tmp_jac);
    //          approximate_jacobians[ii][0] = tmp_jac;
    //        }
    //        // check
    //        for (size_t ii = 0; ii < basis.size(); ++ii)
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
}; // struct RtSpace


template <class G, int p>
struct RtSpaceOnSimplicialLeafView
  : public RtSpace<typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
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

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct RtSpaceOnSimplicialLeafView


using SimplicialGrids = ::testing::Types<YASP_1D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_2D_SIMPLEX_CONFORMING,
                                         ALU_2D_SIMPLEX_NONCONFORMING
#endif
// UG does not work until we have the virtual interfaces for the finite elements.
//#if HAVE_DUNE_UGGRID || HAVE_UG
//                                         ,
//                                         UG_2D
//#endif
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_3D_SIMPLEX_CONFORMING,
                                         ALU_3D_SIMPLEX_NONCONFORMING
#endif
                                         // s.a.
                                         //#if HAVE_DUNE_UGGRID || HAVE_UG
                                         //                                         ,
                                         //                                         UG_3D
                                         //#endif
                                         >;


template <class G>
using Order0SimplicialRtSpace = RtSpaceOnSimplicialLeafView<G, 0>;
TYPED_TEST_CASE(Order0SimplicialRtSpace, SimplicialGrids);
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


template <class G, int p>
struct RtSpaceOnCubicLeafView
  : public RtSpace<typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  RtSpaceOnCubicLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    Dune::FieldVector<D, d> lower_left(-1.5); //  (i) negative coordinates and not the same as the reference element
    Dune::FieldVector<D, d> upper_right(-1.);
    std::array<unsigned int, d> num_elements; // (ii) at least 3 elements to have fully inner ones
    std::fill(num_elements.begin(), num_elements.end(), 3);
    grid_provider = std::make_shared<GridProviderType>(
        Dune::StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
  }

  ~RtSpaceOnCubicLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct RtSpaceOnCubicLeafView


using CubicGrids = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_2D_CUBE
#endif
                                    //// s.a.
                                    //#if HAVE_DUNE_UGGRID || HAVE_UG
                                    //                                    ,
                                    //                                    UG_2D
                                    //#endif
                                    ,
                                    YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_3D_CUBE
#endif
                                    //// s.a.
                                    //#if HAVE_DUNE_UGGRID || HAVE_UG
                                    //                                    ,
                                    //                                    UG_3D
                                    //#endif
                                    >;

template <class G>
using Order0CubicRtSpace = RtSpaceOnCubicLeafView<G, 0>;
TYPED_TEST_CASE(Order0CubicRtSpace, CubicGrids);
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


//// The space cannot handle mixed views (yet)!
// template <class G, int p>
// struct RtSpaceOnMixedLeafView : public RtSpace<typename Dune::XT::Grid::GridProvider<G,
// Dune::XT::Grid::none_t>::LeafGridViewType, p>
//{
//  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
//  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

//  std::shared_ptr<GridProviderType> grid_provider;
//  std::shared_ptr<LeafGridViewType> leaf_view;

//  RtSpaceOnMixedLeafView()
//  {
//    using D = typename G::ctype;
//    static const constexpr size_t d = G::dimension;
//    switch (d) {
//      case 1: {
//        // cannot use ASSERT_... in a ctor
//        EXPECT_TRUE(false) << "Does not make sense in 1d (all cubes are simplices)!\n"
//                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
//        grid_provider = nullptr;
//        leaf_view = nullptr;
//        break;
//      }
//      case 2: {
//        Dune::GridFactory<G> factory;
//        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25})}) {
//          factory.insertVertex(vertex);
//        }
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {3, 0, 4, 1});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {4, 1, 5, 2});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 6, 3});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 5, 6});
//        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
//        grid_provider->global_refine(1);
//        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
//        break;
//      }
//      case 3: {
//        Dune::GridFactory<G> factory;
//        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.5}),
//                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25, -1.})}) {
//          factory.insertVertex(vertex);
//        }
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {3, 0, 4, 1, 9, 6, 10, 7});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {4, 1, 5, 2, 10, 7, 11, 8});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 12, 3, 10});
//        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 5, 12, 10});
//        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
//        grid_provider->global_refine(1);
//        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
//        break;
//      }
//      default: {
//        // cannot use ASSERT_... in a ctor
//        EXPECT_TRUE(false) << "Not implemented yet for dimension " << d << "!\n"
//                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
//        grid_provider = nullptr;
//        leaf_view = nullptr;
//      }
//    }
//  } // RtSpaceOnMixedLeafView(...)

//  ~RtSpaceOnMixedLeafView() = default;

//  std::shared_ptr<LeafGridViewType> grid_layer() override final
//  {
//    return leaf_view;
//  }
//}; // struct RtSpaceOnMixedLeafView


// using MixedGrids = ::testing::Types<
//#ifHAVE_DUNE_UGGRID || HAVE_UG
//    UG_2D,
//    UG_3D
//#endif
//    >;


// template <class G>
// using Order0MixedRtSpace = RtSpaceOnMixedLeafView<G, 0>;
// TYPED_TEST_CASE(Order0MixedRtSpace, MixedGrids);
// TYPED_TEST(Order0MixedRtSpace, basis_exists_on_each_element_with_correct_size)
//{
//  this->basis_exists_on_each_element_with_correct_size();
//}
// TYPED_TEST(Order0MixedRtSpace, basis_exists_on_each_element_with_correct_order)
//{
//  this->basis_exists_on_each_element_with_correct_order();
//}
// TYPED_TEST(Order0MixedRtSpace, mapper_reports_correct_num_DoFs_on_each_element)
//{
//  this->mapper_reports_correct_num_DoFs_on_each_element();
//}
// TYPED_TEST(Order0MixedRtSpace, mapper_reports_correct_max_num_DoFs)
//{
//  this->mapper_reports_correct_max_num_DoFs();
//}
// TYPED_TEST(Order0MixedRtSpace, mapper_maps_correctly)
//{
//  this->mapper_maps_correctly();
//}
// TYPED_TEST(Order0MixedRtSpace, local_DoF_indices_exist_on_each_element_with_correct_size)
//{
//  this->local_DoF_indices_exist_on_each_element_with_correct_size();
//}
// TYPED_TEST(Order0MixedRtSpace, basis_is_rt_basis)
//{
//  this->basis_is_rt_basis();
//}
// TYPED_TEST(Order0MixedRtSpace, basis_jacobians_seem_to_be_correct)
//{
//  this->basis_jacobians_seem_to_be_correct();
//}
