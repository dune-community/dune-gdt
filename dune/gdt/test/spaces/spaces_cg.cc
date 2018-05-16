// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

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

#include <dune/gdt/spaces/cg.hh>


template <class GridLayerType, int p>
struct ContinuousLagrangeSpace : public ::testing::Test
{
  using SpaceType = Dune::GDT::ContinuousLagrangeSpace<GridLayerType, p>;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;

  virtual std::shared_ptr<GridLayerType> grid_layer() = 0;

  std::shared_ptr<SpaceType> space;

  ~ContinuousLagrangeSpace() = default;

  void SetUp() override final
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    space = std::shared_ptr<SpaceType>(new SpaceType(*grid_layer()));
  }

  void TearDown() override final
  {
    space.reset();
  }

  void basis_exists_on_each_element_with_correct_size()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer()))
      EXPECT_EQ(Dune::numLagrangePoints(element.geometry().type().id(), d, p),
                space->base_function_set(element).size());
  }

  void basis_exists_on_each_element_with_correct_order()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer()))
      EXPECT_EQ(p, space->base_function_set(element).order());
  }

  void mapper_reports_correct_num_DoFs_on_each_element()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer()))
      EXPECT_EQ(Dune::numLagrangePoints(element.geometry().type().id(), d, p), space->mapper().numDofs(element));
  }

  void mapper_reports_correct_max_num_DoFs()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    size_t max_num_dofs = 0;
    for (auto&& element : elements(*grid_layer()))
      max_num_dofs = std::max(max_num_dofs, space->mapper().numDofs(element));
    EXPECT_LE(max_num_dofs, space->mapper().maxNumDofs());
  }

  void mapper_maps_correctly()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    // collect all global ids that are associated with a global lagrange point
    std::map<Dune::FieldVector<D, d>, std::set<size_t>, Dune::XT::Common::FieldVectorLess>
        global_lagrange_point_to_global_indices_map;
    for (auto&& element : elements(*grid_layer())) {
      const auto global_indices = space->mapper().globalIndices(element);
      EXPECT_LE(space->mapper().numDofs(element), global_indices.size());
      const auto lagrange_points = space->lagrange_points(element);
      EXPECT_EQ(lagrange_points.size(), space->mapper().numDofs(element));
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto global_lagrange_point = element.geometry().global(lagrange_points[ii]);
        const auto global_index = space->mapper().mapToGlobal(element, ii);
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
    EXPECT_EQ(global_DoF_indices.size(), space->mapper().size());
  } // ... mapper_maps_correctly(...)

  void lagrange_points_exist_on_each_element_with_correct_size()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer()))
      EXPECT_EQ(Dune::numLagrangePoints(element.geometry().type().id(), d, p), space->lagrange_points(element).size());
  }

  void basis_is_lagrange_basis()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer())) {
      const auto basis = space->base_function_set(element);
      const auto lagrange_points = space->lagrange_points(element);
      EXPECT_EQ(lagrange_points.size(), basis.size());
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto values = basis.evaluate(lagrange_points[ii]);
        for (size_t jj = 0; jj < basis.size(); ++jj) {
          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(values[jj][0], ii == jj ? 1. : 0.))
              << "lagrange_points[" << ii << "] = " << lagrange_points[ii]
              << "\nbasis.evaluate(lagrange_points[ii]) = " << values;
        }
      }
    }
  } // ... basis_is_lagrange_basis(...)

  void basis_jacobians_seem_to_be_correct()
  {
    ASSERT_TRUE(grid_layer() != nullptr && grid_layer() != 0);
    ASSERT_TRUE(space != nullptr && space != 0);
    for (auto&& element : elements(*grid_layer())) {
      const auto& reference_element = Dune::ReferenceElements<D, d>::general(element.geometry().type());
      const auto basis = space->base_function_set(element);
      const double h = 1e-6;
      assert(basis.order() <= std::numeric_limits<int>::max());
      for (const auto& quadrature_point :
           Dune::QuadratureRules<D, d>::rule(element.geometry().type(), static_cast<int>(basis.order()))) {
        const auto& xx = quadrature_point.position();
        const auto& J_inv_T = element.geometry().jacobianInverseTransposed(xx);
        const auto jacobians = basis.jacobian(xx);
        EXPECT_EQ(basis.size(), jacobians.size());
        const auto values_xx = basis.evaluate(xx);
        EXPECT_EQ(basis.size(), values_xx.size());
        auto approximate_jacobians = jacobians;
        // compute approximate partial derivatives
        for (size_t dd = 0; dd < d; ++dd) {
          // try to find a suitable x + h
          auto xx_plus_h = xx;
          xx_plus_h[dd] += h;
          if (!reference_element.checkInside(xx_plus_h)) {
            xx_plus_h[dd] -= 2. * h;
          }
          ASSERT_TRUE(reference_element.checkInside(xx_plus_h)) << "xx_plus_h = " << xx_plus_h
                                                                << " is not inside the reference element!";
          const auto values_xx_plus_h = basis.evaluate(xx_plus_h);
          EXPECT_EQ(basis.size(), values_xx_plus_h.size());
          for (size_t ii = 0; ii < basis.size(); ++ii) {
            approximate_jacobians[ii][0][dd] = (values_xx_plus_h[ii] - values_xx[ii]) / (xx_plus_h[dd] - xx[dd]);
            if (xx_plus_h[dd] - xx[dd] < 0)
              approximate_jacobians[ii][0][dd] *= -1.;
          }
        }
        // transform
        auto tmp_jac = approximate_jacobians[0][0];
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          J_inv_T.mv(approximate_jacobians[ii][0], tmp_jac);
          approximate_jacobians[ii][0] = tmp_jac;
        }
        // check
        for (size_t ii = 0; ii < basis.size(); ++ii)
          EXPECT_TRUE(Dune::XT::Common::FloatCmp::eq(jacobians[ii][0], approximate_jacobians[ii][0], 1e-4, 1e-4))
              << "ii = " << ii << "\njacobians[ii][0] = " << jacobians[ii][0] << "\n"
              << "approximate_jacobians[ii][0] = " << approximate_jacobians[ii][0] << "\n"
              << "absolue L_infty error: " << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() << "\n"
              << "relative L_infty error: "
              << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() / jacobians[ii][0].infinity_norm();
      }
    }
  } // ... basis_jacobians_seem_to_be_correct(...)
}; // struct ContinuousLagrangeSpace


template <class G, int p>
struct ContinuousLagrangeSpaceOnSimplicialLeafView
    : public ContinuousLagrangeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  GridProviderType grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  ContinuousLagrangeSpaceOnSimplicialLeafView() //        (i) negative coordinates and not the same as the reference
      : grid_provider(Dune::XT::Grid::make_cube_grid<G>(-1.5, -1, 3).grid_ptr()) //                          element,
  { //                                                   (ii) at least 3 elements to have fully inner ones,
    grid_provider.global_refine(1); //                  (iii) refine at least once to obtain all kinds of orientations
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider.leaf_view());
  }

  ~ContinuousLagrangeSpaceOnSimplicialLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct ContinuousLagrangeSpaceOnSimplicialLeafView


using SimplicialGrids = ::testing::Types<ONED_1D,
                                         YASP_1D_EQUIDISTANT_OFFSET
//#if HAVE_DUNE_ALUGRID
//                                         ,
//                                         ALU_2D_SIMPLEX_CONFORMING,
//                                         ALU_2D_SIMPLEX_NONCONFORMING
//#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_2D
#endif
//#if HAVE_DUNE_ALUGRID
//                                         ,
//                                         ALU_3D_SIMPLEX_CONFORMING,
//                                         ALU_3D_SIMPLEX_NONCONFORMING
//#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_3D
#endif
                                         >;


template <class G>
using Order1SimplicialContinuousLagrangeSpace = ContinuousLagrangeSpaceOnSimplicialLeafView<G, 1>;
TYPED_TEST_CASE(Order1SimplicialContinuousLagrangeSpace, SimplicialGrids);
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1SimplicialContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2SimplicialContinuousLagrangeSpace = ContinuousLagrangeSpaceOnSimplicialLeafView<G, 2>;
TYPED_TEST_CASE(Order2SimplicialContinuousLagrangeSpace, SimplicialGrids);
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2SimplicialContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct ContinuousLagrangeSpaceOnCubicLeafView
    : public ContinuousLagrangeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  ContinuousLagrangeSpaceOnCubicLeafView()
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

  ~ContinuousLagrangeSpaceOnCubicLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct ContinuousLagrangeSpaceOnCubicLeafView


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
using Order1CubicContinuousLagrangeSpace = ContinuousLagrangeSpaceOnCubicLeafView<G, 1>;
TYPED_TEST_CASE(Order1CubicContinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order1CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1CubicContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2CubicContinuousLagrangeSpace = ContinuousLagrangeSpaceOnCubicLeafView<G, 2>;
TYPED_TEST_CASE(Order2CubicContinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct ContinuousLagrangeSpaceOnPrismLeafView
    : public ContinuousLagrangeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  ContinuousLagrangeSpaceOnPrismLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    if (d == 3) {
      Dune::GridFactory<G> factory;
      for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                            Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.})}) {
        factory.insertVertex(vertex);
      }
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism, 3), {0, 1, 2, 3, 4, 5});
      grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
      grid_provider->global_refine(1);
      leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
    } else {
      // cannot use ASSERT_... in a ctor
      EXPECT_TRUE(false) << "Does not make sense in " << d << "d!\n"
                         << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
      grid_provider = nullptr;
      leaf_view = nullptr;
    }
  } // ContinuousLagrangeSpaceOnPrismLeafView(...)

  ~ContinuousLagrangeSpaceOnPrismLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct ContinuousLagrangeSpaceOnPrismLeafView


#if HAVE_DUNE_UGGRID || HAVE_UG


using PrismGrids = ::testing::Types<UG_3D>;

template <class G>
using Order1PrismContinuousLagrangeSpace = ContinuousLagrangeSpaceOnPrismLeafView<G, 1>;
TYPED_TEST_CASE(Order1PrismContinuousLagrangeSpace, PrismGrids);
TYPED_TEST(Order1PrismContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1PrismContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2PrismContinuousLagrangeSpace = ContinuousLagrangeSpaceOnPrismLeafView<G, 2>;
TYPED_TEST_CASE(Order2PrismContinuousLagrangeSpace, PrismGrids);
TYPED_TEST(Order2PrismContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2PrismContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct ContinuousLagrangeSpaceOnMixedLeafView
    : public ContinuousLagrangeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  ContinuousLagrangeSpaceOnMixedLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
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
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {3, 0, 4, 1});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {4, 1, 5, 2});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 6, 3});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 5, 6});
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
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {3, 0, 4, 1, 9, 6, 10, 7});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {4, 1, 5, 2, 10, 7, 11, 8});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 12, 3, 10});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 5, 12, 10});
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
  } // ContinuousLagrangeSpaceOnMixedLeafView(...)

  ~ContinuousLagrangeSpaceOnMixedLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct ContinuousLagrangeSpaceOnMixedLeafView


// The mapper does not work in 3d.
using MixedGrids = ::testing::Types<UG_2D>;

template <class G>
using Order1MixedContinuousLagrangeSpace = ContinuousLagrangeSpaceOnMixedLeafView<G, 1>;
TYPED_TEST_CASE(Order1MixedContinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order1MixedContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1MixedContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2MixedContinuousLagrangeSpace = ContinuousLagrangeSpaceOnMixedLeafView<G, 2>;
TYPED_TEST_CASE(Order2MixedContinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order2MixedContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2MixedContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


#endif // HAVE_DUNE_UGGRID || HAVE_UG
