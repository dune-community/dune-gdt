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

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/discretefunction/default.hh>


template <class GridLayerType, int p>
struct DiscontinuousLagrangeSpace : public ::testing::Test
{
  using SpaceType = Dune::GDT::DiscontinuousLagrangeSpace<GridLayerType, p>;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;

  virtual std::shared_ptr<GridLayerType> grid_layer() = 0;

  std::shared_ptr<SpaceType> space;

  ~DiscontinuousLagrangeSpace() = default;

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
    // we want to check that the numbering is consecutive and that each global index exists only once
    std::set<size_t> global_indices;
    // we test both call variants
    std::set<size_t> map_to_global;
    for (auto&& element : Dune::elements(*grid_layer())) {
      for (const auto& global_index : space->mapper().globalIndices(element))
        global_indices.insert(global_index);
      for (size_t ii = 0; ii < space->mapper().numDofs(element); ++ii)
        map_to_global.insert(space->mapper().mapToGlobal(element, ii));
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
    double tolerance = 1e-15;
    for (auto&& element : elements(*grid_layer())) {
      const auto basis = space->base_function_set(element);
      const auto lagrange_points = space->lagrange_points(element);
      EXPECT_EQ(lagrange_points.size(), basis.size());
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto values = basis.evaluate(lagrange_points[ii]);
        for (size_t jj = 0; jj < basis.size(); ++jj) {
          ASSERT_TRUE(Dune::XT::Common::FloatCmp::eq(values[jj][0], ii == jj ? 1. : 0., tolerance, tolerance))
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
        double tolerance = 1e-4;
        for (size_t ii = 0; ii < basis.size(); ++ii)
          EXPECT_TRUE(
              Dune::XT::Common::FloatCmp::eq(jacobians[ii][0], approximate_jacobians[ii][0], tolerance, tolerance))
              << "ii = " << ii << "\njacobians[ii][0] = " << jacobians[ii][0] << "\n"
              << "approximate_jacobians[ii][0] = " << approximate_jacobians[ii][0] << "\n"
              << "absolue L_infty error: " << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() << "\n"
              << "relative L_infty error: "
              << (jacobians[ii][0] - approximate_jacobians[ii][0]).infinity_norm() / jacobians[ii][0].infinity_norm();
      }
    }
  } // ... basis_jacobians_seem_to_be_correct(...)
}; // struct DiscontinuousLagrangeSpace


template <class G, int p>
struct DiscontinuousLagrangeSpaceOnSimplicialLeafView
    : public DiscontinuousLagrangeSpace<
          typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType,
          p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  GridProviderType grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  DiscontinuousLagrangeSpaceOnSimplicialLeafView() //     (i) negative coordinates and not the same as the reference
      : grid_provider(Dune::XT::Grid::make_cube_grid<G>(-1.5, -1, 3).grid_ptr()) //                          element,
  { //                                                   (ii) at least 3 elements to have fully inner ones,
    grid_provider.global_refine(1); //                  (iii) refine at least once to obtain all kinds of orientations
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider.leaf_view());
  }

  ~DiscontinuousLagrangeSpaceOnSimplicialLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct DiscontinuousLagrangeSpaceOnSimplicialLeafView


using SimplicialGrids = ::testing::Types<YASP_1D_EQUIDISTANT_OFFSET
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
using Order0SimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafView<G, 0>;
TYPED_TEST_CASE(Order0SimplicialDiscontinuousLagrangeSpace, SimplicialGrids);
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0SimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order1SimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafView<G, 1>;
TYPED_TEST_CASE(Order1SimplicialDiscontinuousLagrangeSpace, SimplicialGrids);
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1SimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2SimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafView<G, 2>;
TYPED_TEST_CASE(Order2SimplicialDiscontinuousLagrangeSpace, SimplicialGrids);
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2SimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct DiscontinuousLagrangeSpaceOnCubicLeafView
    : public DiscontinuousLagrangeSpace<
          typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType,
          p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  DiscontinuousLagrangeSpaceOnCubicLeafView()
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

  ~DiscontinuousLagrangeSpaceOnCubicLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct DiscontinuousLagrangeSpaceOnCubicLeafView


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
using Order0CubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafView<G, 0>;
TYPED_TEST_CASE(Order0CubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0CubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order1CubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafView<G, 1>;
TYPED_TEST_CASE(Order1CubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1CubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2CubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafView<G, 2>;
TYPED_TEST_CASE(Order2CubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2CubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct DiscontinuousLagrangeSpaceOnPrismLeafView
    : public DiscontinuousLagrangeSpace<
          typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType,
          p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  DiscontinuousLagrangeSpaceOnPrismLeafView()
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
  } // DiscontinuousLagrangeSpaceOnPrismLeafView(...)

  ~DiscontinuousLagrangeSpaceOnPrismLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct DiscontinuousLagrangeSpaceOnPrismLeafView


using PrismGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D
#endif
    >;


template <class G>
using Order0PrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafView<G, 0>;
TYPED_TEST_CASE(Order0PrismDiscontinuousLagrangeSpace, PrismGrids);
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0PrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order1PrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafView<G, 1>;
TYPED_TEST_CASE(Order1PrismDiscontinuousLagrangeSpace, PrismGrids);
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1PrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2PrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafView<G, 2>;
TYPED_TEST_CASE(Order2PrismDiscontinuousLagrangeSpace, PrismGrids);
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2PrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, int p>
struct DiscontinuousLagrangeSpaceOnMixedLeafView
    : public DiscontinuousLagrangeSpace<
          typename Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>::LeafGridViewType,
          p>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  DiscontinuousLagrangeSpaceOnMixedLeafView()
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
  } // DiscontinuousLagrangeSpaceOnMixedLeafView(...)

  ~DiscontinuousLagrangeSpaceOnMixedLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_layer() override final
  {
    return leaf_view;
  }
}; // struct DiscontinuousLagrangeSpaceOnMixedLeafView


using MixedGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D
//  , UG_3D
#endif
    >;


template <class G>
using Order0MixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafView<G, 0>;
TYPED_TEST_CASE(Order0MixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0MixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order1MixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafView<G, 1>;
TYPED_TEST_CASE(Order1MixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1MixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G>
using Order2MixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafView<G, 2>;
TYPED_TEST_CASE(Order2MixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2MixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}

GTEST_TEST(bindingslike, main)
{
  using namespace Dune;
  using namespace Dune::GDT;

  using Grid = GDT_BINDINGS_GRID;
  using Spc = DgSpaceProvider<Grid, XT::Grid::Layers::leaf, Backends::gdt, 1, double, 1>;
  using Vector = XT::LA::IstlDenseVector<>;
  auto gp = XT::Grid::make_cube_grid<Grid>();
  auto space = Spc::create(gp);
  auto df = make_discrete_function<Vector>(space, "STATE");
  df.visualize("foo");
}