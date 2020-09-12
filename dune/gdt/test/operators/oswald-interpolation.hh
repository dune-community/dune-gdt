// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_OPERATORS_OSWALD_INTERPOLATION_HH
#define DUNE_GDT_TEST_OPERATORS_OSWALD_INTERPOLATION_HH

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/test/common.hh>

#include <dune/xt/la/container/istl.hh>

#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/boundaryinfo/normalbased.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/structuredgridfactory.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/localizable-functional.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/abs.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/localizable-operator.hh>
#include <dune/gdt/operators/oswald-interpolation.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct OswaldInterpolationOperatorOnLeafViewTest : public ::testing::Test
{
  static_assert(XT::Grid::is_grid<G>::value, "");

  using GV = typename G::LeafGridView;
  using D = typename GV::ctype;
  static constexpr size_t d = GV::dimension;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using M = XT::LA::IstlRowMajorSparseMatrix<double>;
  using V = XT::LA::IstlDenseVector<double>;

  std::shared_ptr<XT::Grid::GridProvider<G>> grid_provider;
  std::shared_ptr<DiscreteFunction<V, GV>> source;
  std::shared_ptr<XT::Grid::NormalBasedBoundaryInfo<I>> boundary_info;
  std::shared_ptr<DiscontinuousLagrangeSpace<GV>> space;

  virtual std::shared_ptr<XT::Grid::GridProvider<G>> make_grid()
  {
    return std::make_shared<XT::Grid::GridProvider<G>>(
        XT::Grid::make_cube_grid<G>(XT::Common::from_string<FieldVector<double, d>>("[0 0 0 0]"),
                                    XT::Common::from_string<FieldVector<double, d>>("[3 1 1 1]"),
                                    XT::Common::from_string<std::array<unsigned int, d>>("[3 1 1 1]")));
  }

  void SetUp() override
  {
    grid_provider = make_grid();
    space = std::make_shared<DiscontinuousLagrangeSpace<GV>>(
        make_discontinuous_lagrange_space(grid_provider->leaf_view(), 1));
    source = std::make_shared<DiscreteFunction<V, GV>>(*space);
    auto local_source = source->local_discrete_function();
    for (auto&& element : elements(space->grid_view())) {
      local_source->bind(element);
      local_source->dofs().assign_from(space->finite_elements()
                                           .get(element.type(), 1)
                                           .interpolation()
                                           .interpolate(
                                               [&](const auto& x_local) {
                                                 const auto center = element.geometry().center();
                                                 const auto x = element.geometry().global(x_local);
                                                 if (center[0] < 1.)
                                                   return x[0] + 1; // going from 1 to 2
                                                 else if (center[0] < 2.)
                                                   return x[0] + 2; // going from 3 to 4
                                                 else
                                                   return x[0] + 3; // going from 5 to 6
                                               },
                                               1));
    }
    boundary_info = std::make_shared<XT::Grid::NormalBasedBoundaryInfo<I>>();
    boundary_info->register_new_normal(XT::Common::from_string<FieldVector<double, d>>("[-1 0 0 0]"),
                                       new XT::Grid::DirichletBoundary());
  } // ... SetUp(...)

  void fulfills_l2_interpolation_estimate()
  {
    auto& self = *this;
    OswaldInterpolationOperator<M, GV> oswald_interpolation(
        self.grid_provider->leaf_view(), *self.space, *self.space, *self.boundary_info);
    oswald_interpolation.assemble();
    auto range = make_discrete_function<V>(*self.space);
    oswald_interpolation.apply(*self.source, range);
    const auto difference = *self.source - range;
    auto fv_space = make_finite_volume_space(self.space->grid_view());
    // assemble elementwise L2 norms
    auto l2_element_indicators = make_discrete_function(fv_space);
    auto walker = XT::Grid::make_walker(self.space->grid_view());
    walker.append([]() {},
                  [&](const auto& element) {
                    auto local_difference = difference.local_function();
                    local_difference->bind(element);
                    auto local_indicator = l2_element_indicators.local_discrete_function(element);
                    local_indicator->dofs()[0] += LocalElementIntegralBilinearForm<E>(LocalProductIntegrand<E>())
                                                      .apply2(*local_difference, *local_difference)[0][0];
                  },
                  []() {});
    // assemble intersection L2 jump norms, attach to vertices
    auto cg_space = make_continuous_lagrange_space(self.space->grid_view(), 1);
    auto vertex_indicators = make_discrete_function(cg_space);
    walker.append(
        [](/*prepare nothing*/) {},
        [&](const auto& intersection, const auto& inside_element, const auto& outside_element) {
          // localize data functions
          DUNE_THROW_IF(!intersection.conforming(), Exceptions::operator_error, "Not implemented yet!");
          auto difference_inside = difference.local_function();
          difference_inside->bind(inside_element);
          auto difference_outside = difference.local_function();
          difference_outside->bind(outside_element);
          const auto intersection_diameter =
              (d > 1) ? XT::Grid::diameter(intersection)
                      : (std::min(XT::Grid::diameter(inside_element), XT::Grid::diameter(outside_element)));
          // compute L2 jump norm
          const auto l2_jump_norm2 =
              LocalIntersectionIntegralFunctional<I>(
                  [](const auto& inside_function, const auto& outside_function, const auto& param) {
                    return std::max(inside_function.order(param), outside_function.order(param));
                  },
                  [&](const auto& inside_basis,
                      const auto& outside_basis,
                      const auto& point_in_reference_intersection,
                      auto& result,
                      const auto& param) {
                    const auto point_in_inside_reference_element =
                        intersection.geometryInInside().global(point_in_reference_intersection);
                    const auto point_in_outside_reference_element =
                        intersection.geometryInOutside().global(point_in_reference_intersection);
                    const auto inside_values = inside_basis.evaluate_set(point_in_inside_reference_element, param);
                    const auto outside_values = outside_basis.evaluate_set(point_in_outside_reference_element, param);
                    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
                      for (size_t jj = 0; jj < outside_basis.size(param); ++jj)
                        result[ii][jj] = std::pow(inside_values[ii] - outside_values[ii], 2);
                  })
                  .apply(intersection, *difference_inside, *difference_outside)[0][0];
          // associate intersection corners with global vertex IDs
          auto local_vertex_indicators_inside = vertex_indicators.local_discrete_function(inside_element);
          const auto& inside_reference_element = ReferenceElements<D, d>::general(inside_element.type());
          for (auto&& ii : XT::Common::value_range(inside_reference_element.size(intersection.indexInInside(), 1, d))) {
            const auto vertex_index_inside = inside_reference_element.subEntity(intersection.indexInInside(), 1, ii, d);
            local_vertex_indicators_inside->dofs()[vertex_index_inside] += intersection_diameter * l2_jump_norm2;
          }
          auto local_vertex_indicators_outside = vertex_indicators.local_discrete_function(outside_element);
          const auto& outside_reference_element = ReferenceElements<D, d>::general(outside_element.type());
          for (auto&& ii :
               XT::Common::value_range(outside_reference_element.size(intersection.indexInOutside(), 1, d))) {
            const auto vertex_index_outside =
                outside_reference_element.subEntity(intersection.indexInOutside(), 1, ii, d);
            local_vertex_indicators_outside->dofs()[vertex_index_outside] += intersection_diameter * l2_jump_norm2;
          }
        },
        [](/*finalize nothing*/) {},
        XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
    walker.append(
        [](/*prepare nothing*/) {},
        [&](const auto& intersection, const auto& inside_element, const auto& /*outside_element*/) {
          // localize data functions
          DUNE_THROW_IF(!intersection.conforming(), Exceptions::operator_error, "Not implemented yet!");
          auto difference_inside = difference.local_function();
          difference_inside->bind(inside_element);
          const auto intersection_diameter =
              (d > 1) ? XT::Grid::diameter(intersection) : XT::Grid::diameter(inside_element);
          // compute L2 jump norm
          const auto l2_jump_norm2 =
              LocalIntersectionIntegralFunctional<I>(
                  [](const auto& inside_function, const auto& /*outside_function*/, const auto& param) {
                    return inside_function.order(param);
                  },
                  [&](const auto& inside_basis,
                      const auto& /*outside_basis*/,
                      const auto& point_in_reference_intersection,
                      auto& result,
                      const auto& param) {
                    const auto point_in_inside_reference_element =
                        intersection.geometryInInside().global(point_in_reference_intersection);
                    const auto inside_values = inside_basis.evaluate_set(point_in_inside_reference_element, param);
                    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
                      result[ii][0] = std::pow(inside_values[ii], 2);
                  })
                  .apply(intersection, *difference_inside, *difference_inside)[0][0];
          // associate intersection corners with global vertex IDs
          auto local_vertex_indicators_inside = vertex_indicators.local_discrete_function(inside_element);
          const auto& inside_reference_element = ReferenceElements<D, d>::general(inside_element.type());
          for (auto&& ii : XT::Common::value_range(inside_reference_element.size(intersection.indexInInside(), 1, d))) {
            const auto vertex_index_inside = inside_reference_element.subEntity(intersection.indexInInside(), 1, ii, d);
            local_vertex_indicators_inside->dofs()[vertex_index_inside] += intersection_diameter * l2_jump_norm2;
          }
        },
        [](/*finalize nothing*/) {},
        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(*self.boundary_info, new XT::Grid::DirichletBoundary()));
    walker.walk();
    // collect jump norms to an element indicator
    auto l2_jump_indicators = make_discrete_function(fv_space);
    walker.append([](/*prepare nothing*/) {},
                  [&](const auto& element) {
                    auto local_indicator = l2_jump_indicators.local_discrete_function(element);
                    for (auto&& global_index : cg_space.mapper().global_indices(element))
                      local_indicator->dofs()[0] += vertex_indicators.dofs().vector()[global_index];
                  },
                  [](/*finalize nothing*/) {});
    walker.walk();
    // now that we have both, we can compute the first constant from ESV2007, Lemma 3.5
    double C = std::numeric_limits<double>::min();
    for (size_t ii = 0; ii < l2_element_indicators.dofs().vector().size(); ++ii)
      C = std::max(C, l2_element_indicators.dofs().vector()[ii] / l2_jump_indicators.dofs().vector()[ii]);
    const auto expected_C = DXTC_TEST_CONFIG_GET("results.constant_from_L2_interpolation_estimate", -1.);
    EXPECT_DOUBLE_EQ(expected_C, C) << "XT::Common::Test::get_unique_test_name() = '"
                                    << XT::Common::Test::get_unique_test_name() << "'";
  } // ... fulfills_l2_interpolation_estimate(...)

  void fulfills_h1_interpolation_estimate()
  {
    auto& self = *this;
    OswaldInterpolationOperator<M, GV> oswald_interpolation(
        self.grid_provider->leaf_view(), *self.space, *self.space, *self.boundary_info);
    oswald_interpolation.assemble();
    auto range = make_discrete_function<V>(*self.space);
    oswald_interpolation.apply(*self.source, range);
    const auto difference = *self.source - range;
    auto fv_space = make_finite_volume_space(self.space->grid_view());
    // assemble elementwise H1 semi norms
    auto h1_element_indicators = make_discrete_function(fv_space);
    auto walker = XT::Grid::make_walker(self.space->grid_view());
    walker.append([]() {},
                  [&](const auto& element) {
                    auto local_difference = difference.local_function();
                    local_difference->bind(element);
                    auto local_indicator = h1_element_indicators.local_discrete_function(element);
                    local_indicator->dofs()[0] += LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(1.))
                                                      .apply2(*local_difference, *local_difference)[0][0];
                  },
                  []() {});
    // assemble intersection L2 jump norms, attach to vertices (code is nearly identycal to
    // fulfills_l2_interpolation_estimate(), except for the power of the intersection_diameter)
    auto cg_space = make_continuous_lagrange_space(self.space->grid_view(), 1);
    auto vertex_indicators = make_discrete_function(cg_space);
    walker.append(
        [](/*prepare nothing*/) {},
        [&](const auto& intersection, const auto& inside_element, const auto& outside_element) {
          // localize data functions
          DUNE_THROW_IF(!intersection.conforming(), Exceptions::operator_error, "Not implemented yet!");
          auto difference_inside = difference.local_function();
          difference_inside->bind(inside_element);
          auto difference_outside = difference.local_function();
          difference_outside->bind(outside_element);
          const auto intersection_diameter =
              (d > 1) ? XT::Grid::diameter(intersection)
                      : (std::min(XT::Grid::diameter(inside_element), XT::Grid::diameter(outside_element)));
          // compute L2 jump norm
          const auto l2_jump_norm2 =
              LocalIntersectionIntegralFunctional<I>(
                  [](const auto& inside_function, const auto& outside_function, const auto& param) {
                    return std::max(inside_function.order(param), outside_function.order(param));
                  },
                  [&](const auto& inside_basis,
                      const auto& outside_basis,
                      const auto& point_in_reference_intersection,
                      auto& result,
                      const auto& param) {
                    const auto point_in_inside_reference_element =
                        intersection.geometryInInside().global(point_in_reference_intersection);
                    const auto point_in_outside_reference_element =
                        intersection.geometryInOutside().global(point_in_reference_intersection);
                    const auto inside_values = inside_basis.evaluate_set(point_in_inside_reference_element, param);
                    const auto outside_values = outside_basis.evaluate_set(point_in_outside_reference_element, param);
                    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
                      for (size_t jj = 0; jj < outside_basis.size(param); ++jj)
                        result[ii][jj] = std::pow(inside_values[ii] - outside_values[ii], 2);
                  })
                  .apply(intersection, *difference_inside, *difference_outside)[0][0];
          // associate intersection corners with global vertex IDs
          auto local_vertex_indicators_inside = vertex_indicators.local_discrete_function(inside_element);
          const auto& inside_reference_element = ReferenceElements<D, d>::general(inside_element.type());
          for (auto&& ii : XT::Common::value_range(inside_reference_element.size(intersection.indexInInside(), 1, d))) {
            const auto vertex_index_inside = inside_reference_element.subEntity(intersection.indexInInside(), 1, ii, d);
            local_vertex_indicators_inside->dofs()[vertex_index_inside] += l2_jump_norm2 / intersection_diameter;
          }
          auto local_vertex_indicators_outside = vertex_indicators.local_discrete_function(outside_element);
          const auto& outside_reference_element = ReferenceElements<D, d>::general(outside_element.type());
          for (auto&& ii :
               XT::Common::value_range(outside_reference_element.size(intersection.indexInOutside(), 1, d))) {
            const auto vertex_index_outside =
                outside_reference_element.subEntity(intersection.indexInOutside(), 1, ii, d);
            local_vertex_indicators_outside->dofs()[vertex_index_outside] += l2_jump_norm2 / intersection_diameter;
          }
        },
        [](/*finalize nothing*/) {},
        XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
    walker.append(
        [](/*prepare nothing*/) {},
        [&](const auto& intersection, const auto& inside_element, const auto& /*outside_element*/) {
          // localize data functions
          DUNE_THROW_IF(!intersection.conforming(), Exceptions::operator_error, "Not implemented yet!");
          auto difference_inside = difference.local_function();
          difference_inside->bind(inside_element);
          const auto intersection_diameter =
              (d > 1) ? XT::Grid::diameter(intersection) : XT::Grid::diameter(inside_element);
          // compute L2 jump norm
          const auto l2_jump_norm2 =
              LocalIntersectionIntegralFunctional<I>(
                  [](const auto& inside_function, const auto& /*outside_function*/, const auto& param) {
                    return inside_function.order(param);
                  },
                  [&](const auto& inside_basis,
                      const auto& /*outside_basis*/,
                      const auto& point_in_reference_intersection,
                      auto& result,
                      const auto& param) {
                    const auto point_in_inside_reference_element =
                        intersection.geometryInInside().global(point_in_reference_intersection);
                    const auto inside_values = inside_basis.evaluate_set(point_in_inside_reference_element, param);
                    for (size_t ii = 0; ii < inside_basis.size(param); ++ii)
                      result[ii][0] = std::pow(inside_values[ii], 2);
                  })
                  .apply(intersection, *difference_inside, *difference_inside)[0][0];
          // associate intersection corners with global vertex IDs
          auto local_vertex_indicators_inside = vertex_indicators.local_discrete_function(inside_element);
          const auto& inside_reference_element = ReferenceElements<D, d>::general(inside_element.type());
          for (auto&& ii : XT::Common::value_range(inside_reference_element.size(intersection.indexInInside(), 1, d))) {
            const auto vertex_index_inside = inside_reference_element.subEntity(intersection.indexInInside(), 1, ii, d);
            local_vertex_indicators_inside->dofs()[vertex_index_inside] += l2_jump_norm2 / intersection_diameter;
          }
        },
        [](/*finalize nothing*/) {},
        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(*self.boundary_info, new XT::Grid::DirichletBoundary()));
    walker.walk();
    // collect jump norms to an element indicator
    auto l2_jump_indicators = make_discrete_function(fv_space);
    walker.append([](/*prepare nothing*/) {},
                  [&](const auto& element) {
                    auto local_indicator = l2_jump_indicators.local_discrete_function(element);
                    for (auto&& global_index : cg_space.mapper().global_indices(element))
                      local_indicator->dofs()[0] += vertex_indicators.dofs().vector()[global_index];
                  },
                  [](/*finalize nothing*/) {});
    walker.walk();
    // now that we have both, we can compute the second constant from ESV2007, Lemma 3.5
    double C = std::numeric_limits<double>::min();
    for (size_t ii = 0; ii < h1_element_indicators.dofs().vector().size(); ++ii)
      C = std::max(C, h1_element_indicators.dofs().vector()[ii] / l2_jump_indicators.dofs().vector()[ii]);
    const auto expected_C = DXTC_TEST_CONFIG_GET("results.constant_from_H1_interpolation_estimate", -1.);
    EXPECT_DOUBLE_EQ(expected_C, C) << "XT::Common::Test::get_unique_test_name() = '"
                                    << XT::Common::Test::get_unique_test_name() << "'";
  } // ... fulfills_h1_interpolation_estimate(...)
}; // struct OswaldInterpolationOperatorOnLeafViewTest


template <class G>
struct OswaldInterpolationOperatorOnCubicLeafViewTest : public OswaldInterpolationOperatorOnLeafViewTest<G>
{
  using BaseType = OswaldInterpolationOperatorOnLeafViewTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::M;
  using typename BaseType::V;

  std::shared_ptr<DiscreteFunction<V, GV>> expected_range_for_cubic_grids;

  std::shared_ptr<XT::Grid::GridProvider<G>> make_grid() override final
  {
    FieldVector<D, d> lower_left(0.);
    auto upper_right = XT::Common::from_string<FieldVector<double, d>>("[3 1 1 1]");
    std::array<unsigned int, d> num_elements;
    std::fill(num_elements.begin(), num_elements.end(), 1);
    num_elements[0] = 3;
    return std::make_shared<XT::Grid::GridProvider<G>>(
        XT::Grid::StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
  }

  void SetUp() override final
  {
    BaseType::SetUp();
    auto& self = *this;
    expected_range_for_cubic_grids = std::make_shared<DiscreteFunction<V, GV>>(*self.space);
    auto local_expected_range_for_cubic_grids = expected_range_for_cubic_grids->local_discrete_function();
    for (auto&& element : elements(self.space->grid_view())) {
      local_expected_range_for_cubic_grids->bind(element);
      local_expected_range_for_cubic_grids->dofs().assign_from(
          self.space->finite_elements()
              .get(element.type(), 1)
              .interpolation()
              .interpolate(
                  [&](const auto& x_local) {
                    const auto center = element.geometry().center();
                    const auto x = element.geometry().global(x_local);
                    if (center[0] < 1.)
                      return 2.5 * x[0]; // going from 0 to 2.5
                    else if (center[0] < 2.)
                      return 2.0 * x[0] + 0.5; // going from 2.5 to 4.5
                    else
                      return 1.5 * x[0] + 1.5; // going from 4.5 to 6
                  },
                  1));
    }
  } // ... SetUp(...)

  void applies_correctly_on_cubic_grids()
  {
    auto& self = *this;
    OswaldInterpolationOperator<M, GV> oswald_interpolation(
        self.grid_provider->leaf_view(), *self.space, *self.space, *self.boundary_info);
    oswald_interpolation.assemble();
    auto range = make_discrete_function<V>(*self.space);
    oswald_interpolation.apply(*self.source, range);
    const auto difference = *self.expected_range_for_cubic_grids - range;
    auto L_infinity_error = make_localizable_functional(self.space->grid_view(), difference);
    L_infinity_error.append(LocalElementIntegralFunctional<E>(LocalElementAbsIntegrand<E>()));
    L_infinity_error.assemble();
    EXPECT_LT(L_infinity_error.result(), 1e-15);
  } // ... applies_correctly_on_cubic_grids(...)
}; // struct OswaldInterpolationOperatorOnCubicLeafViewTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_OSWALD_INTERPOLATION_HH
