// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_INTERPOLATIONS_BOUNDARY_HH
#define DUNE_GDT_TEST_INTERPOLATIONS_BOUNDARY_HH

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/test/common.hh>
#include <dune/xt/grid/boundaryinfo/normalbased.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/structuredgridfactory.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/norms.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct BoundaryInterpolationOnLeafViewTest : public ::testing::Test
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
  std::shared_ptr<XT::Functions::ConstantFunction<d>> source;
  std::shared_ptr<XT::Grid::NormalBasedBoundaryInfo<I>> boundary_info;
  std::shared_ptr<ContinuousLagrangeSpace<GV>> space;
  std::shared_ptr<DiscreteFunction<V, GV>> range;
  std::shared_ptr<XT::Functions::GenericFunction<d>> expected_range;

  virtual std::shared_ptr<XT::Grid::GridProvider<G>> make_grid()
  {
    return std::make_shared<XT::Grid::GridProvider<G>>(
        XT::Grid::make_cube_grid<G>(XT::Common::from_string<FieldVector<double, d>>("[-1 0 0 0]"),
                                    XT::Common::from_string<FieldVector<double, d>>("[1 1 1 1]"),
                                    XT::Common::from_string<std::array<unsigned int, d>>("[2 1 1 1]")));
  }

  void SetUp() override
  {
    grid_provider = make_grid();
    space =
        std::make_shared<ContinuousLagrangeSpace<GV>>(make_continuous_lagrange_space(grid_provider->leaf_view(), 1));
    source = std::make_shared<XT::Functions::ConstantFunction<d>>(1.);
    range = std::make_shared<DiscreteFunction<V, GV>>(*space);
    expected_range = std::make_shared<XT::Functions::GenericFunction<d>>([](const auto&) { return 1; },
                                                                         [](const auto& x, const auto&) {
                                                                           if (x[0] > 0)
                                                                             return x[0];
                                                                           else
                                                                             return 0.;
                                                                         });
    boundary_info = std::make_shared<XT::Grid::NormalBasedBoundaryInfo<I>>();
    boundary_info->register_new_normal(XT::Common::from_string<FieldVector<double, d>>("[1 0 0 0]"),
                                       new XT::Grid::DirichletBoundary());
  } // ... SetUp(...)

  void interpolates_correctly()
  {
    boundary_interpolation(XT::Functions::make_grid_function<E>(*source),
                           *range,
                           space->grid_view(),
                           *boundary_info,
                           XT::Grid::DirichletBoundary());
    const auto expected_L2_error = 1e-15;
    const auto l2_error = l2_norm(space->grid_view(), XT::Functions::make_grid_function<E>(*expected_range) - *range);
    EXPECT_LT(l2_error, expected_L2_error)
        << "XT::Common::Test::get_unique_test_name() = '" << XT::Common::Test::get_unique_test_name() << "'";

  } // ... interpolates_correctly(...)
}; // struct BoundaryInterpolationOnLeafViewTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INTERPOLATIONS_BOUNDARY_HH
