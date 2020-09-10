// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_INTERPOLATIONS_DEFAULT_HH
#define DUNE_GDT_TEST_INTERPOLATIONS_DEFAULT_HH

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/fmatrix.hh>
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
struct DefaultInterpolationOnLeafViewTest : public ::testing::Test
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
  std::shared_ptr<XT::Functions::GenericFunction<d>> source;
  std::shared_ptr<ContinuousLagrangeSpace<GV>> space;
  std::shared_ptr<DiscreteFunction<V, GV>> range;

  virtual std::shared_ptr<XT::Grid::GridProvider<G>> make_grid()
  {
    return std::make_shared<XT::Grid::GridProvider<G>>(
        XT::Grid::make_cube_grid<G>(XT::Common::from_string<FieldVector<double, d>>("[-1 0 0 0]"),
                                    XT::Common::from_string<FieldVector<double, d>>("[1 1 1 1]"),
                                    XT::Common::from_string<std::array<unsigned int, d>>("[5 5 5 2]")));
  }

  void SetUp() override
  {
    grid_provider = make_grid();
    space =
        std::make_shared<ContinuousLagrangeSpace<GV>>(make_continuous_lagrange_space(grid_provider->leaf_view(), 2));
    source = std::make_shared<XT::Functions::GenericFunction<d>>(
        [](const auto&) { return 2; },
        [](const auto& x, const auto&) {
          const auto x_dependent = 2 * std::pow(x[0], 2) - x[0] + 3;
          D xy_dependent;
          if constexpr (d >= 2)
            xy_dependent = x_dependent + x[0] * x[1] + 0.5 * x[1] - std::pow(x[1], 2);
          if constexpr (d == 1)
            return x_dependent;
          else if constexpr (d == 2)
            return xy_dependent;
          else
            return xy_dependent + 0.5 * std::pow(x[2], 2) + x[2] * x[0] - 3 * x[2] * x[1];
        },
        "second order polynomial",
        XT::Common::ParameterType{},
        [](const auto& x, const auto&) {
          const std::vector<double> x_dependent_jacobian{4 * x[0] - 1, 0, 0, 0};
          std::vector<double> y_dependent_jacobian, z_dependent_jacobian;
          if constexpr (d >= 2)
            y_dependent_jacobian = {x[1], x[0] + 0.5 - 2 * x[1], 0, 0};
          if constexpr (d >= 3)
            z_dependent_jacobian = {x[2], -3 * x[2], x[2] + x[0] - 3 * x[1], 0};
          XT::Common::FieldMatrix<double, 1, d> jacobian;
          for (size_t ii = 0; ii < d; ++ii) {
            jacobian[0][ii] = x_dependent_jacobian[ii];
            if constexpr (d >= 2)
              jacobian[0][ii] += y_dependent_jacobian[ii];
            if constexpr (d >= 3)
              jacobian[0][ii] += z_dependent_jacobian[ii];
          }
          return jacobian;
        });
    range = std::make_shared<DiscreteFunction<V, GV>>(*space);
  } // ... SetUp(...)

  void interpolates_correctly(const double expected_l2_error = 1e-14)
  {
    default_interpolation(*source, *range, space->grid_view());
    const auto l2_error = l2_norm(space->grid_view(), source->template as_grid_function<E>() - *range);
    EXPECT_LT(l2_error, expected_l2_error)
        << "XT::Common::Test::get_unique_test_name() = '" << XT::Common::Test::get_unique_test_name() << "'";
    const auto local_range = range->local_discrete_function();
    for (auto&& element : Dune::elements(space->grid_view())) {
      local_range->bind(element);
      const auto center = element.geometry().center();
      const auto range_jacobian = local_range->jacobian(element.geometry().local(center));
      const auto expected_jacobian = source->jacobian(center);
      EXPECT_LT((range_jacobian - expected_jacobian)[0].two_norm(), expected_l2_error);
    }
  } // ... interpolates_correctly(...)
}; // struct DefaultInterpolationOnLeafViewTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INTERPOLATIONS_DEFAULT_HH
