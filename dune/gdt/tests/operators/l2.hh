// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_L2_HH

#include <dune/stuff/common/string.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/l2.hh>

#include "weighted-l2.hh"

namespace Dune {
namespace GDT {
namespace Tests {


template <class SpaceType>
struct L2LocalizableProductTest : public WeightedL2ProductBase<SpaceType>, public LocalizableProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  using typename LocalizableBaseType::GridViewType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using typename LocalizableBaseType::ScalarFunctionType;
  using typename LocalizableBaseType::RangeFieldType;

  L2LocalizableProductTest()
    : WeightedL2BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    typedef L2LocalizableProduct<GridViewType, ScalarFunctionType, ScalarFunctionType, double> CtorTestProductType;
    CtorTestProductType DUNE_UNUSED(wo_over_integrate)(grid_view, range, source);
    CtorTestProductType DUNE_UNUSED(with_over_integrate)(1, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    auto DUNE_UNUSED(wo_over_integrate) = make_l2_localizable_product(grid_view, range, source);
    auto DUNE_UNUSED(with_over_integrate) = make_l2_localizable_product(grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_view = this->space_.grid_view();

    auto product      = make_l2_localizable_product(grid_view, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_l2_localizable_product(grid_view, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    EXPECT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_localizable_product()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    auto product = make_l2_localizable_product(grid_view, range, source);
    this->localizable_product_test(*product);
  } // ... is_localizable_product(...)
}; // struct L2LocalizableProductTest


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_L2_HH
