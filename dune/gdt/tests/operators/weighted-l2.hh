// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH

#include <dune/stuff/common/string.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/operators/weighted-l2.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \note The values in correct_for_constant_arguments(), etc., are valid for the d-dimendional unit cube.
 */
template <class SpaceType>
struct WeightedL2ProductBase
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1>
      ExpressionFunctionType;

  static const int weight_value = 42;

  WeightedL2ProductBase() // linker error if int(...) is missing, at least with clang
      : weight_("x", DSC::toString(int(weight_value)), 0),
        constant_("x", "1.0", 0),
        linear_("x", "x[0] - 1.0", 1),
        quadratic_("x", "x[0]*x[0]", 2)
  {
  }

  virtual ~WeightedL2ProductBase() = default;

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const = 0;

  void correct_for_constant_arguments() const
  {
    check(compute(constant_), weight_value * 1.0);
  }

  void correct_for_linear_arguments() const
  {
    check(compute(linear_), weight_value * (1.0 / 3.0));
  }

  void correct_for_quadratic_arguments() const
  {
    check(compute(quadratic_), weight_value * (1.0 / 5.0));
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected,
             const RangeFieldType epsilon = 2.14e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon) << "result:     " << result << "\n"
                              << "expected:   " << expected << "\n"
                              << "difference: " << std::scientific << error;
  } // ... check(...)

  const ExpressionFunctionType weight_;
  const ExpressionFunctionType constant_;
  const ExpressionFunctionType linear_;
  const ExpressionFunctionType quadratic_;
}; // struct WeightedL2ProductBase


template <class SpaceType>
struct WeightedL2LocalizableProductTest : public WeightedL2ProductBase<SpaceType>,
                                          public LocalizableProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> EllipticBaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  typedef typename LocalizableBaseType::GridViewType GridViewType;
  typedef typename EllipticBaseType::ExpressionFunctionType ExpressionFunctionType;
  typedef typename LocalizableBaseType::ScalarFunctionType ScalarFunctionType;
  typedef typename LocalizableBaseType::TensorFunctionType TensorFunctionType;
  typedef typename LocalizableBaseType::RangeFieldType RangeFieldType;

  void constructible_by_ctor()
  {
    const auto& weight    = this->weight_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    typedef WeightedL2LocalizableProduct<ExpressionFunctionType,
                                         GridViewType,
                                         ScalarFunctionType,
                                         ScalarFunctionType,
                                         double> CtorTestProductType;
    CtorTestProductType DUNE_UNUSED(wo_over_integrate)(weight, grid_view, range, source);
    CtorTestProductType DUNE_UNUSED(with_over_integrate)(1, weight, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& weight    = this->weight_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    auto DUNE_UNUSED(wo_over_integrate) = make_weighted_l2_localizable_product(weight, grid_view, range, source);
    auto DUNE_UNUSED(with_over_integrate) = make_weighted_l2_localizable_product(weight, grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& weight    = this->weight_;
    const auto& grid_view = this->space_.grid_view();

    auto product      = make_weighted_l2_localizable_product(weight, grid_view, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_weighted_l2_localizable_product(weight, grid_view, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    EXPECT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_localizable_product()
  {
    const auto& weight    = this->weight_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source    = this->scalar_function_;
    const auto& range     = this->scalar_function_;

    auto product = make_weighted_l2_localizable_product(weight, grid_view, range, source);
    this->localizable_product_test(*product);
  } // ... is_localizable_product(...)
}; // struct WeightedL2LocalizableProductTest


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH
