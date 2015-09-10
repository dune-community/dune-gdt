// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH

#include <dune/gdt/operators/elliptic.hh>

#include "../operators.hh"

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \note The values in correct_for_constant_arguments(), etc., are valid for the d-dimendional unit cube.
 */
template< class SpaceType >
struct EllipticProductBase
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid      GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  typedef Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > ExpressionFunctionType;

  EllipticProductBase()
   : one_("x", "1.0", 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
   , constant_gradient_("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
   , linear_gradient_("x", "fake_value", 2, "affine gradient", {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}})
   , quadratic_gradient_("x", "fake_value", 3, ", quadratic gradient", {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}})
  {}

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const = 0;

  void correct_for_constant_arguments() const
  {
    check(compute(constant_gradient_), dimDomain*1.0);
  }

  void correct_for_linear_arguments() const
  {
    check(compute(linear_gradient_), dimDomain*(1.0/3.0));
  }

  void correct_for_quadratic_arguments() const
  {
    check(compute(quadratic_gradient_), dimDomain*(1.0/5.0));
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon = 1e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon)
        << "result:     " << result << "\n"
        << "expected:   " << expected << "\n"
        << "difference: " << std::scientific << error;
  } // ... check(...)

  const ExpressionFunctionType one_;
  const ExpressionFunctionType constant_gradient_;
  const ExpressionFunctionType linear_gradient_;
  const ExpressionFunctionType quadratic_gradient_;
}; // struct EllipticProductBase


template< class SpaceType >
struct EllipticLocalizableProductTest
  : public EllipticProductBase< SpaceType >
  , public LocalizableProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > EllipticBaseType;
  typedef LocalizableProductBase< SpaceType > LocalizableBaseType;
  typedef typename LocalizableBaseType::GridViewType   GridViewType;
  typedef typename EllipticBaseType::ExpressionFunctionType ExpressionFunctionType;
  typedef typename LocalizableBaseType::ScalarFunctionType     ScalarFunctionType;
  typedef typename LocalizableBaseType::TensorFunctionType     TensorFunctionType;
  typedef typename LocalizableBaseType::RangeFieldType RangeFieldType;

  void constructible_by_ctor()
  {
    const auto& diffusion_factor = this->one_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef EllipticLocalizableProduct
        < ExpressionFunctionType, void, GridViewType, ScalarFunctionType, ScalarFunctionType, double >
      OnlyFactorType;
    OnlyFactorType DUNE_UNUSED(factor_only)(                       diffusion_factor, grid_view, range, source);
    OnlyFactorType DUNE_UNUSED(factor_only_with_over_integrate)(1, diffusion_factor, grid_view, range, source);

    typedef EllipticLocalizableProduct
        < TensorFunctionType, void, GridViewType, ExpressionFunctionType, ExpressionFunctionType, double >
      OnlyTensorType;
    OnlyTensorType DUNE_UNUSED(tensor_only)(                       diffusion_tensor, grid_view, range, source);
    OnlyTensorType DUNE_UNUSED(tensor_only_with_over_integrate)(1, diffusion_tensor, grid_view, range, source);

    typedef EllipticLocalizableProduct
        < ExpressionFunctionType, TensorFunctionType, GridViewType, ExpressionFunctionType, ExpressionFunctionType, double >
      BothType;
    BothType DUNE_UNUSED(both)(                       diffusion_factor, diffusion_tensor, grid_view, range, source);
    BothType DUNE_UNUSED(both_with_over_integrate)(1, diffusion_factor, diffusion_tensor, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->one_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto DUNE_UNUSED(factor_only) = make_elliptic_localizable_product(diffusion_factor, grid_view, range, source);
    auto DUNE_UNUSED(factor_only_with_over_integrate)
                                  = make_elliptic_localizable_product(diffusion_factor, grid_view, range, source, 1);

    auto DUNE_UNUSED(tensor_only) = make_elliptic_localizable_product(diffusion_tensor, grid_view, range, source);
    auto DUNE_UNUSED(tensor_only_with_over_integrate)
                                  = make_elliptic_localizable_product(diffusion_tensor, grid_view, range, source, 1);

    auto DUNE_UNUSED(both)
        = make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source);
    auto DUNE_UNUSED(both_with_over_integrate)
        = make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    return make_elliptic_localizable_product(this->one_, this->tensor_function_, this->space_.grid_view(), function,
                                             function)->apply2();
  }

  void is_localizable_product()
  {
    const auto& diffusion_factor = this->one_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto product = make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source);
    this->localizable_product_test(*product);
  }
}; // struct EllipticLocalizableProductTest


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
