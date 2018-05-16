// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/test/float_cmp.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/elliptic.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {

struct EllipticDefaultTolerances
{
#ifndef NDEBUG
  static constexpr const double constant = 5.05e-13;
  static constexpr const double linear = 1.75e-13;
  static constexpr const double quadratic = 9.83e-13;
#else
  static constexpr const double constant = 6.54e-13;
  static constexpr const double linear = 2.65e-13;
  static constexpr const double quadratic = 9.83e-13;
#endif // #ifndef NDEBUG
};
/**
 * \note The values in correct_for_constant_arguments(), etc., are valid for the d-dimendional unit cube.
 */
template <class SpaceType>
struct EllipticProductBase
{
  typedef typename SpaceType::GridLayerType GridLayerType;
  typedef XT::Grid::extract_grid_t<GridLayerType> GridType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1>
      ExpressionFunctionType;

  EllipticProductBase(const double factor_value = 42.)
    : factor_value_(factor_value)
    , factor_("x", Dune::XT::Common::to_string(factor_value_), 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
    , constant_gradient_("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
    , linear_gradient_("x", "fake_value", 2, "affine gradient", {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}})
    , quadratic_gradient_("x", "fake_value", 3, "quadratic gradient", {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}})
  {
  }

  virtual ~EllipticProductBase() = default;

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const = 0;

  void correct_for_constant_arguments(const RangeFieldType epsilon = 1e-14) const
  {
    check(compute(constant_gradient_), factor_value_ * dimDomain * 1.0, epsilon);
  }

  void correct_for_linear_arguments(const RangeFieldType epsilon = 1e-15) const
  {
    check(compute(linear_gradient_), factor_value_ * dimDomain * (1.0 / 3.0), epsilon);
  }

  void correct_for_quadratic_arguments(const RangeFieldType epsilon = 1e-15) const
  {
    check(compute(quadratic_gradient_), factor_value_ * dimDomain * (1.0 / 5.0), epsilon);
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon) const
  {
    //! might be off, since atol was used before
    DXTC_EXPECT_FLOAT_LE(expected, result, epsilon);
  } // ... check(...)

  const double factor_value_;
  const ExpressionFunctionType factor_;
  const ExpressionFunctionType constant_gradient_;
  const ExpressionFunctionType linear_gradient_;
  const ExpressionFunctionType quadratic_gradient_;
}; // struct EllipticProductBase


template <class SpaceType>
struct EllipticLocalizableProductTest : public EllipticProductBase<SpaceType>, public LocalizableProductBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> EllipticBaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  typedef typename LocalizableBaseType::GridLayerType GridLayerType;
  typedef typename EllipticBaseType::ExpressionFunctionType ExpressionFunctionType;
  typedef typename LocalizableBaseType::ScalarFunctionType ScalarFunctionType;
  typedef typename LocalizableBaseType::TensorFunctionType TensorFunctionType;
  typedef typename LocalizableBaseType::RangeFieldType RangeFieldType;
  using EllipticBaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef EllipticLocalizableProduct<ExpressionFunctionType,
                                       void,
                                       GridLayerType,
                                       ScalarFunctionType,
                                       ScalarFunctionType,
                                       double>
        OnlyFactorType;
    DUNE_UNUSED OnlyFactorType factor_only(diffusion_factor, grid_layer, range, source);
    DUNE_UNUSED OnlyFactorType factor_only_with_over_integrate(1, diffusion_factor, grid_layer, range, source);

    typedef EllipticLocalizableProduct<TensorFunctionType,
                                       void,
                                       GridLayerType,
                                       ExpressionFunctionType,
                                       ExpressionFunctionType,
                                       double>
        OnlyTensorType;
    DUNE_UNUSED OnlyTensorType tensor_only(diffusion_tensor, grid_layer, range, source);
    DUNE_UNUSED OnlyTensorType tensor_only_with_over_integrate(1, diffusion_tensor, grid_layer, range, source);

    typedef EllipticLocalizableProduct<ExpressionFunctionType,
                                       TensorFunctionType,
                                       GridLayerType,
                                       ExpressionFunctionType,
                                       ExpressionFunctionType,
                                       double>
        BothType;
    DUNE_UNUSED BothType both(diffusion_factor, diffusion_tensor, grid_layer, range, source);
    DUNE_UNUSED BothType both_with_over_integrate(1, diffusion_factor, diffusion_tensor, grid_layer, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto factor_only DUNE_UNUSED = make_elliptic_localizable_product(diffusion_factor, grid_layer, range, source);
    auto factor_only_with_over_integrate DUNE_UNUSED =
        make_elliptic_localizable_product(diffusion_factor, grid_layer, range, source, 1);

    auto tensor_only DUNE_UNUSED = make_elliptic_localizable_product(diffusion_tensor, grid_layer, range, source);
    auto tensor_only_with_over_integrate DUNE_UNUSED =
        make_elliptic_localizable_product(diffusion_tensor, grid_layer, range, source, 1);

    auto both DUNE_UNUSED =
        make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_layer, range, source);
    auto both_with_over_integrate DUNE_UNUSED =
        make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_layer, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    auto product = make_elliptic_localizable_product(
        this->factor_, this->tensor_function_, this->space_.grid_layer(), function, function);
    const auto result = product->apply2();
    auto product_tbb = make_elliptic_localizable_product(
        this->factor_, this->tensor_function_, this->space_.grid_layer(), function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();
    DXTC_EXPECT_FLOAT_EQ(result_tbb, result, 1e-14);
    return result;
  }

  void is_localizable_product()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto product = make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_layer, range, source);
    this->localizable_product_test(*product);
  }
}; // struct EllipticLocalizableProductTest

/**
 * \note Assumes that Operators::Projection does the right thing!
 */
template <class SpaceType>
struct EllipticMatrixOperatorTest : public EllipticProductBase<SpaceType>, public MatrixOperatorBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> EllipticBaseType;
  typedef MatrixOperatorBase<SpaceType> MatrixBaseType;
  typedef typename MatrixBaseType::GridLayerType GridLayerType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename MatrixBaseType::DiscreteFunctionType;
  using typename MatrixBaseType::ScalarFunctionType;
  using typename MatrixBaseType::TensorFunctionType;
  typedef typename MatrixBaseType::RangeFieldType RangeFieldType;
  using typename MatrixBaseType::MatrixType;

  EllipticMatrixOperatorTest(const double factor_value = 42.)
    : EllipticBaseType(factor_value)
  {
  }

  void constructible_by_ctor()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();

    // only diffusion factor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    using FourArgsOp = EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType>;
    DUNE_UNUSED FourArgsOp factor_no_matrix_1(diffusion_factor, space);
    DUNE_UNUSED FourArgsOp factor_no_matrix_2(diffusion_factor, space, grid_layer);
    DUNE_UNUSED FourArgsOp factor_no_matrix_3(diffusion_factor, space, space, grid_layer);
    //       full argument list
    using EightArgsOp =
        EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridLayerType, SpaceType, double>;
    EightArgsOp DUNE_UNUSED factor_no_matrix_4(diffusion_factor, space);
    EightArgsOp DUNE_UNUSED factor_no_matrix_5(diffusion_factor, space, grid_layer);
    EightArgsOp DUNE_UNUSED factor_no_matrix_6(diffusion_factor, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED FourArgsOp factor_no_matrix_7(1, diffusion_factor, space);
    DUNE_UNUSED FourArgsOp factor_no_matrix_8(1, diffusion_factor, space, grid_layer);
    DUNE_UNUSED FourArgsOp factor_no_matrix_9(1, diffusion_factor, space, space, grid_layer);
    //       full argument list
    EightArgsOp DUNE_UNUSED factor_no_matrix_10(1, diffusion_factor, space);
    EightArgsOp DUNE_UNUSED factor_no_matrix_11(1, diffusion_factor, space, grid_layer);
    EightArgsOp DUNE_UNUSED factor_no_matrix_12(1, diffusion_factor, space, space, grid_layer);
    //   with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //     without over_integrate
    //       simplified argument list
    DUNE_UNUSED FourArgsOp factor_matrix_1(diffusion_factor, matrix, space);
    DUNE_UNUSED FourArgsOp factor_matrix_2(diffusion_factor, matrix, space, grid_layer);
    DUNE_UNUSED FourArgsOp factor_matrix_3(diffusion_factor, matrix, space, space, grid_layer);
    //       full argument list
    EightArgsOp DUNE_UNUSED factor_matrix_4(diffusion_factor, matrix, space);
    EightArgsOp DUNE_UNUSED factor_matrix_5(diffusion_factor, matrix, space, grid_layer);
    EightArgsOp DUNE_UNUSED factor_matrix_6(diffusion_factor, matrix, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED FourArgsOp factor_matrix_7(1, diffusion_factor, matrix, space);
    DUNE_UNUSED FourArgsOp factor_matrix_8(1, diffusion_factor, matrix, space, grid_layer);
    DUNE_UNUSED FourArgsOp factor_matrix_9(1, diffusion_factor, matrix, space, space, grid_layer);
    //       full argument list
    EightArgsOp DUNE_UNUSED factor_matrix_10(1, diffusion_factor, matrix, space);
    EightArgsOp DUNE_UNUSED factor_matrix_11(1, diffusion_factor, matrix, space, grid_layer);
    EightArgsOp DUNE_UNUSED factor_matrix_12(1, diffusion_factor, matrix, space, space, grid_layer);

    // only diffusion tensor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    using TensorFourArgsOp = EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType>;
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_1(diffusion_tensor, space);
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_2(diffusion_tensor, space, grid_layer);
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_3(diffusion_tensor, space, space, grid_layer);
    //       full argument list
    using TensorSevenArgsOp =
        EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridLayerType, SpaceType, double>;
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_4(diffusion_tensor, space);
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_5(diffusion_tensor, space, grid_layer);
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_6(diffusion_tensor, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_7(1, diffusion_tensor, space);
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_8(1, diffusion_tensor, space, grid_layer);
    DUNE_UNUSED TensorFourArgsOp tensor_no_matrix_9(1, diffusion_tensor, space, space, grid_layer);
    //       full argument list
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_10(1, diffusion_tensor, space);
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_11(1, diffusion_tensor, space, grid_layer);
    TensorSevenArgsOp DUNE_UNUSED tensor_no_matrix_12(1, diffusion_tensor, space, space, grid_layer);
    //   with matrix
    //     without over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_1(diffusion_tensor, matrix, space);
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_2(diffusion_tensor, matrix, space, grid_layer);
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_3(diffusion_tensor, matrix, space, space, grid_layer);
    //       full argument list
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_4(diffusion_tensor, matrix, space);
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_5(diffusion_tensor, matrix, space, grid_layer);
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_6(diffusion_tensor, matrix, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_7(1, diffusion_tensor, matrix, space);
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_8(1, diffusion_tensor, matrix, space, grid_layer);
    DUNE_UNUSED TensorFourArgsOp tensor_matrix_9(1, diffusion_tensor, matrix, space, space, grid_layer);
    //       full argument list
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_10(1, diffusion_tensor, matrix, space);
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_11(1, diffusion_tensor, matrix, space, grid_layer);
    TensorSevenArgsOp DUNE_UNUSED tensor_matrix_12(1, diffusion_tensor, matrix, space, space, grid_layer);

    // both diffusion factor and tensor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    using TensorB_FourArgsOp =
        EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType, MatrixType>;
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_1(diffusion_factor, diffusion_tensor, space);
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_2(diffusion_factor, diffusion_tensor, space, grid_layer);
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_3(diffusion_factor, diffusion_tensor, space, space, grid_layer);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_4(diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_5(diffusion_factor, diffusion_tensor, space, grid_layer);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_6(diffusion_factor, diffusion_tensor, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_7(1, diffusion_factor, diffusion_tensor, space);
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_8(1, diffusion_factor, diffusion_tensor, space, grid_layer);
    DUNE_UNUSED TensorB_FourArgsOp both_no_matrix_9(1, diffusion_factor, diffusion_tensor, space, space, grid_layer);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_10(1, diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_11(1, diffusion_factor, diffusion_tensor, space, grid_layer);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_no_matrix_12(1, diffusion_factor, diffusion_tensor, space, space, grid_layer);
    //   with matrix
    //     without over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_1(diffusion_factor, diffusion_tensor, matrix, space);
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_2(diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_3(diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_4(diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_5(diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_6(diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer);
    //     with over_integrate
    //       simplified argument list
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_7(1, diffusion_factor, diffusion_tensor, matrix, space);
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_8(1, diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
    DUNE_UNUSED TensorB_FourArgsOp both_matrix_9(
        1, diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_10(1, diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_11(1, diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridLayerType,
                           SpaceType,
                           double>
        DUNE_UNUSED both_matrix_12(1, diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // both diffusion factor and tensor
    //   without matrix
    auto op_01 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    auto op_02 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, 1);
    auto op_03 DUNE_UNUSED =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, grid_layer);
    auto op_04 DUNE_UNUSED =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, grid_layer, 1);
    auto op_05 DUNE_UNUSED =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, space, grid_layer);
    auto op_06 DUNE_UNUSED =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, space, grid_layer, 1);
    //   with matrix
    auto op_07 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space);
    auto op_08 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, 1);
    auto op_09 DUNE_UNUSED =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
    auto op_10 DUNE_UNUSED =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, grid_layer, 1);
    auto op_11 DUNE_UNUSED =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer);
    auto op_12 DUNE_UNUSED =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, space, grid_layer, 1);

    // single diffusion factor
    //   without matrix
    auto op_13 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space);
    auto op_14 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, 1);
    auto op_15 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, grid_layer);
    auto op_16 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, grid_layer, 1);
    auto op_17 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, space, grid_layer);
    auto op_18 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, space, grid_layer, 1);
    //   with matrix
    auto op_19 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space);
    auto op_20 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, 1);
    auto op_21 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_layer);
    auto op_22 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_layer, 1);
    auto op_23 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_layer);
    auto op_24 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_layer, 1);

    // single diffusion tensor
    //   without matrix
    auto op_25 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space);
    auto op_26 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, 1);
    auto op_27 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, grid_layer);
    auto op_28 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, grid_layer, 1);
    auto op_29 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, space, grid_layer);
    auto op_30 DUNE_UNUSED = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, space, grid_layer, 1);
    //   with matrix
    auto op_31 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space);
    auto op_32 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, 1);
    auto op_33 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_layer);
    auto op_34 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_layer, 1);
    auto op_35 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_layer);
    auto op_36 DUNE_UNUSED = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_layer, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space = this->space_;
    auto op = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(space.grid_layer(), function, discrete_function, 2);
    // compute product
    const auto result = op->apply2(discrete_function, discrete_function);
    auto op_tbb = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    op_tbb->assemble(true);
    const auto result_tbb = op_tbb->apply2(discrete_function, discrete_function);
    DXTC_EXPECT_FLOAT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_constant_arguments(const double tolerance = EllipticDefaultTolerances::constant) const
  {
    const ExpressionFunctionType constant_gradient("x", "x[0]", 1, "constant gradient", {{"1.0", "0.0", "0.0"}});
    this->check(compute(constant_gradient), factor_value_ * 1.0, tolerance);
  }

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_linear_arguments(const double tolerance = EllipticDefaultTolerances::linear) const
  {
    const ExpressionFunctionType linear_gradient(
        "x", "0.5 * x[0] * x[0] - x[0]", 2, "affine gradient", {{"x[0] - 1.0", "0.0", "0.0"}});

    this->check(compute(linear_gradient), factor_value_ * 1.0 / 3.0, tolerance);
  }

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_quadratic_arguments(const double tolerance = EllipticDefaultTolerances::constant) const
  {
    const ExpressionFunctionType quadratic_gradient(
        "x", "(1.0/3.0) * x[0] * x[0] * x[0]", 3, ", quadratic gradient", {{"x[0]*x[0]", "0.0", "0.0"}});
    this->check(compute(quadratic_gradient), factor_value_ * 1.0 / 5.0, tolerance);
  }

  void is_matrix_operator()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space = this->space_;

    auto op = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    this->matrix_operator_test(*op);
  } // ... is_matrix_operator(...)

  using EllipticBaseType::factor_value_;
}; // struct EllipticMatrixOperatorTest


template <class SpaceType>
struct EllipticOperatorTest : public EllipticProductBase<SpaceType>, public OperatorBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> EllipticBaseType;
  typedef OperatorBase<SpaceType> OperatorBaseType;
  using typename OperatorBaseType::GridLayerType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename OperatorBaseType::DiscreteFunctionType;
  using typename OperatorBaseType::ScalarFunctionType;
  using typename OperatorBaseType::TensorFunctionType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::MatrixType;
  using typename OperatorBaseType::VectorType;
  using EllipticBaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();

    DUNE_UNUSED EllipticOperator<ExpressionFunctionType, void, GridLayerType> only_factor(grid_layer, diffusion_factor);
    DUNE_UNUSED EllipticOperator<ExpressionFunctionType, void, GridLayerType> only_factor_w_over(
        1, grid_layer, diffusion_factor);

    DUNE_UNUSED EllipticOperator<TensorFunctionType, void, GridLayerType> only_tensor(grid_layer, diffusion_tensor);
    DUNE_UNUSED EllipticOperator<TensorFunctionType, void, GridLayerType> only_tensor_w_over(
        1, grid_layer, diffusion_tensor);

    DUNE_UNUSED EllipticOperator<ExpressionFunctionType, TensorFunctionType, GridLayerType> both(
        grid_layer, diffusion_factor, diffusion_tensor);
    DUNE_UNUSED EllipticOperator<ExpressionFunctionType, TensorFunctionType, GridLayerType> both_w_over(
        1, grid_layer, diffusion_factor, diffusion_tensor);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();

    auto only_factor DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_factor);
    auto only_factor_w_over DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_factor, 1);

    auto only_tensor DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_tensor);
    auto only_tensor_w_over DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_tensor, 1);

    auto both DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_factor, diffusion_tensor);
    auto both_w_over DUNE_UNUSED = make_elliptic_operator(grid_layer, diffusion_factor, diffusion_tensor, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();

    return make_elliptic_operator(grid_layer, diffusion_factor, diffusion_tensor)->apply2(function, function);
  } // ... compute(...)

  void apply_is_callable()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_layer = this->space_.grid_layer();
    auto& source = this->discrete_function_;
    auto range = make_discrete_function<VectorType>(this->space_);

    auto op = make_elliptic_operator(grid_layer, diffusion_factor, diffusion_tensor);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct EllipticOperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
