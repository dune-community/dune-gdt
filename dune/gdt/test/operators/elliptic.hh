// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/elliptic.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \note The values in correct_for_constant_arguments(), etc., are valid for the d-dimendional unit cube.
 */
template <class SpaceType>
struct EllipticProductBase
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

  EllipticProductBase(const double factor_value = 42.)
    : factor_value_(factor_value)
    , factor_("x", DSC::to_string(factor_value_), 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
    , constant_gradient_("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
    , linear_gradient_("x", "fake_value", 2, "affine gradient", {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}})
    , quadratic_gradient_("x", "fake_value", 3, "quadratic gradient", {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}})
  {
  }

  virtual ~EllipticProductBase() = default;

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const = 0;

  void correct_for_constant_arguments(const RangeFieldType epsilon = 1e-15) const
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
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon) << "result:     " << result << "\n"
                              << "expected:   " << expected << "\n"
                              << "difference: " << std::scientific << error;
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
  typedef typename LocalizableBaseType::GridViewType GridViewType;
  typedef typename EllipticBaseType::ExpressionFunctionType ExpressionFunctionType;
  typedef typename LocalizableBaseType::ScalarFunctionType ScalarFunctionType;
  typedef typename LocalizableBaseType::TensorFunctionType TensorFunctionType;
  typedef typename LocalizableBaseType::RangeFieldType RangeFieldType;
  using EllipticBaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();
    const auto& source           = this->scalar_function_;
    const auto& range            = this->scalar_function_;

    typedef EllipticLocalizableProduct<ExpressionFunctionType,
                                       void,
                                       GridViewType,
                                       ScalarFunctionType,
                                       ScalarFunctionType,
                                       double> OnlyFactorType;
    OnlyFactorType DUNE_UNUSED(factor_only)(diffusion_factor, grid_view, range, source);
    OnlyFactorType DUNE_UNUSED(factor_only_with_over_integrate)(1, diffusion_factor, grid_view, range, source);

    typedef EllipticLocalizableProduct<TensorFunctionType,
                                       void,
                                       GridViewType,
                                       ExpressionFunctionType,
                                       ExpressionFunctionType,
                                       double> OnlyTensorType;
    OnlyTensorType DUNE_UNUSED(tensor_only)(diffusion_tensor, grid_view, range, source);
    OnlyTensorType DUNE_UNUSED(tensor_only_with_over_integrate)(1, diffusion_tensor, grid_view, range, source);

    typedef EllipticLocalizableProduct<ExpressionFunctionType,
                                       TensorFunctionType,
                                       GridViewType,
                                       ExpressionFunctionType,
                                       ExpressionFunctionType,
                                       double> BothType;
    BothType DUNE_UNUSED(both)(diffusion_factor, diffusion_tensor, grid_view, range, source);
    BothType DUNE_UNUSED(both_with_over_integrate)(1, diffusion_factor, diffusion_tensor, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();
    const auto& source           = this->scalar_function_;
    const auto& range            = this->scalar_function_;

    auto DUNE_UNUSED(factor_only) = make_elliptic_localizable_product(diffusion_factor, grid_view, range, source);
    auto DUNE_UNUSED(factor_only_with_over_integrate) =
        make_elliptic_localizable_product(diffusion_factor, grid_view, range, source, 1);

    auto DUNE_UNUSED(tensor_only) = make_elliptic_localizable_product(diffusion_tensor, grid_view, range, source);
    auto DUNE_UNUSED(tensor_only_with_over_integrate) =
        make_elliptic_localizable_product(diffusion_tensor, grid_view, range, source, 1);

    auto DUNE_UNUSED(both) =
        make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source);
    auto DUNE_UNUSED(both_with_over_integrate) =
        make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    auto product = make_elliptic_localizable_product(
        this->factor_, this->tensor_function_, this->space_.grid_view(), function, function);
    const auto result = product->apply2();
    auto product_tbb = make_elliptic_localizable_product(
        this->factor_, this->tensor_function_, this->space_.grid_view(), function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();
    EXPECT_EQ(result_tbb, result);
    return result;
  }

  void is_localizable_product()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();
    const auto& source           = this->scalar_function_;
    const auto& range            = this->scalar_function_;

    auto product = make_elliptic_localizable_product(diffusion_factor, diffusion_tensor, grid_view, range, source);
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
  typedef typename MatrixBaseType::GridViewType GridViewType;
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
    const auto& space            = this->space_;
    const auto& grid_view        = this->space_.grid_view();

    // only diffusion factor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_1)(diffusion_factor,
                                                                                                    space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_2)(
        diffusion_factor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_3)(
        diffusion_factor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_4)(diffusion_factor, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_5)(diffusion_factor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_6)(diffusion_factor, space, space, grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_7)(
        1, diffusion_factor, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_8)(
        1, diffusion_factor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_no_matrix_9)(
        1, diffusion_factor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_10)(1, diffusion_factor, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_11)(1, diffusion_factor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_no_matrix_12)(1, diffusion_factor, space, space, grid_view);
    //   with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_1)(
        diffusion_factor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_2)(
        diffusion_factor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_3)(
        diffusion_factor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_4)(diffusion_factor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_5)(diffusion_factor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_6)(diffusion_factor, matrix, space, space, grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_7)(
        1, diffusion_factor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_8)(
        1, diffusion_factor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType> DUNE_UNUSED(factor_matrix_9)(
        1, diffusion_factor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_10)(1, diffusion_factor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_11)(1, diffusion_factor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(factor_matrix_12)(1, diffusion_factor, matrix, space, space, grid_view);

    // only diffusion tensor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_1)(diffusion_tensor,
                                                                                                space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_2)(
        diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_3)(
        diffusion_tensor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_4)(diffusion_tensor, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_5)(diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_6)(diffusion_tensor, space, space, grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_7)(
        1, diffusion_tensor, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_8)(
        1, diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_no_matrix_9)(
        1, diffusion_tensor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_10)(1, diffusion_tensor, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_11)(1, diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_no_matrix_12)(1, diffusion_tensor, space, space, grid_view);
    //   with matrix
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_1)(
        diffusion_tensor, matrix, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_2)(
        diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_3)(
        diffusion_tensor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_4)(diffusion_tensor, matrix, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_5)(diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_6)(diffusion_tensor, matrix, space, space, grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_7)(
        1, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_8)(
        1, diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType> DUNE_UNUSED(tensor_matrix_9)(
        1, diffusion_tensor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_10)(1, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_11)(1, diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<TensorFunctionType, void, SpaceType, MatrixType, GridViewType, SpaceType, double>
        DUNE_UNUSED(tensor_matrix_12)(1, diffusion_tensor, matrix, space, space, grid_view);

    // both diffusion factor and tensor
    //   without matrix
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_1)(
        diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_2)(
        diffusion_factor, diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_3)(
        diffusion_factor, diffusion_tensor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_4)(diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_5)(diffusion_factor, diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_6)(diffusion_factor,
                                                                 diffusion_tensor,
                                                                 space,
                                                                 space,
                                                                 grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_7)(
        1, diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_8)(
        1, diffusion_factor, diffusion_tensor, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_no_matrix_9)(
        1, diffusion_factor, diffusion_tensor, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_10)(1, diffusion_factor, diffusion_tensor, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_11)(1,
                                                                  diffusion_factor,
                                                                  diffusion_tensor,
                                                                  space,
                                                                  grid_view);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_no_matrix_12)(1,
                                                                  diffusion_factor,
                                                                  diffusion_tensor,
                                                                  space,
                                                                  space,
                                                                  grid_view);
    //   with matrix
    //     without over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_1)(
        diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_2)(
        diffusion_factor, diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_3)(
        diffusion_factor, diffusion_tensor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_matrix_4)(diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_matrix_5)(diffusion_factor,
                                                              diffusion_tensor,
                                                              matrix,
                                                              space,
                                                              grid_view);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_matrix_6)(diffusion_factor,
                                                              diffusion_tensor,
                                                              matrix,
                                                              space,
                                                              space,
                                                              grid_view);
    //     with over_integrate
    //       simplified argument list
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_7)(
        1, diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_8)(
        1, diffusion_factor, diffusion_tensor, matrix, space, grid_view);
    EllipticMatrixOperator<ExpressionFunctionType, TensorFunctionType, SpaceType> DUNE_UNUSED(both_matrix_9)(
        1, diffusion_factor, diffusion_tensor, matrix, space, space, grid_view);
    //       full argument list
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_matrix_10)(1, diffusion_factor, diffusion_tensor, matrix, space);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double> DUNE_UNUSED(both_matrix_11)(1,
                                                               diffusion_factor,
                                                               diffusion_tensor,
                                                               matrix,
                                                               space,
                                                               grid_view);
    EllipticMatrixOperator<ExpressionFunctionType,
                           TensorFunctionType,
                           SpaceType,
                           MatrixType,
                           GridViewType,
                           SpaceType,
                           double>
        DUNE_UNUSED(both_matrix_12)(1, diffusion_factor, diffusion_tensor, matrix, space, space, grid_view);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space            = this->space_;
    const auto& grid_view = this->space_.grid_view();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // both diffusion factor and tensor
    //   without matrix
    auto DUNE_UNUSED(op_01) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    auto DUNE_UNUSED(op_02) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, 1);
    auto DUNE_UNUSED(op_03) =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, grid_view);
    auto DUNE_UNUSED(op_04) =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, grid_view, 1);
    auto DUNE_UNUSED(op_05) =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, space, grid_view);
    auto DUNE_UNUSED(op_06) =
        make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space, space, grid_view, 1);
    //   with matrix
    auto DUNE_UNUSED(op_07) = make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space);
    auto DUNE_UNUSED(op_08) = make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, 1);
    auto DUNE_UNUSED(op_09) =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, grid_view);
    auto DUNE_UNUSED(op_10) =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, grid_view, 1);
    auto DUNE_UNUSED(op_11) =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, space, grid_view);
    auto DUNE_UNUSED(op_12) =
        make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, space, grid_view, 1);

    // single diffusion factor
    //   without matrix
    auto DUNE_UNUSED(op_13) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space);
    auto DUNE_UNUSED(op_14) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, 1);
    auto DUNE_UNUSED(op_15) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, grid_view);
    auto DUNE_UNUSED(op_16) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, grid_view, 1);
    auto DUNE_UNUSED(op_17) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, space, grid_view);
    auto DUNE_UNUSED(op_18) = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, space, space, grid_view, 1);
    //   with matrix
    auto DUNE_UNUSED(op_19) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space);
    auto DUNE_UNUSED(op_20) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, 1);
    auto DUNE_UNUSED(op_21) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_view);
    auto DUNE_UNUSED(op_22) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_view, 1);
    auto DUNE_UNUSED(op_23) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_view);
    auto DUNE_UNUSED(op_24) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_view, 1);

    // single diffusion tensor
    //   without matrix
    auto DUNE_UNUSED(op_25) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space);
    auto DUNE_UNUSED(op_26) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, 1);
    auto DUNE_UNUSED(op_27) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, grid_view);
    auto DUNE_UNUSED(op_28) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, grid_view, 1);
    auto DUNE_UNUSED(op_29) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, space, grid_view);
    auto DUNE_UNUSED(op_30) = make_elliptic_matrix_operator<MatrixType>(diffusion_tensor, space, space, grid_view, 1);
    //   with matrix
    auto DUNE_UNUSED(op_31) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space);
    auto DUNE_UNUSED(op_32) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, 1);
    auto DUNE_UNUSED(op_33) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_view);
    auto DUNE_UNUSED(op_34) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, grid_view, 1);
    auto DUNE_UNUSED(op_35) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_view);
    auto DUNE_UNUSED(op_36) = make_elliptic_matrix_operator(diffusion_tensor, matrix, space, space, grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space            = this->space_;
    auto op                      = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(space.grid_view(), function, discrete_function, 2);
    // compute product
    const auto result = op->apply2(discrete_function, discrete_function);
    auto op_tbb = make_elliptic_matrix_operator<MatrixType>(diffusion_factor, diffusion_tensor, space);
    op_tbb->assemble(true);
    const auto result_tbb = op_tbb->apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_constant_arguments() const
  {
    const ExpressionFunctionType constant_gradient("x", "x[0]", 1, "constant gradient", {{"1.0", "0.0", "0.0"}});
    this->check(compute(constant_gradient), factor_value_ * 1.0, 5.05e-13);
  }

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_linear_arguments() const
  {
    const ExpressionFunctionType linear_gradient(
        "x", "0.5 * x[0] * x[0] - x[0]", 2, "affine gradient", {{"x[0] - 1.0", "0.0", "0.0"}});
    this->check(compute(linear_gradient), factor_value_ * 1.0 / 3.0, 1.71e-13);
  }

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void correct_for_quadratic_arguments() const
  {
    const ExpressionFunctionType quadratic_gradient(
        "x", "(1.0/3.0) * x[0] * x[0] * x[0]", 3, ", quadratic gradient", {{"x[0]*x[0]", "0.0", "0.0"}});
    this->check(compute(quadratic_gradient), factor_value_ * 1.0 / 5.0, 5.33e-13);
  }

  void is_matrix_operator()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& space            = this->space_;

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
  using typename OperatorBaseType::GridViewType;
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
    const auto& grid_view        = this->space_.grid_view();

    EllipticOperator<ExpressionFunctionType, void, GridViewType> DUNE_UNUSED(only_factor)(grid_view, diffusion_factor);
    EllipticOperator<ExpressionFunctionType, void, GridViewType> DUNE_UNUSED(only_factor_w_over)(
        1, grid_view, diffusion_factor);

    EllipticOperator<TensorFunctionType, void, GridViewType> DUNE_UNUSED(only_tensor)(grid_view, diffusion_tensor);
    EllipticOperator<TensorFunctionType, void, GridViewType> DUNE_UNUSED(only_tensor_w_over)(
        1, grid_view, diffusion_tensor);

    EllipticOperator<ExpressionFunctionType, TensorFunctionType, GridViewType> DUNE_UNUSED(both)(
        grid_view, diffusion_factor, diffusion_tensor);
    EllipticOperator<ExpressionFunctionType, TensorFunctionType, GridViewType> DUNE_UNUSED(both_w_over)(
        1, grid_view, diffusion_factor, diffusion_tensor);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();

    auto DUNE_UNUSED(only_factor) = make_elliptic_operator(grid_view, diffusion_factor);
    auto DUNE_UNUSED(only_factor_w_over) = make_elliptic_operator(grid_view, diffusion_factor, 1);

    auto DUNE_UNUSED(only_tensor) = make_elliptic_operator(grid_view, diffusion_tensor);
    auto DUNE_UNUSED(only_tensor_w_over) = make_elliptic_operator(grid_view, diffusion_tensor, 1);

    auto DUNE_UNUSED(both) = make_elliptic_operator(grid_view, diffusion_factor, diffusion_tensor);
    auto DUNE_UNUSED(both_w_over) = make_elliptic_operator(grid_view, diffusion_factor, diffusion_tensor, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();

    return make_elliptic_operator(grid_view, diffusion_factor, diffusion_tensor)->apply2(function, function);
  } // ... compute(...)

  void apply_is_callable()
  {
    const auto& diffusion_factor = this->factor_;
    const auto& diffusion_tensor = this->tensor_function_;
    const auto& grid_view        = this->space_.grid_view();
    auto& source                 = this->discrete_function_;
    auto range                   = make_discrete_function<VectorType>(this->space_);

    auto op = make_elliptic_operator(grid_view, diffusion_factor, diffusion_tensor);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct EllipticOperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_ELLIPTIC_HH
