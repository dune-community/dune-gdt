// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/test/float_cmp.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/weighted-l2.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \note The values in correct_for_constant_arguments(), etc., are valid for the d-dimendional unit cube.
 */
template <class SpaceType>
struct WeightedL2ProductBase
{
  typedef typename SpaceType::GridLayerType GridLayerType;
  typedef XT::Grid::extract_grid_t<GridLayerType> GridType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1>
      ExpressionFunctionType;

  WeightedL2ProductBase(const double weight_value = 42)
    : weight_value_(weight_value)
    , weight_("x", Dune::XT::Common::to_string(double(weight_value_)), 0) // linker error if double(...) is missing
    , constant_("x", "1.0", 0)
    , linear_("x", "x[0] - 1.0", 1)
    , quadratic_("x", "x[0]*x[0]", 2)
  {}

  virtual ~WeightedL2ProductBase() = default;

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const = 0;

  void correct_for_constant_arguments(const RangeFieldType epsilon = 1e-15) const
  {
    check(compute(constant_), weight_value_ * 1.0, epsilon);
  }

  void correct_for_linear_arguments(const RangeFieldType epsilon = 1e-15) const
  {
    check(compute(linear_), weight_value_ * (1.0 / 3.0), epsilon);
  }

  void correct_for_quadratic_arguments(const RangeFieldType epsilon = 1e-15) const
  {
    check(compute(quadratic_), weight_value_ * (1.0 / 5.0), epsilon);
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon) const
  {
    //! might be off, since atol was used before
    DXTC_EXPECT_FLOAT_LE(expected, result, epsilon);
  } // ... check(...)

  const double weight_value_;
  const ExpressionFunctionType weight_;
  const ExpressionFunctionType constant_;
  const ExpressionFunctionType linear_;
  const ExpressionFunctionType quadratic_;
}; // struct WeightedL2ProductBase


template <class SpaceType>
struct WeightedL2LocalizableProductTest
  : public WeightedL2ProductBase<SpaceType>
  , public LocalizableProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  using typename LocalizableBaseType::GridLayerType;
  using typename LocalizableBaseType::RangeFieldType;
  using typename LocalizableBaseType::ScalarFunctionType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using WeightedL2BaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef WeightedL2LocalizableProduct<ExpressionFunctionType,
                                         GridLayerType,
                                         ScalarFunctionType,
                                         ScalarFunctionType,
                                         double>
        CtorTestProductType;
    DUNE_UNUSED CtorTestProductType wo_over_integrate(weight, grid_layer, range, source);
    DUNE_UNUSED CtorTestProductType with_over_integrate(1, weight, grid_layer, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto wo_over_integrate DUNE_UNUSED = make_weighted_l2_localizable_product(weight, grid_layer, range, source);
    auto with_over_integrate DUNE_UNUSED = make_weighted_l2_localizable_product(weight, grid_layer, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();

    auto product = make_weighted_l2_localizable_product(weight, grid_layer, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_weighted_l2_localizable_product(weight, grid_layer, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    DXTC_EXPECT_FLOAT_EQ(result_tbb, result, 1e-14);
    return result;
  } // ... compute(...)

  void is_localizable_product()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto product = make_weighted_l2_localizable_product(weight, grid_layer, range, source);
    this->localizable_product_test(*product);
  } // ... is_localizable_product(...)
}; // struct WeightedL2LocalizableProductTest


/**
 * \note Assumes that Operators::Projection does the right thing!
 */
template <class SpaceType>
struct WeightedL2MatrixOperatorTest
  : public WeightedL2ProductBase<SpaceType>
  , public MatrixOperatorBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef MatrixOperatorBase<SpaceType> MatrixBaseType;
  using typename MatrixBaseType::DiscreteFunctionType;
  using typename MatrixBaseType::GridLayerType;
  using typename MatrixBaseType::MatrixType;
  using typename MatrixBaseType::RangeFieldType;
  using typename MatrixBaseType::ScalarFunctionType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using WeightedL2BaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& weight = this->weight_;
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();

    // without matrix
    //   without over_integrate
    //     simplified argument list
    using ThreeArgsOp = WeightedL2MatrixOperator<ExpressionFunctionType, SpaceType, MatrixType>;
    DUNE_UNUSED ThreeArgsOp no_matrix_1(weight, space);
    DUNE_UNUSED ThreeArgsOp no_matrix_2(weight, space, grid_layer);
    DUNE_UNUSED ThreeArgsOp no_matrix_3(weight, space, space, grid_layer);
    //     full argument list
    using SixArgsOp =
        WeightedL2MatrixOperator<ExpressionFunctionType, SpaceType, MatrixType, GridLayerType, SpaceType, double>;
    DUNE_UNUSED SixArgsOp no_matrix_4(weight, space);
    SixArgsOp DUNE_UNUSED no_matrix_5(weight, space, grid_layer);
    SixArgsOp DUNE_UNUSED no_matrix_6(weight, space, space, grid_layer);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED ThreeArgsOp no_matrix_7(1, weight, space);
    DUNE_UNUSED ThreeArgsOp no_matrix_8(1, weight, space, grid_layer);
    DUNE_UNUSED ThreeArgsOp no_matrix_9(1, weight, space, space, grid_layer);
    //     full argument list
    SixArgsOp DUNE_UNUSED no_matrix_10(1, weight, space);
    SixArgsOp DUNE_UNUSED no_matrix_11(1, weight, space, grid_layer);
    SixArgsOp DUNE_UNUSED no_matrix_12(1, weight, space, space, grid_layer);
    // with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //   without over_integrate
    //     simplified argument list
    DUNE_UNUSED ThreeArgsOp matrix_1(weight, matrix, space);
    DUNE_UNUSED ThreeArgsOp matrix_2(weight, matrix, space, grid_layer);
    DUNE_UNUSED ThreeArgsOp matrix_3(weight, matrix, space, space, grid_layer);
    //     full argument list
    SixArgsOp DUNE_UNUSED matrix_4(weight, matrix, space);
    SixArgsOp DUNE_UNUSED matrix_5(weight, matrix, space, grid_layer);
    SixArgsOp DUNE_UNUSED matrix_6(weight, matrix, space, space, grid_layer);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED ThreeArgsOp matrix_7(1, weight, matrix, space);
    DUNE_UNUSED ThreeArgsOp matrix_8(1, weight, matrix, space, grid_layer);
    DUNE_UNUSED ThreeArgsOp matrix_9(1, weight, matrix, space, space, grid_layer);
    //     full argument list
    SixArgsOp DUNE_UNUSED matrix_10(1, weight, matrix, space);
    SixArgsOp DUNE_UNUSED matrix_11(1, weight, matrix, space, grid_layer);
    SixArgsOp DUNE_UNUSED matrix_12(1, weight, matrix, space, space, grid_layer);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& weight = this->weight_;
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // without matrix
    auto op_01 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space);
    auto op_02 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space, 1);
    auto op_03 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space, grid_layer);
    auto op_04 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space, grid_layer, 1);
    auto op_05 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space, space, grid_layer);
    auto op_06 DUNE_UNUSED = make_weighted_l2_matrix_operator<MatrixType>(weight, space, space, grid_layer, 1);
    // with matrix
    auto op_07 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space);
    auto op_08 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space, 1);
    auto op_09 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space, grid_layer);
    auto op_10 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space, grid_layer, 1);
    auto op_11 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space, space, grid_layer);
    auto op_12 DUNE_UNUSED = make_weighted_l2_matrix_operator(weight, matrix, space, space, grid_layer, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& weight = this->weight_;
    const auto& space = this->space_;
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(function, discrete_function);
    // compute product
    auto product = make_weighted_l2_matrix_operator<MatrixType>(weight, space);
    const auto result = product->apply2(discrete_function, discrete_function);

    auto product_tbb = make_weighted_l2_matrix_operator<MatrixType>(weight, space);
    product_tbb->assemble(true);
    const auto result_tbb = product_tbb->apply2(discrete_function, discrete_function);
    DXTC_EXPECT_FLOAT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_matrix_operator()
  {
    const auto& weight = this->weight_;
    const auto& space = this->space_;

    auto op = make_weighted_l2_matrix_operator<MatrixType>(weight, space);
    this->matrix_operator_test(*op);
  } // ... is_matrix_operator(...)
}; // struct WeightedL2MatrixOperatorTest


template <class SpaceType>
struct WeightedL2OperatorTest
  : public WeightedL2ProductBase<SpaceType>
  , public OperatorBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef OperatorBase<SpaceType> OperatorBaseType;
  using typename OperatorBaseType::DiscreteFunctionType;
  using typename OperatorBaseType::GridLayerType;
  using typename OperatorBaseType::MatrixType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::ScalarFunctionType;
  using typename OperatorBaseType::VectorType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using WeightedL2BaseType::dimDomain;

  void constructible_by_ctor()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();

    DUNE_UNUSED WeightedL2Operator<ExpressionFunctionType, GridLayerType, MatrixType> wo_over_integrate(weight,
                                                                                                        grid_layer);
    DUNE_UNUSED WeightedL2Operator<ExpressionFunctionType, GridLayerType, MatrixType> with_over_integrate(
        weight, grid_layer, 1);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();

    auto wo_over_integrate DUNE_UNUSED = make_weighted_l2_operator(grid_layer, weight);
    auto with_over_integrate DUNE_UNUSED = make_weighted_l2_operator(grid_layer, weight, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();

    return make_weighted_l2_operator(grid_layer, weight)->apply2(function, function);
  }

  void apply_is_callable()
  {
    const auto& weight = this->weight_;
    const auto& grid_layer = this->space_.grid_layer();
    auto& source = this->discrete_function_;
    auto range = make_discrete_function<VectorType>(this->space_);

    auto op = make_weighted_l2_operator(grid_layer, weight);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct WeightedL2OperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_WEIGHTED_L2_HH
