// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_OPERATORS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_L2_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/test/float_cmp.hh>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/l2.hh>

#include "weighted-l2.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct L2LocalizableProductTest : public WeightedL2ProductBase<SpaceType>, public LocalizableProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  using typename LocalizableBaseType::GridLayerType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using typename LocalizableBaseType::ScalarFunctionType;
  using typename LocalizableBaseType::RangeFieldType;

  L2LocalizableProductTest()
    : WeightedL2BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef L2LocalizableProduct<GridLayerType, ScalarFunctionType, ScalarFunctionType, double> CtorTestProductType;
    DUNE_UNUSED CtorTestProductType wo_over_integrate(grid_layer, range, source);
    DUNE_UNUSED CtorTestProductType with_over_integrate(1, grid_layer, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto wo_over_integrate DUNE_UNUSED = make_l2_localizable_product(grid_layer, range, source);
    auto with_over_integrate DUNE_UNUSED = make_l2_localizable_product(grid_layer, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_layer = this->space_.grid_layer();

    auto product = make_l2_localizable_product(grid_layer, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_l2_localizable_product(grid_layer, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    DXTC_EXPECT_FLOAT_EQ(result_tbb, result, 1e-14);
    return result;
  } // ... compute(...)

  void is_localizable_product()
  {
    const auto& grid_layer = this->space_.grid_layer();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto product = make_l2_localizable_product(grid_layer, range, source);
    this->localizable_product_test(*product);
  } // ... is_localizable_product(...)
}; // struct L2LocalizableProductTest


/**
 * \note Assumes that Operators::Projection does the right thing!
 */
template <class SpaceType>
struct L2MatrixOperatorTest : public WeightedL2ProductBase<SpaceType>, public MatrixOperatorBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef MatrixOperatorBase<SpaceType> MatrixBaseType;
  using typename MatrixBaseType::GridLayerType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using typename MatrixBaseType::DiscreteFunctionType;
  using typename MatrixBaseType::RangeFieldType;
  using typename MatrixBaseType::MatrixType;

  L2MatrixOperatorTest()
    : WeightedL2BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();

    // without matrix
    //   without over_integrate
    //     simplified argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_1(space);
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_2(space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_3(space, space, grid_layer);
    //     full argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_4(space);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_5(space,
                                                                                                      grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_6(
        space, space, grid_layer);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_7(1, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_8(1, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType> no_matrix_9(1, space, space, grid_layer);
    //     full argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_10(1, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_11(
        1, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> no_matrix_12(
        1, space, space, grid_layer);
    // with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //   without over_integrate
    //     simplified argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_1(matrix, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_2(matrix, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_3(matrix, space, space, grid_layer);
    //     full argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_4(matrix, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_5(
        matrix, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_6(
        matrix, space, space, grid_layer);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_7(1, matrix, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_8(1, matrix, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType> matrix_9(1, matrix, space, space, grid_layer);
    //     full argument list
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_10(1, matrix, space);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_11(
        1, matrix, space, grid_layer);
    DUNE_UNUSED L2MatrixOperator<SpaceType, MatrixType, GridLayerType, SpaceType, double> matrix_12(
        1, matrix, space, space, grid_layer);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& space = this->space_;
    const auto& grid_layer = this->space_.grid_layer();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // without matrix
    auto op_01 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space);
    auto op_02 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space, 1);
    auto op_03 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space, grid_layer);
    auto op_04 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space, grid_layer, 1);
    auto op_05 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space, space, grid_layer);
    auto op_06 DUNE_UNUSED = make_l2_matrix_operator<MatrixType>(space, space, grid_layer, 1);
    // with matrix
    auto op_07 DUNE_UNUSED = make_l2_matrix_operator(matrix, space);
    auto op_08 DUNE_UNUSED = make_l2_matrix_operator(matrix, space, 1);
    auto op_09 DUNE_UNUSED = make_l2_matrix_operator(matrix, space, grid_layer);
    auto op_10 DUNE_UNUSED = make_l2_matrix_operator(matrix, space, grid_layer, 1);
    auto op_11 DUNE_UNUSED = make_l2_matrix_operator(matrix, space, space, grid_layer);
    auto op_12 DUNE_UNUSED = make_l2_matrix_operator(matrix, space, space, grid_layer, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& space = this->space_;
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(function, discrete_function);
    // compute product
    auto product = make_l2_matrix_operator<MatrixType>(space);
    const auto result = product->apply2(discrete_function, discrete_function);

    auto product_tbb = make_l2_matrix_operator<MatrixType>(space);
    product_tbb->assemble(true);
    const auto result_tbb = product_tbb->apply2(discrete_function, discrete_function);
    DXTC_EXPECT_FLOAT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_matrix_operator()
  {
    const auto& space = this->space_;

    auto op = make_l2_matrix_operator<MatrixType>(space);
    this->matrix_operator_test(*op);
  } // ... is_matrix_operator(...)
}; // struct L2MatrixOperatorTest


template <class SpaceType>
struct L2OperatorTest : public WeightedL2ProductBase<SpaceType>, public OperatorBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef OperatorBase<SpaceType> OperatorBaseType;
  using typename OperatorBaseType::GridLayerType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::VectorType;

  L2OperatorTest()
    : WeightedL2BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_layer = this->space_.grid_layer();

    DUNE_UNUSED L2Operator<GridLayerType> wo_over_integrate(grid_layer);
    DUNE_UNUSED L2Operator<GridLayerType> with_over_integrate(grid_layer, 1);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_layer = this->space_.grid_layer();

    auto wo_over_integrate DUNE_UNUSED = make_l2_operator(grid_layer);
    auto with_over_integrate DUNE_UNUSED = make_l2_operator(grid_layer, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_layer = this->space_.grid_layer();

    return make_l2_operator(grid_layer)->apply2(function, function);
  }

  void apply_is_callable()
  {
    const auto& grid_layer = this->space_.grid_layer();
    auto& source = this->discrete_function_;
    auto range = make_discrete_function<VectorType>(this->space_);

    auto op = make_l2_operator(grid_layer);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct L2OperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_L2_HH
