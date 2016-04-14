// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_OPERATORS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_L2_HH

#include <dune/stuff/common/string.hh>
#include <dune/stuff/test/gtest/gtest.h>

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


/**
 * \note Assumes that Operators::Projection does the right thing!
 */
template <class SpaceType>
struct L2MatrixOperatorTest : public WeightedL2ProductBase<SpaceType>, public MatrixOperatorBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> WeightedL2BaseType;
  typedef MatrixOperatorBase<SpaceType> MatrixBaseType;
  using typename MatrixBaseType::GridViewType;
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
    const auto& space     = this->space_;
    const auto& grid_view = this->space_.grid_view();

    // without matrix
    //   without over_integrate
    //     simplified argument list
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_1)(space);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_2)(space, grid_view);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_3)(space, space, grid_view);
    //     full argument list
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_4)(space);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_5)(space, grid_view);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_6)(
        space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_7)(1, space);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_8)(1, space, grid_view);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(no_matrix_9)(1, space, space, grid_view);
    //     full argument list
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_10)(1, space);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_11)(
        1, space, grid_view);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(no_matrix_12)(
        1, space, space, grid_view);
    // with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //   without over_integrate
    //     simplified argument list
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_1)(matrix, space);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_2)(matrix, space, grid_view);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_3)(matrix, space, space, grid_view);
    //     full argument list
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_4)(matrix, space);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_5)(
        matrix, space, grid_view);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_6)(
        matrix, space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_7)(1, matrix, space);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_8)(1, matrix, space, grid_view);
    L2MatrixOperator<SpaceType> DUNE_UNUSED(matrix_9)(1, matrix, space, space, grid_view);
    //     full argument list
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_10)(1, matrix, space);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_11)(
        1, matrix, space, grid_view);
    L2MatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> DUNE_UNUSED(matrix_12)(
        1, matrix, space, space, grid_view);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& space     = this->space_;
    const auto& grid_view = this->space_.grid_view();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // without matrix
    auto DUNE_UNUSED(op_01) = make_l2_matrix_operator<MatrixType>(space);
    auto DUNE_UNUSED(op_02) = make_l2_matrix_operator<MatrixType>(space, 1);
    auto DUNE_UNUSED(op_03) = make_l2_matrix_operator<MatrixType>(space, grid_view);
    auto DUNE_UNUSED(op_04) = make_l2_matrix_operator<MatrixType>(space, grid_view, 1);
    auto DUNE_UNUSED(op_05) = make_l2_matrix_operator<MatrixType>(space, space, grid_view);
    auto DUNE_UNUSED(op_06) = make_l2_matrix_operator<MatrixType>(space, space, grid_view, 1);
    // with matrix
    auto DUNE_UNUSED(op_07) = make_l2_matrix_operator(matrix, space);
    auto DUNE_UNUSED(op_08) = make_l2_matrix_operator(matrix, space, 1);
    auto DUNE_UNUSED(op_09) = make_l2_matrix_operator(matrix, space, grid_view);
    auto DUNE_UNUSED(op_10) = make_l2_matrix_operator(matrix, space, grid_view, 1);
    auto DUNE_UNUSED(op_11) = make_l2_matrix_operator(matrix, space, space, grid_view);
    auto DUNE_UNUSED(op_12) = make_l2_matrix_operator(matrix, space, space, grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& space = this->space_;
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(function, discrete_function);
    // compute product
    auto product      = make_l2_matrix_operator<MatrixType>(space);
    const auto result = product->apply2(discrete_function, discrete_function);

    auto product_tbb = make_l2_matrix_operator<MatrixType>(space);
    product_tbb->assemble(true);
    const auto result_tbb = product_tbb->apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
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
  using typename OperatorBaseType::GridViewType;
  using typename WeightedL2BaseType::ExpressionFunctionType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::VectorType;

  L2OperatorTest()
    : WeightedL2BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();

    L2Operator<GridViewType> DUNE_UNUSED(wo_over_integrate)(grid_view);
    L2Operator<GridViewType> DUNE_UNUSED(with_over_integrate)(grid_view, 1);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();

    auto DUNE_UNUSED(wo_over_integrate) = make_l2_operator(grid_view);
    auto DUNE_UNUSED(with_over_integrate) = make_l2_operator(grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_view = this->space_.grid_view();

    return make_l2_operator(grid_view)->apply2(function, function);
  }

  void apply_is_callable()
  {
    const auto& grid_view = this->space_.grid_view();
    auto& source          = this->discrete_function_;
    auto range            = make_discrete_function<VectorType>(this->space_);

    auto op = make_l2_operator(grid_view);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct L2OperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_L2_HH
