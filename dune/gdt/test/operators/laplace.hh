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

#ifndef DUNE_GDT_TEST_OPERATORS_LAPLACE_HH
#define DUNE_GDT_TEST_OPERATORS_LAPLACE_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/test/float_cmp.hh>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/laplace.hh>

#include "elliptic.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
struct LaplaceLocalizableProductTest : public EllipticProductBase<SpaceType>, public LocalizableProductBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> EllipticBaseType;
  typedef LocalizableProductBase<SpaceType> LocalizableBaseType;
  using typename LocalizableBaseType::GridViewType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename LocalizableBaseType::ScalarFunctionType;
  using typename LocalizableBaseType::RangeFieldType;

  LaplaceLocalizableProductTest()
    : EllipticBaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef LaplaceLocalizableProduct<GridViewType, ScalarFunctionType, ScalarFunctionType, double> CtorTestProductType;
    DUNE_UNUSED CtorTestProductType wo_over_integrate(grid_view, range, source);
    DUNE_UNUSED CtorTestProductType with_over_integrate(1, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto wo_over_integrate DUNE_UNUSED = make_laplace_localizable_product(grid_view, range, source);
    auto with_over_integrate DUNE_UNUSED = make_laplace_localizable_product(grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_view = this->space_.grid_view();

    auto product = make_laplace_localizable_product(grid_view, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_laplace_localizable_product(grid_view, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    DXTC_EXPECT_FLOAT_EQ(result_tbb, result, 1e-14);
    return result;
  } // ... compute(...)

  void is_localizable_product()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto product = make_laplace_localizable_product(grid_view, range, source);
    this->localizable_product_test(*product);
  } // ... is_localizable_product(...)
}; // struct LaplaceLocalizableProductTest


/**
 * \note Assumes that Operators::Projection does the right thing!
 */
template <class SpaceType>
struct LaplaceMatrixOperatorTest : public EllipticMatrixOperatorTest<SpaceType>
{
  typedef EllipticMatrixOperatorTest<SpaceType> BaseType;
  using typename BaseType::GridViewType;
  using typename BaseType::ExpressionFunctionType;
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::MatrixType;

  LaplaceMatrixOperatorTest()
    : BaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& space = this->space_;
    const auto& grid_view = this->space_.grid_view();

    // without matrix
    //   without over_integrate
    //     simplified argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_1(space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_2(space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_3(space, space, grid_view);
    //     full argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_4(space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_5(space,
                                                                                                          grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_6(
        space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_7(1, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_8(1, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> no_matrix_9(1, space, space, grid_view);
    //     full argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_10(1, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_11(
        1, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> no_matrix_12(
        1, space, space, grid_view);
    // with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //   without over_integrate
    //     simplified argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_1(matrix, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_2(matrix, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_3(matrix, space, space, grid_view);
    //     full argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_4(matrix, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_5(
        matrix, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_6(
        matrix, space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_7(1, matrix, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_8(1, matrix, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType> matrix_9(1, matrix, space, space, grid_view);
    //     full argument list
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_10(
        1, matrix, space);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_11(
        1, matrix, space, grid_view);
    DUNE_UNUSED LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType, SpaceType, double> matrix_12(
        1, matrix, space, space, grid_view);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& space = this->space_;
    const auto& grid_view = this->space_.grid_view();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // without matrix
    auto op_01 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space);
    auto op_02 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space, 1);
    auto op_03 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space, grid_view);
    auto op_04 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space, grid_view, 1);
    auto op_05 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space, space, grid_view);
    auto op_06 DUNE_UNUSED = make_laplace_matrix_operator<MatrixType>(space, space, grid_view, 1);
    // with matrix
    auto op_07 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space);
    auto op_08 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space, 1);
    auto op_09 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space, grid_view);
    auto op_10 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space, grid_view, 1);
    auto op_11 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space, space, grid_view);
    auto op_12 DUNE_UNUSED = make_laplace_matrix_operator(matrix, space, space, grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& space = this->space_;
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(function, discrete_function);
    // compute product
    auto product = make_laplace_matrix_operator<MatrixType>(space);
    const auto result = product->apply2(discrete_function, discrete_function);

    auto product_tbb = make_laplace_matrix_operator<MatrixType>(space);
    product_tbb->assemble(true);
    const auto result_tbb = product_tbb->apply2(discrete_function, discrete_function);
    DXTC_EXPECT_FLOAT_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_matrix_operator()
  {
    const auto& space = this->space_;

    auto op = make_laplace_matrix_operator<MatrixType>(space);
    this->matrix_operator_test(*op);
  } // ... is_matrix_operator(...)
}; // struct LaplaceMatrixOperatorTest


template <class SpaceType>
struct LaplaceOperatorTest : public EllipticProductBase<SpaceType>, public OperatorBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> EllipticBaseType;
  typedef OperatorBase<SpaceType> OperatorBaseType;
  using typename OperatorBaseType::GridViewType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::VectorType;

  LaplaceOperatorTest()
    : EllipticBaseType(1.)
  {
  }

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();

    DUNE_UNUSED LaplaceOperator<GridViewType> wo_over_integrate(grid_view);
    DUNE_UNUSED LaplaceOperator<GridViewType> with_over_integrate(grid_view, 1);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();

    auto wo_over_integrate DUNE_UNUSED = make_laplace_operator(grid_view);
    auto with_over_integrate DUNE_UNUSED = make_laplace_operator(grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_view = this->space_.grid_view();

    return make_laplace_operator(grid_view)->apply2(function, function);
  }

  void apply_is_callable()
  {
    const auto& grid_view = this->space_.grid_view();
    auto& source = this->discrete_function_;
    auto range = make_discrete_function<VectorType>(this->space_);

    auto op = make_laplace_operator(grid_view);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct LaplaceOperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_LAPLACE_HH
