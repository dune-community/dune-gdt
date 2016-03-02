// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_LAPLACE_HH
#define DUNE_GDT_TEST_OPERATORS_LAPLACE_HH

#include <dune/stuff/common/string.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/laplace.hh>

#include "elliptic.hh"

namespace Dune {
namespace GDT {
namespace Test {


template< class SpaceType >
struct LaplaceLocalizableProductTest
  : public EllipticProductBase< SpaceType >
  , public LocalizableProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType >    EllipticBaseType;
  typedef LocalizableProductBase< SpaceType > LocalizableBaseType;
  using typename LocalizableBaseType::GridViewType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename LocalizableBaseType::ScalarFunctionType;
  using typename LocalizableBaseType::RangeFieldType;

  LaplaceLocalizableProductTest()
    : EllipticBaseType(1.)
  {}

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    typedef LaplaceLocalizableProduct< GridViewType, ScalarFunctionType, ScalarFunctionType, double >
        CtorTestProductType;
    CtorTestProductType DUNE_UNUSED(wo_over_integrate)(     grid_view, range, source);
    CtorTestProductType DUNE_UNUSED(with_over_integrate)(1, grid_view, range, source);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();
    const auto& source = this->scalar_function_;
    const auto& range = this->scalar_function_;

    auto DUNE_UNUSED(wo_over_integrate)   = make_laplace_localizable_product(grid_view, range, source   );
    auto DUNE_UNUSED(with_over_integrate) = make_laplace_localizable_product(grid_view, range, source, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& grid_view = this->space_.grid_view();

    auto product = make_laplace_localizable_product(grid_view, function, function);
    const auto result = product->apply2();

    auto product_tbb = make_laplace_localizable_product(grid_view, function, function);
    product_tbb->walk(true);
    const auto result_tbb = product_tbb->apply2();

    EXPECT_EQ(result_tbb, result);
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
template< class SpaceType >
struct LaplaceMatrixOperatorTest
  : public EllipticMatrixOperatorTest< SpaceType >
{
  typedef EllipticMatrixOperatorTest< SpaceType > BaseType;
  using typename BaseType::GridViewType;
  using typename BaseType::ExpressionFunctionType;
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::MatrixType;

  LaplaceMatrixOperatorTest()
    : BaseType(1.)
  {}

  void constructible_by_ctor()
  {
    const auto& space = this->space_;
    const auto& grid_view = this->space_.grid_view();

    // without matrix
    //   without over_integrate
    //     simplified argument list
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_1)(space);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_2)(space, grid_view);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_3)(space, space, grid_view);
    //     full argument list
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_4)(space);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_5)(space, grid_view);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_6)(space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_7)(1, space);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_8)(1, space, grid_view);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(no_matrix_9)(1, space, space, grid_view);
    //     full argument list
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_10)(1, space);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_11)(1, space, grid_view);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(no_matrix_12)(1, space, space, grid_view);
    // with matrix
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());
    //   without over_integrate
    //     simplified argument list
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_1)(matrix, space);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_2)(matrix, space, grid_view);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_3)(matrix, space, space, grid_view);
    //     full argument list
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_4)(matrix, space);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_5)(matrix, space, grid_view);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_6)(matrix, space, space, grid_view);
    //   with over_integrate
    //     simplified argument list
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_7)(1, matrix, space);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_8)(1, matrix, space, grid_view);
    LaplaceMatrixOperator< SpaceType > DUNE_UNUSED(matrix_9)(1, matrix, space, space, grid_view);
    //     full argument list
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_10)(1, matrix, space);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_11)(1, matrix, space, grid_view);
    LaplaceMatrixOperator< SpaceType, MatrixType, GridViewType, SpaceType, double >
        DUNE_UNUSED(matrix_12)(1, matrix, space, space, grid_view);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& space = this->space_;
    const auto& grid_view = this->space_.grid_view();
    MatrixType matrix(space.mapper().size(), space.mapper().size(), space.compute_volume_pattern());

    // without matrix
    auto DUNE_UNUSED(op_01) = make_laplace_matrix_operator< MatrixType >(space   );
    auto DUNE_UNUSED(op_02) = make_laplace_matrix_operator< MatrixType >(space, 1);
    auto DUNE_UNUSED(op_03) = make_laplace_matrix_operator< MatrixType >(space, grid_view   );
    auto DUNE_UNUSED(op_04) = make_laplace_matrix_operator< MatrixType >(space, grid_view, 1);
    auto DUNE_UNUSED(op_05) = make_laplace_matrix_operator< MatrixType >(space, space, grid_view   );
    auto DUNE_UNUSED(op_06) = make_laplace_matrix_operator< MatrixType >(space, space, grid_view, 1);
    // with matrix
    auto DUNE_UNUSED(op_07) = make_laplace_matrix_operator(matrix, space   );
    auto DUNE_UNUSED(op_08) = make_laplace_matrix_operator(matrix, space, 1);
    auto DUNE_UNUSED(op_09) = make_laplace_matrix_operator(matrix, space, grid_view   );
    auto DUNE_UNUSED(op_10) = make_laplace_matrix_operator(matrix, space, grid_view, 1);
    auto DUNE_UNUSED(op_11) = make_laplace_matrix_operator(matrix, space, space, grid_view   );
    auto DUNE_UNUSED(op_12) = make_laplace_matrix_operator(matrix, space, space, grid_view, 1);
  } // ... constructible_by_factory()

  virtual RangeFieldType compute(const ExpressionFunctionType& function) const override final
  {
    const auto& space = this->space_;
    // project the function
    DiscreteFunctionType discrete_function(space);
    project(function, discrete_function);
    // compute product
    auto product = make_laplace_matrix_operator< MatrixType >(space);
    const auto result = product->apply2(discrete_function, discrete_function);

    auto product_tbb = make_laplace_matrix_operator< MatrixType >(space);
    product_tbb->assemble(true);
    const auto result_tbb = product_tbb->apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void is_matrix_operator()
  {
    const auto& space = this->space_;

    auto op = make_laplace_matrix_operator< MatrixType >(space);
    this->matrix_operator_test(*op);
  } // ... is_matrix_operator(...)
}; // struct LaplaceMatrixOperatorTest


template< class SpaceType >
struct LaplaceOperatorTest
  : public EllipticProductBase< SpaceType >
  , public OperatorBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > EllipticBaseType;
  typedef OperatorBase< SpaceType > OperatorBaseType;
  using typename OperatorBaseType::GridViewType;
  using typename EllipticBaseType::ExpressionFunctionType;
  using typename OperatorBaseType::RangeFieldType;
  using typename OperatorBaseType::VectorType;

  LaplaceOperatorTest()
    : EllipticBaseType(1.)
  {}

  void constructible_by_ctor()
  {
    const auto& grid_view = this->space_.grid_view();

    LaplaceOperator< GridViewType > DUNE_UNUSED(wo_over_integrate)(  grid_view   );
    LaplaceOperator< GridViewType > DUNE_UNUSED(with_over_integrate)(grid_view, 1);
  } // ... constructible_by_ctor(...)

  void constructible_by_factory()
  {
    const auto& grid_view = this->space_.grid_view();

    auto DUNE_UNUSED(wo_over_integrate)   = make_laplace_operator(grid_view   );
    auto DUNE_UNUSED(with_over_integrate) = make_laplace_operator(grid_view, 1);
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
    auto range = make_discrete_function< VectorType >(this->space_);

    auto op = make_laplace_operator(grid_view);
    op->apply(source, range);
  } // ... apply_is_callable(...)
}; // struct LaplaceOperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_LAPLACE_HH
