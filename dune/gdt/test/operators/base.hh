// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_OPERATORS_BASE_HH
#define DUNE_GDT_TEST_OPERATORS_BASE_HH

#include <dune/common/unused.hh>

#include <dune/xt/functions/expression.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/tools.hh>

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


template <class SpaceType>
struct OperatorBaseTraits
{
  static_assert(is_space<SpaceType>::value, "");
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const size_t dimRange       = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1>
      ScalarFunctionType;
  typedef Dune::XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::XT::Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain,
                                                dimDomain>
      TensorFunctionType;
  typedef typename XT::LA::Container<RangeFieldType, XT::LA::default_backend>::MatrixType MatrixType;
  typedef typename XT::LA::Container<RangeFieldType, XT::LA::default_backend>::VectorType VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
}; // class OperatorBaseTraits


} // namespace internal


template <class SpaceType>
struct OperatorBase : public ::testing::Test
{
  typedef internal::OperatorBaseTraits<SpaceType> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::XT::Grid::GridProvider<GridType> GridProviderType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::ScalarFunctionType ScalarFunctionType;
  typedef typename Traits::FunctionType FunctionType;
  typedef typename Traits::TensorFunctionType TensorFunctionType;
  typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::VectorType VectorType;
  static const size_t dimDomain = Traits::dimDomain;

  OperatorBase()
    : grid_provider_(XT::Grid::make_cube_grid<GridType>(0.0, 1.0, 3u))
    , space_(Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid_provider_.grid()))
    , scalar_function_("x", "x[0]", 1, "scalar function", {{"1.0", "0.0", "0.0"}})
    , function_("x", {"x[0]", "0", "0"}, 1)
    , tensor_function_(XT::Functions::internal::UnitMatrix<RangeFieldType, dimDomain>::value())
    , discrete_function_(space_)
  {
  }

  GridProviderType grid_provider_;
  const SpaceType space_;
  const ScalarFunctionType scalar_function_;
  const FunctionType function_;
  const TensorFunctionType tensor_function_;
  DiscreteFunctionType discrete_function_;
}; // class OperatorBase


template <class SpaceType>
struct LocalizableProductBase : public OperatorBase<SpaceType>
{
  typedef OperatorBase<SpaceType> BaseType;
  using typename BaseType::GridViewType;

  template <class ProductImp>
  void localizable_product_test(ProductImp& prod)
  {
    const auto& DUNE_UNUSED(source)    = prod.source();
    const auto& DUNE_UNUSED(range)     = prod.range();
    auto& DUNE_UNUSED(non_const_range) = prod.range();

    XT::Grid::Walker<GridViewType> walker(this->space_.grid_view());
    walker.add(prod);
    walker.walk();

    auto DUNE_UNUSED(result) = prod.apply2();
  } // ... localizable_product_test(...)
}; // class LocalizableProductBase


template <class SpaceType>
struct MatrixOperatorBase : public OperatorBase<SpaceType>
{
  typedef OperatorBase<SpaceType> BaseType;
  using typename BaseType::GridViewType;
  using typename BaseType::MatrixType;

  template <class OperatorImp>
  void matrix_operator_test(OperatorImp& op)
  {
    const auto& DUNE_UNUSED(matrix)       = op.matrix();
    auto& DUNE_UNUSED(non_const_matrix)   = op.matrix();
    const auto& DUNE_UNUSED(source_space) = op.source_space();
    const auto& DUNE_UNUSED(range_space)  = op.range_space();

    XT::Grid::Walker<GridViewType> walker(this->space_.grid_view());
    walker.add(op);
    walker.walk();
  } // ... matrix_operator_test(...)
}; // class LocalizableProductBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_BASE_HH
