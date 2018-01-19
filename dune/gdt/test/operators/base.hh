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

#ifndef DUNE_GDT_TEST_OPERATORS_BASE_HH
#define DUNE_GDT_TEST_OPERATORS_BASE_HH

#include <dune/common/unused.hh>

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/la/container.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/expression.hh>
#include <dune/xt/functions/constant.hh>

#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


template <class SpaceType>
struct OperatorBaseTraits
{
  static_assert(is_space<SpaceType>::value, "");
  typedef typename SpaceType::GridLayerType GridLayerType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const size_t dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1>
      ScalarFunctionType;
  typedef Dune::XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::XT::Functions::
      ConstantFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain>
          TensorFunctionType;
  typedef typename XT::LA::Container<RangeFieldType, XT::LA::Backends::istl_sparse>::MatrixType MatrixType;
  typedef typename XT::LA::Container<RangeFieldType, XT::LA::Backends::istl_sparse>::VectorType VectorType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
}; // class OperatorBaseTraits


} // namespace internal


template <class SpaceType>
struct OperatorBase : public ::testing::Test
{
  typedef internal::OperatorBaseTraits<SpaceType> Traits;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef XT::Grid::extract_grid_t<GridLayerType> GridType;
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
    : grid_provider_(XT::Grid::make_cube_grid<GridType>(0.0, 1.0, 6u))
    , space_(grid_provider_.template layer<XT::Grid::Layers::leaf, SpaceType::layer_backend>())
    , scalar_function_("x", "x[0]", 1, "scalar function", {"1.0", "0.0", "0.0"})
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
  using typename BaseType::GridLayerType;

  template <class ProductImp>
  void localizable_product_test(ProductImp& prod)
  {
    const auto& source DUNE_UNUSED = prod.source();
    const auto& range DUNE_UNUSED = prod.range();
    auto& non_const_range DUNE_UNUSED = prod.range();

    XT::Grid::Walker<GridLayerType> walker(this->space_.grid_layer());
    walker.append(prod);
    walker.walk();

    auto result DUNE_UNUSED = prod.apply2();
  } // ... localizable_product_test(...)
}; // class LocalizableProductBase


template <class SpaceType>
struct MatrixOperatorBase : public OperatorBase<SpaceType>
{
  typedef OperatorBase<SpaceType> BaseType;
  using typename BaseType::GridLayerType;
  using typename BaseType::MatrixType;

  template <class OperatorImp>
  void matrix_operator_test(OperatorImp& op)
  {
    const auto& matrix DUNE_UNUSED = op.matrix();
    auto& non_const_matrix DUNE_UNUSED = op.matrix();
    const auto& source_space DUNE_UNUSED = op.source_space();
    const auto& range_space DUNE_UNUSED = op.range_space();

    XT::Grid::Walker<GridLayerType> walker(this->space_.grid_layer());
    walker.append(op);
    walker.walk();
  } // ... matrix_operator_test(...)
}; // class LocalizableProductBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_BASE_HH
