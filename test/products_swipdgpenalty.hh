// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_SWIPDGPENALTY_HH
#define DUNE_GDT_TEST_PRODUCTS_SWIPDGPENALTY_HH

#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/localevaluation/swipdg.hh>
#include <dune/gdt/playground/products/swipdgpenalty.hh>

#include "products.hh"

using namespace Dune;
using namespace Dune::GDT;


template< class SpaceType >
struct SwipdgPenaltyProductBase
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > ScalarType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain >                     TensorType;
  typedef Stuff::Functions::Checkerboard
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > JumpFunctionType;
  typedef Stuff::LocalizableFunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  static DSC::Configuration jump_config()
  {
    DSC::Configuration cfg = JumpFunctionType::default_config();
    cfg["num_elements"] = "[2 1 1]";
    cfg["values"] = "[0 1]";
    return cfg;
  }

  SwipdgPenaltyProductBase()
    : grid_provider_(0.0, 1.0, 2u)
    , space_(grid_provider_.template leaf< SpaceType::part_view_type >())
    , one_("x", "1", 0)
    , unit_matrix_(Stuff::Functions::internal::unit_matrix< RangeFieldType, dimDomain >())
    , continuous_("x", "x[0] - 1", 1)
    , jump_(JumpFunctionType::create(jump_config()))
  {}

  virtual RangeFieldType compute(const FunctionType& function) const = 0;

  void continuous_arguments() const
  {
    // there is no jump in the middle, but at the left boundary we are not zero
    // the 1 is the order of continuous_
    check(compute(continuous_), LocalEvaluation::SIPDG::internal::boundary_sigma(1), 1e-13);
  }

  void discontinuous_arguments() const
  {
    // there is a jump in the middle and at the right boundary we are not zero
    // the 0 is the order of *jump_
    // *0.5 because the value of *jump_ is 0 left and 1 right
    check(compute(*jump_),
          0.5*LocalEvaluation::SIPDG::internal::inner_sigma(0) + LocalEvaluation::SIPDG::internal::boundary_sigma(0));
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon = 1e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon)
        << "result:   " << result << "\n"
        << "expected: " << expected << "\n"
        << "difference: " << std::scientific << error;
  } // ... check(...)

  GridProviderType grid_provider_;
  const SpaceType space_;
  const ScalarType one_;
  const TensorType unit_matrix_;
  const ScalarType continuous_;
  const std::unique_ptr< const JumpFunctionType > jump_;
}; // struct SwipdgPenaltyProductBase


template< class SpaceType >
struct SwipdgPenaltyLocalizableProduct
  : public SwipdgPenaltyProductBase< SpaceType >
{
  typedef SwipdgPenaltyProductBase< SpaceType > BaseType;
  typedef typename BaseType::GridViewType   GridViewType;
  typedef typename BaseType::ScalarType     ScalarType;
  typedef typename BaseType::TensorType     TensorType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    return Products::SwipdgPenaltyLocalizable< GridViewType, ScalarType, TensorType, FunctionType, FunctionType >
        (this->space_.grid_view(), function, function, this->one_, this->unit_matrix_).apply2();
  }

  void fulfills_interface() const
  {
    typedef Products::SwipdgPenaltyLocalizable< GridViewType, ScalarType, TensorType, ScalarType, ScalarType >
      ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_, this->one_, this->unit_matrix_);
    LocalizableProductBase< SpaceType, ProductType >::fulfills_interface(product);
  }
}; // struct SwipdgPenaltyLocalizableProduct


template< class SpaceType >
struct SwipdgPenaltyAssemblableProduct
  : public SwipdgPenaltyProductBase< SpaceType >
{
  typedef SwipdgPenaltyProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::ScalarType      ScalarType;
  typedef typename BaseType::TensorType      TensorType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< RangeFieldType > MatrixType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    // create the product
    typedef Products::SwipdgPenaltyAssemblable< MatrixType, ScalarType, TensorType,
                                                SpaceType, GridViewType, SpaceType > Product;
    Product product(this->space_, this->one_, this->unit_matrix_);
    product.assemble(false);
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    const auto result = product.apply2(discrete_function, discrete_function);
    Product product_tbb(this->space_, this->one_, this->unit_matrix_);
    product_tbb.assemble(true);
    const auto result_tbb = product_tbb.apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::SwipdgPenaltyAssemblable
        < MatrixType, ScalarType, TensorType, SpaceType, GridViewType, SpaceType > ProductType;
    ProductType product(this->space_, this->one_, this->unit_matrix_);
    AssemblableProductBase< SpaceType, ProductType, VectorType >::fulfills_interface(product);
  }
}; // struct SwipdgPenaltyAssemblableProduct


template< class SpaceType >
struct SwipdgPenaltyProduct
  : public SwipdgPenaltyProductBase< SpaceType >
{
  typedef SwipdgPenaltyProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::ScalarType      ScalarType;
  typedef typename BaseType::TensorType      TensorType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    Products::SwipdgPenalty< GridViewType, FunctionType, TensorType >
        product(this->space_.grid_view(), this->one_, this->unit_matrix_);
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::SwipdgPenalty< GridViewType, FunctionType, TensorType > ProductType;
    ProductBase< SpaceType, ProductType >::fulfills_interface(ProductType(this->space_.grid_view(),
                                                                          this->one_, this->unit_matrix_));
  }
}; // struct SwipdgPenaltyProduct


#endif // DUNE_GDT_TEST_PRODUCTS_SWIPDGPENALTY_HH
