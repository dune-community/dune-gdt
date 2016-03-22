// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH
#define DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/test/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/boundaryl2.hh>

#include "products.hh"

using namespace Dune;
using namespace Dune::GDT;


template< class SpaceType >
struct BoundaryL2ProductBase
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid      GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const size_t                         dimRange = SpaceType::dimRange;
  typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  BoundaryL2ProductBase()
   : grid_(GridProviderType(0.0, 1.0, Stuff::Test::grid_elements()).grid_ptr())
   , space_(Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid_))
   , one_("x", "1.0", 0)
  {}

  virtual RangeFieldType compute(const FunctionType& function) const = 0;

  void constant_arguments() const
  {
    check(compute(one_), 2.0*dimDomain);
  } // ... constant_arguments(...)

  void linear_arguments() const
  {
    const FunctionType linear("x", "x[0] - 1.0", 1);
    check(compute(linear), 1.0 + (2.0*(dimDomain - 1.0))/3.0);
  }

  void quadratic_arguments() const
  {
    const FunctionType quadratic("x", "x[0]*x[0]", 2);
    check(compute(quadratic), 1.0 + (2.0*(dimDomain - 1.0))/5.0);
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon = 1e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon)
        << "result:   " << result << "\n"
        << "expected: " << expected << "\n"
        << "difference: " << std::scientific << error;
  } // ... check(...)

  std::shared_ptr< GridType > grid_;
  const SpaceType space_;
  const FunctionType one_;
}; // struct BoundaryL2ProductBase


template< class SpaceType >
struct BoundaryL2LocalizableProduct
  : public BoundaryL2ProductBase< SpaceType >
{
  typedef BoundaryL2ProductBase< SpaceType > BaseType;
  typedef typename BaseType::GridViewType   GridViewType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    return Products::BoundaryL2Localizable< GridViewType, FunctionType, FunctionType >
        (this->space_.grid_view(), function, function).apply2();
  }

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2Localizable< GridViewType, FunctionType, FunctionType > ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_);
    LocalizableProductBase< SpaceType, ProductType >::fulfills_interface(product);
  }
}; // struct BoundaryL2LocalizableProduct


template< class SpaceType >
struct BoundaryL2AssemblableProduct
  : public BoundaryL2ProductBase< SpaceType >
{
  typedef BoundaryL2ProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< RangeFieldType > MatrixType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    // create the product
    typedef Products::BoundaryL2Assemblable< MatrixType, SpaceType, GridViewType, SpaceType > Product;
    Product product(this->space_);
    product.assemble(false);
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    const auto result = product.apply2(discrete_function, discrete_function);
    Product product_tbb(this->space_);
    product_tbb.assemble(true);
    const auto result_tbb = product_tbb.apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2Assemblable< MatrixType, SpaceType, GridViewType, SpaceType > ProductType;
    ProductType product(this->space_);
    AssemblableProductBase< SpaceType, ProductType, VectorType >::fulfills_interface(product);
  }
}; // struct BoundaryL2AssemblableProduct


template< class SpaceType >
struct BoundaryL2Product
  : public BoundaryL2ProductBase< SpaceType >
{
  typedef BoundaryL2ProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    typedef typename DSG::Intersection< GridViewType >::Type IT;
    const auto view = this->space_.grid_view();
    const Products::BoundaryL2< GridViewType > product(view);
    const auto only_first_half = [=](const GridViewType&, const IT& it){ return it.boundary() && it.geometry().center()[0] < 0.5; };
    const auto only_second_half = [&](const GridViewType&, const IT& it){ return it.boundary() && it.geometry().center()[0] >= 0.5;};
    const Products::BoundaryL2< GridViewType > product_first_half(view, only_first_half);
    const Products::BoundaryL2< GridViewType > product_second_half(view, only_second_half);
    const auto full = product.apply2(function, function);
    const auto first_half = product_first_half.apply2(function, function);
    const auto second_half = product_second_half.apply2(function, function);
    EXPECT_TRUE(DSC::FloatCmp::eq(first_half + second_half, full));
    return full;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2< GridViewType > ProductType;
    ProductBase< SpaceType, ProductType >::fulfills_interface(ProductType(this->space_.grid_view()));
  }
}; // struct BoundaryL2Product


#endif // DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH
