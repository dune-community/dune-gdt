// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_ELLIPTIC_HH
#define DUNE_GDT_TEST_PRODUCTS_ELLIPTIC_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/elliptic.hh>

#include "products.hh"

using namespace Dune;
using namespace Dune::GDT;


template< class SpaceType >
struct EllipticProductBase
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
  typedef Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >             FunctionType;
  typedef Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain > TensorType;

  EllipticProductBase()
   : grid_(GridProviderType(0.0, 1.0, boost::numeric_cast< size_t >(dsc_grid_elements())).grid_ptr())
   , leaf_view_(Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid_))
   , space_(leaf_view_)
   , one_("x", "1.0", 1, "constant gradient", {{"1.0", "1.0", "1.0"}})
   , unit_matrix_(Stuff::Functions::internal::unit_matrix< RangeFieldType, dimDomain >())
  {}

  virtual RangeFieldType compute(const FunctionType& function) const = 0;

  void constant_arguments() const
  {
    const FunctionType constant_gradient("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    check(compute(constant_gradient), dimDomain*1.0);
  } // ... constant_arguments(...)

  void linear_arguments() const
  {
    const FunctionType linear_gradient("x", "fake_value", 2, "affine gradient",
                                       {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    check(compute(linear_gradient), dimDomain*(1.0/3.0));
  }

  void quadratic_arguments() const
  {
    const FunctionType quadratic_gradient("x", "fake_value", 3, ", quadratic gradient",
                                          {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    check(compute(quadratic_gradient), dimDomain*(1.0/5.0));
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
  typename Dune::GDT::SpaceTools::GridPartView< SpaceType >::LeafGridViewType leaf_view_;
  const SpaceType space_;
  const FunctionType one_;
  const TensorType unit_matrix_;
}; // struct EllipticProductBase


template< class SpaceType >
struct EllipticLocalizableProduct
  : public EllipticProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > BaseType;
  typedef typename BaseType::GridViewType   GridViewType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::TensorType     TensorType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    return Products::EllipticLocalizable
        < GridViewType, FunctionType, FunctionType, FunctionType, RangeFieldType, TensorType >
        (this->space_.grid_view(), function, function, this->one_, this->unit_matrix_).apply2();
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::EllipticLocalizable
        < GridViewType, FunctionType, FunctionType, FunctionType, RangeFieldType, TensorType > ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_, this->one_, this->unit_matrix_);
    LocalizableProductBase< SpaceType, ProductType >::fulfills_interface(product);
  }
}; // struct EllipticLocalizableProduct


template< class SpaceType >
struct EllipticAssemblableProduct
  : public EllipticProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::TensorType      TensorType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< RangeFieldType > MatrixType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    // create the product
    typedef Products::EllipticAssemblable< MatrixType, FunctionType, SpaceType, GridViewType,
                                           SpaceType, RangeFieldType, TensorType > Product;
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

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void constant_arguments() const
  {
    const FunctionType constant_gradient("x", "x[0]", 1, "constant gradient", {{"1.0", "0.0", "0.0"}});
    this->check(compute(constant_gradient), 1.0, 1e-13);
  } // ... constant_arguments(...)

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void linear_arguments() const
  {
    const FunctionType linear_gradient("x", "0.5 * x[0] * x[0] - x[0]", 2, "affine gradient",
                                       {{"x[0] - 1.0", "0.0", "0.0"}});
    this->check(compute(linear_gradient), 1.0/3.0, 1e-13);
  }

  /**
   * \note we can not use the base variant bc. of the projection in compute()
   */
  void quadratic_arguments() const
  {
    const FunctionType quadratic_gradient("x", "(1.0/3.0) * x[0] * x[0] * x[0]", 3, ", quadratic gradient",
                                          {{"x[0]*x[0]", "0.0", "0.0"}});
    this->check(compute(quadratic_gradient), 1.0/5.0);
  }

  void fulfills_interface() const
  {
    typedef Products::EllipticAssemblable
        < MatrixType, FunctionType, SpaceType, GridViewType, SpaceType, RangeFieldType, TensorType > ProductType;
    ProductType product(this->space_, this->one_, this->unit_matrix_);
    AssemblableProductBase< SpaceType, ProductType, VectorType >::fulfills_interface(product);
  }
}; // EllipticAssemblableProduct


template< class SpaceType >
struct EllipticProduct
  : public EllipticProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::TensorType      TensorType;
  typedef typename BaseType::GridViewType    GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    Products::Elliptic< GridViewType, FunctionType, RangeFieldType, TensorType >
        product(this->space_.grid_view(), this->one_, this->unit_matrix_);
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::Elliptic< GridViewType, FunctionType, RangeFieldType, TensorType > ProductType;
    ProductBase< SpaceType, ProductType >::fulfills_interface(ProductType(this->space_.grid_view(),
                                                                          this->one_,
                                                                          this->unit_matrix_));
  }
}; // struct EllipticProduct


template< class SpaceType >
struct SimplifiedEllipticLocalizableProduct
  : public EllipticProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType > BaseType;
  typedef typename BaseType::GridViewType   GridViewType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::TensorType     TensorType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    return Products::EllipticLocalizable< GridViewType, FunctionType, FunctionType, FunctionType >
        (this->space_.grid_view(), function, function, this->one_).apply2();
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::EllipticLocalizable< GridViewType, FunctionType, FunctionType, FunctionType > ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_, this->one_);
    LocalizableProductBase< SpaceType, ProductType >::fulfills_interface(product);
  }
}; // struct SimplifiedEllipticLocalizableProduct


template< class SpaceType >
struct SimplifiedEllipticAssemblableProduct
  : public EllipticAssemblableProduct< SpaceType >
{
  typedef EllipticAssemblableProduct< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< RangeFieldType > MatrixType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    // create the product
    typedef Products::EllipticAssemblable< MatrixType, FunctionType, SpaceType, GridViewType, SpaceType > Product;
    Product product(this->space_, this->one_);
    product.assemble(false);
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    const auto result = product.apply2(discrete_function, discrete_function);
    Product product_tbb(this->space_, this->one_);
    product_tbb.assemble(true);
    const auto result_tbb = product_tbb.apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::EllipticAssemblable< MatrixType, FunctionType, SpaceType, GridViewType, SpaceType > ProductType;
    ProductType product(this->space_, this->one_);
    AssemblableProductBase< SpaceType, ProductType, VectorType >::fulfills_interface(product);
  }
}; // struct SimplifiedEllipticAssemblableProduct


template< class SpaceType >
struct SimplifiedEllipticProduct
  : public EllipticProductBase< SpaceType >
{
  typedef EllipticProductBase< SpaceType >  BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::TensorType     TensorType;
  typedef typename BaseType::GridViewType   GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    Products::Elliptic< GridViewType, FunctionType > product(this->space_.grid_view(), this->one_);
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::Elliptic< GridViewType, FunctionType > ProductType;
    ProductBase< SpaceType, ProductType >::fulfills_interface(ProductType(this->space_.grid_view(), this->one_));
  }
}; // struct SimplifiedEllipticProduct


#endif // DUNE_GDT_TEST_PRODUCTS_ELLIPTIC_HH
