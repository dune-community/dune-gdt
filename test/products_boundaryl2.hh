// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH
#define DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH

#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/boundaryl2.hh>

#include "products.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class SpaceType>
struct BoundaryL2ProductBase : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange> FunctionType;

  BoundaryL2ProductBase()
    : grid_(GridProviderType(0.0, 1.0, 3u).grid_ptr())
    , space_(Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_))
    , one_("x", "1.0", 0)
  {
  }

  virtual RangeFieldType compute(const FunctionType& function) const = 0;

  void constant_arguments() const
  {
    check(compute(one_), 2.0 * dimDomain);
  } // ... constant_arguments(...)

  void linear_arguments() const
  {
    const FunctionType linear("x", "x[0] - 1.0", 1);
    check(compute(linear), 1.0 + (2.0 * (dimDomain - 1.0)) / 3.0);
  }

  void quadratic_arguments() const
  {
    const FunctionType quadratic("x", "x[0]*x[0]", 2);
    check(compute(quadratic), 1.0 + (2.0 * (dimDomain - 1.0)) / 5.0);
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon = 1e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon) << "result:   " << result << "\n"
                              << "expected: " << expected << "\n"
                              << "difference: " << std::scientific << error;
  } // ... check(...)

  std::shared_ptr<GridType> grid_;
  const SpaceType space_;
  const FunctionType one_;
}; // struct BoundaryL2ProductBase


template <class SpaceType>
struct BoundaryL2LocalizableProduct : public BoundaryL2ProductBase<SpaceType>
{
  typedef BoundaryL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    return Products::BoundaryL2Localizable<GridViewType, FunctionType, FunctionType>(
               this->space_.grid_view(), function, function)
        .apply2();
  }

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2Localizable<GridViewType, FunctionType, FunctionType> ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_);
    LocalizableProductBase<SpaceType, ProductType>::fulfills_interface(product);
  }
}; // struct BoundaryL2LocalizableProduct


template <class SpaceType>
struct BoundaryL2AssemblableProduct : public BoundaryL2ProductBase<SpaceType>
{
  typedef BoundaryL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<RangeFieldType> MatrixType;
  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection<GridViewType> ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    // create the product
    Products::BoundaryL2Assemblable<MatrixType, SpaceType, GridViewType, SpaceType> product(this->space_);
    product.assemble();
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    return product.apply2(discrete_function, discrete_function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2Assemblable<MatrixType, SpaceType, GridViewType, SpaceType> ProductType;
    ProductType product(this->space_);
    AssemblableProductBase<SpaceType, ProductType, VectorType>::fulfills_interface(product);
  }
}; // struct BoundaryL2AssemblableProduct


template <class SpaceType>
struct BoundaryL2Product : public BoundaryL2ProductBase<SpaceType>
{
  typedef BoundaryL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::GridViewType GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    const Products::BoundaryL2<GridViewType> product(this->space_.grid_view());
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::BoundaryL2<GridViewType> ProductType;
    ProductBase<SpaceType, ProductType>::fulfills_interface(ProductType(this->space_.grid_view()));
  }
}; // struct BoundaryL2Product


#endif // DUNE_GDT_TEST_PRODUCTS_BOUNDARYL2_HH
