// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_H1_HH
#define DUNE_GDT_TEST_PRODUCTS_H1_HH

#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/h1.hh>

#include "products_elliptic.hh"


template <class SpaceType>
struct H1SemiLocalizableProduct : public EllipticProductBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> BaseType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::TensorType TensorType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    return Products::H1SemiLocalizable<GridViewType, FunctionType, FunctionType>(
               this->space_.grid_view(), function, function)
        .apply2();
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::H1SemiLocalizable<GridViewType, FunctionType, FunctionType> ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_);
    LocalizableProductBase<SpaceType, ProductType>::fulfills_interface(product);
  }
}; // struct H1SemiLocalizableProduct


template <class SpaceType>
struct H1SemiAssemblableProduct : public EllipticAssemblableProduct<SpaceType>
{
  typedef EllipticAssemblableProduct<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<RangeFieldType> MatrixType;
  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection<GridViewType> ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    // create the product
    Products::H1SemiAssemblable<MatrixType, SpaceType, GridViewType, SpaceType> product(this->space_);
    product.assemble(false);
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    const auto result = product.apply2(discrete_function, discrete_function);
    product.assemble(true);
    const auto tbb_result = product.apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(tbb_result, result);
    return result;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::H1SemiAssemblable<MatrixType, SpaceType, GridViewType, SpaceType> ProductType;
    ProductType product(this->space_);
    AssemblableProductBase<SpaceType, ProductType, VectorType>::fulfills_interface(product);
  }
}; // struct H1SemiAssemblableProduct


template <class SpaceType>
struct H1SemiProduct : public EllipticProductBase<SpaceType>
{
  typedef EllipticProductBase<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::TensorType TensorType;
  typedef typename BaseType::GridViewType GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const
  {
    Products::H1Semi<GridViewType> product(this->space_.grid_view());
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::H1Semi<GridViewType> ProductType;
    ProductBase<SpaceType, ProductType>::fulfills_interface(ProductType(this->space_.grid_view()));
  }
}; // struct H1SemiProduct


#endif // DUNE_GDT_TEST_PRODUCTS_H1_HH
