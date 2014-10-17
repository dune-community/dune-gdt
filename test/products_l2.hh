// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_L2_HH
#define DUNE_GDT_TEST_PRODUCTS_L2_HH

#include <dune/gdt/products/l2.hh>

#include "products_weightedl2.hh"


template< class SpaceType >
struct L2LocalizableProduct
  : public WeightedL2LocalizableProduct< SpaceType >
{
  typedef WeightedL2LocalizableProduct< SpaceType > BaseType;
  typedef typename BaseType::FunctionType   FunctionType;
  typedef typename BaseType::GridViewType   GridViewType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    return Products::L2Localizable< GridViewType, FunctionType, FunctionType >(this->space_.grid_view(),
                                                                               function,
                                                                               function).apply2();
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::L2Localizable< GridViewType, FunctionType, FunctionType > ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_);
    LocalizableProductBase< SpaceType, ProductType >::fulfills_interface(product);
  }
}; // struct L2LocalizableProduct


template< class SpaceType >
struct L2AssemblableProduct
  : public WeightedL2ProductBase< SpaceType >
{
  typedef WeightedL2ProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< RangeFieldType > MatrixType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    // create the product
    Products::L2Assemblable< MatrixType, SpaceType, GridViewType, SpaceType > product(this->space_);
    product.assemble();
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    return product.apply2(discrete_function, discrete_function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::L2Assemblable< MatrixType, SpaceType, GridViewType, SpaceType > ProductType;
    ProductType product(this->space_);
    AssemblableProductBase< SpaceType, ProductType, VectorType >::fulfills_interface(product);
  }
}; // struct L2AssemblableProduct



template< class SpaceType >
struct L2Product
  : public WeightedL2ProductBase< SpaceType >
{
  typedef WeightedL2ProductBase< SpaceType > BaseType;
  typedef typename BaseType::RangeFieldType  RangeFieldType;
  typedef typename BaseType::FunctionType    FunctionType;
  typedef typename BaseType::GridViewType    GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const /*DS_OVERIDE DS_FINAL*/
  {
    Products::L2< GridViewType > product(this->space_.grid_view());
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::L2< GridViewType > ProductType;
    ProductBase< SpaceType, ProductType >::fulfills_interface(ProductType(this->space_.grid_view()));
  }
}; // struct L2Product


#endif // DUNE_GDT_TEST_PRODUCTS_L2_HH
