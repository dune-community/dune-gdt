// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_HH
#define DUNE_GDT_TEST_PRODUCTS_HH

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/products/interfaces.hh>
#include <dune/gdt/spaces/tools.hh>


template< class SpaceType, class ProductType >
struct ProductBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const size_t                         dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef typename Dune::GDT::ProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface(const ProductType& product)
  {
    // static tests
    // * of the derived type
    typedef typename ProductType::Traits        Traits;
    typedef typename ProductType::GridViewType  D_GridViewType;
    typedef typename ProductType::FieldType     D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type  derived_type;
    typedef typename InterfaceType::GridViewType  I_GridViewType;
    typedef typename InterfaceType::FieldType     I_FieldType;
    static_assert(std::is_base_of< derived_type, ProductType >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    const FunctionType function("x", "1.0", 0, "function", {{"1.0", "1.0", "1.0"}}); // <- it is on purpose that the
    // dynamic tests                                                                 //    gradient (1) is not the real
    // * of the derived type                                                         //    gradient of the function
    const D_GridViewType& d_gv = product.grid_view();
    D_FieldType d_a = product.apply2(function, function);
    // * of the derived type as the interface
    const InterfaceType& i_product = static_cast< const InterfaceType& >(product);
    const I_GridViewType& i_gv = i_product.grid_view();
    EXPECT_EQ(&i_gv, &d_gv);
    I_FieldType i_a = i_product.apply2(function, function);
    EXPECT_EQ(i_a, d_a);
  }
}; // struct ProductBase


template< class SpaceType, class ProductType >
struct LocalizableProductBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const size_t                         dimRange = SpaceType::dimRange;
  typedef typename Dune::GDT::LocalizableProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface(ProductType& product)
  {
    // static tests
    // * of the derived type
    typedef typename ProductType::Traits        Traits;
    typedef typename ProductType::GridViewType  D_GridViewType;
    typedef typename ProductType::RangeType     D_RangeType;
    typedef typename ProductType::SourceType    D_SourceType;
    typedef typename ProductType::FieldType     D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridViewType I_GridViewType;
    typedef typename InterfaceType::RangeType    I_RangeType;
    typedef typename InterfaceType::SourceType   I_SourceType;
    typedef typename InterfaceType::FieldType    I_FieldType;
    static_assert(std::is_base_of< derived_type, ProductType >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_RangeType, I_RangeType >::value, "");
    static_assert(std::is_same< D_SourceType, I_SourceType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gv = product.grid_view();
    const D_RangeType& d_r = product.range();
    const D_SourceType& d_s = product.source();
    D_FieldType d_a = product.apply2();
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    EXPECT_EQ(&i_gp, &d_gv);
    const I_RangeType& i_r = i_product.range();
    EXPECT_EQ(&i_r, &d_r);
    const I_SourceType& i_s = i_product.source();
    EXPECT_EQ(&i_s, &d_s);
    I_FieldType i_a = i_product.apply2();
    EXPECT_EQ(i_a, d_a);
  }
}; // struct LocalizableProductBase


template< class SpaceType, class ProductType, class VectorType >
struct AssemblableProductBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const size_t                         dimRange = SpaceType::dimRange;
  typedef typename Dune::GDT::AssemblableProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface(ProductType& product)
  {
    // static tests
    // * of the derived type
    typedef typename ProductType::Traits          Traits;
    typedef typename ProductType::GridViewType    D_GridViewType;
    typedef typename ProductType::RangeSpaceType  D_RangeSpaceType;
    typedef typename ProductType::SourceSpaceType D_SourceSpaceType;
    typedef typename ProductType::MatrixType      D_MatrixType;
    typedef typename ProductType::FieldType       D_FieldType;
    typedef typename ProductType::PatternType     D_PatternType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type    derived_type;
    typedef typename InterfaceType::GridViewType    I_GridViewType;
    typedef typename InterfaceType::RangeSpaceType  I_RangeSpaceType;
    typedef typename InterfaceType::SourceSpaceType I_SourceSpaceType;
    typedef typename InterfaceType::MatrixType      I_MatrixType;
    typedef typename InterfaceType::FieldType       I_FieldType;
    typedef typename InterfaceType::PatternType     I_PatternType;
    static_assert(std::is_base_of< derived_type, ProductType >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_RangeSpaceType, I_RangeSpaceType >::value, "");
    static_assert(std::is_same< D_SourceSpaceType, I_SourceSpaceType >::value, "");
    static_assert(std::is_same< D_MatrixType, D_MatrixType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    static_assert(std::is_same< D_PatternType, I_PatternType >::value, "");
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gv = product.grid_view();
    const D_RangeSpaceType& d_r = product.range_space();
    const D_SourceSpaceType& d_s = product.source_space();
    D_PatternType d_p = ProductType::pattern(d_r, d_s, d_gv);
    VectorType vector(d_r.mapper().size(), 1);
    D_MatrixType& d_m = product.matrix();
    D_FieldType d_a = product.apply2(vector, vector);
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    I_PatternType i_p = product.pattern(d_r, d_s, d_gv);
    EXPECT_EQ(i_p, d_p);
    const I_GridViewType& i_gv = i_product.grid_view();
    EXPECT_EQ(&i_gv, &d_gv);
    const I_RangeSpaceType& i_r = i_product.range_space();
    EXPECT_EQ(&i_r, &d_r);
    const I_SourceSpaceType& i_s = i_product.source_space();
    EXPECT_EQ(&i_s, &d_s);
    I_MatrixType& i_m = i_product.matrix();
    EXPECT_EQ(&i_m, &d_m);
    I_FieldType i_a = i_product.apply2(vector, vector);
    EXPECT_EQ(i_a, d_a);
  }
}; // struct AssemblableProductBase


#endif // DUNE_GDT_TEST_PRODUCTS_HH
