// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products.hh"


template <class SpaceType, class ProductType>
struct LocalizableProduct
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef typename Dune::GDT::LocalizableProductInterface<typename ProductType::Traits> InterfaceType;

  static void fulfills_interface()
  {
// static tests
// * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType D_GridViewType;
    typedef typename ProductType::RangeType D_RangeType;
    typedef typename ProductType::SourceType D_SourceType;
    typedef typename ProductType::FieldType D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridViewType I_GridViewType;
    typedef typename InterfaceType::RangeType I_RangeType;
    typedef typename InterfaceType::SourceType I_SourceType;
    typedef typename InterfaceType::FieldType I_FieldType;
    static_assert(std::is_same<ProductType, derived_type>::value, "");
    static_assert(std::is_same<D_GridViewType, I_GridViewType>::value, "");
    static_assert(std::is_same<D_RangeType, I_RangeType>::value, "");
    static_assert(std::is_same<D_SourceType, I_SourceType>::value, "");
    static_assert(std::is_same<D_FieldType, D_FieldType>::value, "");
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "1.0", 0, "function", {{"1.0", "1.0", "1.0"}});
    ProductType product(*(space.grid_view()), function, function);
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view())))
      DUNE_THROW(Dune::Exception, "");
    const D_RangeType& d_r = product.range();
    if (&d_r != &function)
      DUNE_THROW(Dune::Exception, "");
    const D_SourceType& d_s = product.source();
    if (&d_s != &function)
      DUNE_THROW(Dune::Exception, "");
    D_FieldType d_a = product.apply2();
    // * of the derived type as the interface
    InterfaceType& i_product   = static_cast<InterfaceType&>(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp)
      DUNE_THROW(Dune::Exception, "");
    const I_RangeType& i_r = i_product.range();
    if (&i_r != &d_r)
      DUNE_THROW(Dune::Exception, "");
    const I_SourceType& i_s = i_product.source();
    if (&i_s != &d_s)
      DUNE_THROW(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2();
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a))
      DUNE_THROW(Dune::Exception, "");
  }
}; // struct LocalizableProduct


template <class SpaceType>
struct L2LocalizableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::GDT::Products::L2Localizable<GridViewType, FunctionType, FunctionType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    ProductType l2_product_operator_1(*(space.grid_view()), function_1, function_1);
    auto l2_product      = l2_product_operator_1.apply2();
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    ProductType l2_product_operator_2(*(space.grid_view()), function_2, function_2);
    l2_product = l2_product_operator_2.apply2();
    error = l2_product - RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 3.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    ProductType l2_product_operator_3(*(space.grid_view()), function_3, function_3);
    l2_product = l2_product_operator_3.apply2();
    error = l2_product - RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 5.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    LocalizableProduct<SpaceType, ProductType>::fulfills_interface();
  }
}; // L2LocalizableProduct


template <class SpaceType>
struct H1SemiLocalizableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::GDT::Products::H1SemiLocalizable<GridViewType, FunctionType, FunctionType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    ProductType h1_semi_product_operator_1(*(space.grid_view()), function_1, function_1);
    auto h1_semi_product = h1_semi_product_operator_1.apply2();
    RangeFieldType error = h1_semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 2 (linear)
    const FunctionType function_2(
        "x", "fake_value", 2, "affine gradient", {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    ProductType h1_semi_product_operator_2(*(space.grid_view()), function_2, function_2);
    h1_semi_product = h1_semi_product_operator_2.apply2();
    error = h1_semi_product - dimDomain * RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0 / 3.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 3 (quadratic)
    const FunctionType function_3(
        "x", "fake_value", 3, ", quadratic gradient", {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    ProductType h1_semi_product_operator_3(*(space.grid_view()), function_3, function_3);
    h1_semi_product = h1_semi_product_operator_3.apply2();
    error = h1_semi_product - dimDomain * RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0 / 5.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    LocalizableProduct<SpaceType, ProductType>::fulfills_interface();
  }
}; // H1SemiLocalizableProduct


TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, produces_correct_results)
{
  this->produces_correct_results();
}


TYPED_TEST_CASE(H1SemiLocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiLocalizableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(H1SemiLocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiLocalizableProduct, produces_correct_results)
{
  this->produces_correct_results();
}


#include <dune/stuff/test/test_main.cxx>
