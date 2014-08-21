// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products.hh"

template< class SpaceType >
struct L2ProductOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef Dune::GDT::Products::L2< GridViewType > ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const ProductType l2_product_operator(*(space.grid_view()));
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    auto l2_product = l2_product_operator.apply2(function_1, function_1);
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    l2_product = l2_product_operator.apply2(function_2, function_2);
    error = l2_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    l2_product = l2_product_operator.apply2(function_3, function_3);
    error = l2_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  }

  void fulfills_interface() const
  {
    ProductBase< SpaceType, ProductType >::fulfills_interface();
  }
}; // L2ProductOperator


template< class SpaceType >
struct H1SemiProductOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef Dune::GDT::Products::H1SemiGeneric< GridViewType > ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const ProductType h1_semi_product_operator(*(space.grid_view()));
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    auto h1_semi_product = h1_semi_product_operator.apply2(function_1, function_1);
    RangeFieldType error = h1_semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "fake_value", 2, "affine gradient",
                                  {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    h1_semi_product = h1_semi_product_operator.apply2(function_2, function_2);
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "fake_value", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    h1_semi_product = h1_semi_product_operator.apply2(function_3, function_3);
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  }

  void fulfills_interface() const
  {
    ProductBase< SpaceType, ProductType >::fulfills_interface();
  }
}; // H1SemiProductOperator


TYPED_TEST_CASE(L2ProductOperator, ProductOperatorSpaceTypes);
TYPED_TEST(L2ProductOperator, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2ProductOperator, ProductOperatorSpaceTypes);
TYPED_TEST(L2ProductOperator, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(H1SemiProductOperator, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiProductOperator, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(H1SemiProductOperator, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiProductOperator, produces_correct_results) {
  this->produces_correct_results();
}


#include <dune/stuff/test/test_main.cxx>
