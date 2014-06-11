// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
#   include <dune/grid/sgrid.hh>
#   include <dune/grid/yaspgrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H


#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/continuouslagrange/fem.hh>
#include <dune/gdt/playground/spaces/discontinuouslagrange/fem.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/operators/prolongations.hh>

class errors_are_not_as_expected
  : public Dune::Exception
{};

typedef Dune::Stuff::LA::EigenDenseVector< double > VectorType;

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

// +-------------------------+
// |  * to test the products |
// +-------------------------+

//      - first the interfaces

template< class SpaceType, class ProductType >
struct ProductBase
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
  typedef typename Dune::GDT::ProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface()
  {
    // static tests
    // * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits        Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType  D_GridViewType;
    typedef typename ProductType::FieldType     D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type  derived_type;
    typedef typename InterfaceType::GridViewType  I_GridViewType;
    typedef typename InterfaceType::FieldType     I_FieldType;
    static_assert(std::is_same< ProductType, derived_type >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "1.0", 0, "function", {{"1.0", "1.0", "1.0"}});
    ProductType product(*(space.grid_view()));
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view()))) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_FieldType d_a = product.apply2(function, function);
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2(function, function);
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a)) DUNE_THROW_COLORFULLY(Dune::Exception, "");
  }
}; // struct ProductBase

template< class SpaceType, class ProductType >
struct LocalizableProduct
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
  typedef typename Dune::GDT::LocalizableProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface()
  {
    // static tests
    // * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits        Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType  D_GridViewType;
    typedef typename ProductType::RangeType     D_RangeType;
    typedef typename ProductType::SourceType    D_SourceType;
    typedef typename ProductType::FieldType     D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type  derived_type;
    typedef typename InterfaceType::GridViewType  I_GridViewType;
    typedef typename InterfaceType::RangeType     I_RangeType;
    typedef typename InterfaceType::SourceType    I_SourceType;
    typedef typename InterfaceType::FieldType     I_FieldType;
    static_assert(std::is_same< ProductType, derived_type >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_RangeType, I_RangeType >::value, "");
    static_assert(std::is_same< D_SourceType, I_SourceType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "1.0", 0, "function", {{"1.0", "1.0", "1.0"}});
    ProductType product(*(space.grid_view()), function, function);
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view()))) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_RangeType& d_r = product.range();
    if (&d_r != &function) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_SourceType& d_s = product.source();
    if (&d_s != &function) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_FieldType d_a = product.apply2();
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_RangeType& i_r = i_product.range();
    if (&i_r != &d_r) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_SourceType& i_s = i_product.source();
    if (&i_s != &d_s) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2();
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a)) DUNE_THROW_COLORFULLY(Dune::Exception, "");
  }
}; // struct LocalizableProduct

template< class SpaceType, class ProductType, class VectorType >
struct AssemblableProduct
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef typename Dune::GDT::AssemblableProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface()
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    // static tests
    // * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits          Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType    D_GridViewType;
    typedef typename ProductType::RangeSpaceType  D_RangeSpaceType;
    typedef typename ProductType::SourceSpaceType D_SourceSpaceType;
    typedef typename ProductType::MatrixType      D_MatrixType;
    typedef typename ProductType::FieldType       D_FieldType;
    typedef typename ProductType::PatternType     D_PatternType;
    D_PatternType DUNE_UNUSED(d_pattern) = ProductType::pattern(space);
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type    derived_type;
    typedef typename InterfaceType::GridViewType    I_GridViewType;
    typedef typename InterfaceType::RangeSpaceType  I_RangeSpaceType;
    typedef typename InterfaceType::SourceSpaceType I_SourceSpaceType;
    typedef typename InterfaceType::MatrixType      I_MatrixType;
    typedef typename InterfaceType::FieldType       I_FieldType;
    typedef typename InterfaceType::PatternType     I_PatternType;
    static_assert(std::is_same< ProductType, derived_type >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_RangeSpaceType, I_RangeSpaceType >::value, "");
    static_assert(std::is_same< D_SourceSpaceType, I_SourceSpaceType >::value, "");
    static_assert(std::is_same< D_MatrixType, D_MatrixType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    static_assert(std::is_same< D_PatternType, I_PatternType >::value, "");
    I_PatternType DUNE_UNUSED(i_pattern) = InterfaceType::pattern(space);
    // dynamic tests
    D_MatrixType product_matrix(space.mapper().size(),
                                space.mapper().size(),
                                ProductType::pattern(space));
    const VectorType vector(space.mapper().size(), 1.0);
    ProductType product(product_matrix, space, *(space.grid_view()), space);
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view()))) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_RangeSpaceType& d_r = product.range_space();
    if (&d_r != &space) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_SourceSpaceType& d_s = product.source_space();
    if (&d_s != &space) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_MatrixType& d_m = product.matrix();
    if (&d_m != &product_matrix) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_FieldType d_a = product.apply2(vector, vector);
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_RangeSpaceType& i_r = i_product.range_space();
    if (&i_r != &d_r) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_SourceSpaceType& i_s = i_product.source_space();
    if (&i_s != &d_s) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_MatrixType& i_m = i_product.matrix();
    if (&i_m != &d_m) DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2(vector, vector);
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a)) DUNE_THROW_COLORFULLY(Dune::Exception, "");
  }
}; // struct AssemblableProduct

//      - then the L2 products

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
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const ProductType l2_product_operator(*(space.grid_view()));
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    auto l2_product = l2_product_operator.apply2(function_1, function_1);
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    l2_product = l2_product_operator.apply2(function_2, function_2);
    error = l2_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    l2_product = l2_product_operator.apply2(function_3, function_3);
    error = l2_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  }

  void fulfills_interface() const
  {
    ProductBase< SpaceType, ProductType >::fulfills_interface();
  }
}; // L2ProductOperator

template< class SpaceType >
struct L2LocalizableProduct
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
  typedef Dune::GDT::Products::L2Localizable< GridViewType, FunctionType, FunctionType > ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    ProductType l2_product_operator_1(*(space.grid_view()), function_1, function_1);
    auto l2_product = l2_product_operator_1.apply2();
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    ProductType l2_product_operator_2(*(space.grid_view()), function_2, function_2);
    l2_product = l2_product_operator_2.apply2();
    error = l2_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    ProductType l2_product_operator_3(*(space.grid_view()), function_3, function_3);
    l2_product = l2_product_operator_3.apply2();
    error = l2_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    LocalizableProduct< SpaceType, ProductType >::fulfills_interface();
  }
}; // L2LocalizableProduct

template< class SpaceType >
struct L2AssemblableProduct
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< double > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< double > MatrixType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;
  typedef Dune::GDT::Products::L2Assemblable< MatrixType, SpaceType, GridViewType, SpaceType > ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    VectorType vector(space.mapper().size());
    DiscreteFunctionType discrete_function(space, vector);
    MatrixType product_matrix(space.mapper().size(),
                              space.mapper().size(),
                              ProductType::pattern(space));
    ProductType product(product_matrix, space, *(space.grid_view()), space);
    product.assemble();
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    projection_operator.apply(function_1, discrete_function);
    auto l2_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    projection_operator.apply(function_2, discrete_function);
    l2_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = l2_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    projection_operator.apply(function_3, discrete_function);
    l2_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = l2_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-2)) // <- since the space can be linear only we make quite some projection error
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    AssemblableProduct< SpaceType, ProductType, VectorType >::fulfills_interface();
  }
}; // L2AssemblableProduct

//      - then the H1 semi products

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
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const ProductType h1_semi_product_operator(*(space.grid_view()));
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    auto h1_semi_product = h1_semi_product_operator.apply2(function_1, function_1);
    RangeFieldType error = h1_semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "fake_value", 2, "affine gradient",
                                  {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    h1_semi_product = h1_semi_product_operator.apply2(function_2, function_2);
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "fake_value", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    h1_semi_product = h1_semi_product_operator.apply2(function_3, function_3);
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  }

  void fulfills_interface() const
  {
    ProductBase< SpaceType, ProductType >::fulfills_interface();
  }
}; // H1SemiProductOperator

template< class SpaceType >
struct H1SemiLocalizableProduct
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
  typedef Dune::GDT::Products::H1SemiLocalizable< GridViewType, FunctionType, FunctionType > ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    ProductType h1_semi_product_operator_1(*(space.grid_view()), function_1, function_1);
    auto h1_semi_product = h1_semi_product_operator_1.apply2();
    RangeFieldType error = h1_semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "fake_value", 2, "affine gradient",
                                  {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    ProductType h1_semi_product_operator_2(*(space.grid_view()), function_2, function_2);
    h1_semi_product = h1_semi_product_operator_2.apply2();
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "fake_value", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    ProductType h1_semi_product_operator_3(*(space.grid_view()), function_3, function_3);
    h1_semi_product = h1_semi_product_operator_3.apply2();
    error = h1_semi_product - dimDomain * RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << dimDomain * RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    LocalizableProduct< SpaceType, ProductType >::fulfills_interface();
  }
}; // H1SemiLocalizableProduct

template< class SpaceType >
struct H1SemiAssemblableProduct
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< double > VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix< double > MatrixType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection< GridViewType > ProjectionOperatorType;
  typedef Dune::GDT::Products::H1SemiAssemblable< MatrixType, SpaceType, GridViewType, SpaceType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    VectorType vector(space.mapper().size());
    DiscreteFunctionType discrete_function(space, vector);
    MatrixType product_matrix(space.mapper().size(),
                              space.mapper().size(),
                              ProductType::pattern(space));
    ProductType product(product_matrix, space, *(space.grid_view()), space);
    product.assemble();
    // test 1 (constant)
    const FunctionType function_1("x", "x[0]", 1, "constant gradient", {{"1.0", "0.0", "0.0"}}); // <- this is not the
                                                                                                 //    same as above
                                                                                                 //    (b.c. of the projection)
    projection_operator.apply(function_1, discrete_function);
    auto h1_semi_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    RangeFieldType error = h1_semi_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-14)) // <- this is not as above (b.c. of the projection)
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << RangeFieldType(1.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "0.5 * x[0] * x[0] - x[0]", 2, "affine gradient",
                                  {{"x[0] - 1.0", "0.0", "0.0"}}); // <- this is not as above (b.c. of the projection)
    projection_operator.apply(function_2, discrete_function);
    h1_semi_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = h1_semi_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << RangeFieldType(1.0/3.0)
                            << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "(1.0/3.0) * x[0] * x[0] * x[0]", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "0.0", "0.0"}}); // <- this is not as above (b.c. of the projection)
    projection_operator.apply(function_3, discrete_function);
    h1_semi_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = h1_semi_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << RangeFieldType(1.0/5.0)
                            << " (difference: " << std::scientific << error << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    AssemblableProduct< SpaceType, ProductType, VectorType >::fulfills_interface();
  }
}; // H1SemiAssemblableProduct

// +-------------------------------------+
// |  * to test the projection operators |
// +-------------------------------------+

template< class SpaceType, class ProjectionOperatorType >
struct ProjectionOperatorBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(*(space.grid_view()));
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // ProjectionOperatorType

template< class SpaceType >
struct L2ProjectionOperator
  : public ProjectionOperatorBase< SpaceType, Dune::GDT::Operators::L2Projection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct LagrangeProjectionOperator
  : public ProjectionOperatorBase< SpaceType,
                                   Dune::GDT::Operators::LagrangeProjection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct ProjectionOperator
  : public ProjectionOperatorBase< SpaceType,
                                   Dune::GDT::Operators::Projection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct DirichletProjectionOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 1u); // this has to be 1, otherwise the projection does not equal
    auto grid = grid_provider.grid();             // x[0] any more!
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    DomainType dirichlet_normal(0);
    dirichlet_normal[0] = DomainFieldType(1);
    const Dune::Stuff::Grid::BoundaryInfos::NormalBased< typename GridViewType::Intersection >
        boundary_info(false, {dirichlet_normal});
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::Operators::DirichletProjection< GridViewType > projection_operator(*(space.grid_view()),
                                                                                       boundary_info);
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(*(space.grid_view()));
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // DirichletProjectionOperator

// +---------------------------------------+
// |  * to test the prolongation operators |
// +---------------------------------------+

template< class CoarseSpaceType, class FineSpaceType, class ProlongationOperatorType >
struct ProlongationOperatorBase
{
  typedef typename FineSpaceType::GridViewType      GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                       dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename FineSpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                       dimRange = FineSpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 2u);
    auto grid = grid_provider.grid();
    grid->globalRefine(1);
    const auto coarse_grid_part_view = Dune::GDT::SpaceTools::GridPartView< CoarseSpaceType >::create_level(*grid, 0);
    assert(grid->maxLevel() > 0);
    const auto fine_grid_part_view = Dune::GDT::SpaceTools::GridPartView< FineSpaceType >::create_level(*grid,
                                                                                                        grid->maxLevel());
    assert(fine_grid_part_view->indexSet().size(0) > coarse_grid_part_view->indexSet().size(0));
    // first, project an anlytical function onto the coarse grid
    const FunctionType function("x", "x[0]", 1, "function");
    const CoarseSpaceType coarse_space(coarse_grid_part_view);
    VectorType coarse_vector(coarse_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< CoarseSpaceType, VectorType > CoarseDiscreteFunctionType;
    CoarseDiscreteFunctionType coarse_discrete_function(coarse_space, coarse_vector, "coarse discrete function");
    const Dune::GDT::Operators::Projection< GridViewType > coarse_projection_operator(*(coarse_space.grid_view()));
    coarse_projection_operator.apply(function, coarse_discrete_function);
    // since the projection operator was tested above we are confident this worked
    // but we check anyway (the L2 product operator was also tested above)
    const Dune::GDT::Products::L2< GridViewType > coarse_l2_product_operator(*(coarse_space.grid_view()));
    const Dune::Stuff::Functions::Difference< FunctionType, CoarseDiscreteFunctionType >
        coarse_difference(function, coarse_discrete_function);
    const auto coarse_l2_error = std::sqrt(coarse_l2_product_operator.apply2(coarse_difference, coarse_difference));
    if (coarse_l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(Dune::Stuff::Exceptions::internal_error,
                            "This should not happen, those operators were tested above!\n"
                            << coarse_l2_error << " vs. " << RangeFieldType(1e-15));
    // now we prolong the discrete function from the coarse to the fine grid part
    const FineSpaceType fine_space(fine_grid_part_view);
    VectorType fine_vector(fine_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< FineSpaceType, VectorType > FineDiscreteFunctionType;
    FineDiscreteFunctionType fine_discrete_function(fine_space, fine_vector, "fine discrete function");
    const ProlongationOperatorType prolongation_operator(*(fine_space.grid_view()));
    prolongation_operator.apply(coarse_discrete_function, fine_discrete_function);
    // and measure the error
    const Dune::GDT::Products::L2< GridViewType > fine_l2_product_operator(*(fine_space.grid_view()));
    const Dune::Stuff::Functions::Difference< FunctionType, FineDiscreteFunctionType >
        fine_difference(function, fine_discrete_function);
    const auto fine_l2_error = std::sqrt(fine_l2_product_operator.apply2(fine_difference, fine_difference));
    if (fine_l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected, "\n" << fine_l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // ProlongationOperatorBase

template< class P >
struct L2ProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::L2Prolongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};

template< class P >
struct LagrangeProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::LagrangeProlongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};

template< class P >
struct ProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::Prolongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

// +-----------------------------------------------------------------+
// | Therefore we need to define all grids and gridparts of interest |
// +-----------------------------------------------------------------+
#define SGRID_TYPES(dim) \
  typedef Dune::SGrid< dim, dim >                                 S ## dim ## dGridType; \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S ## dim ## dGridType, false >::Type  S ## dim ## dLeafGridPartType; \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< S ## dim ## dGridType, false >::Type S ## dim ## dLevelGridPartType;
SGRID_TYPES(1)
SGRID_TYPES(2)
SGRID_TYPES(3)
#undef SGRID_TYPES

#define YASPGRID_TYPES(dim) \
  typedef Dune::YaspGrid< dim >                                     Yasp ## dim ## dGridType; \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp ## dim ## dGridType, false >::Type   Yasp ## dim ## dLeafGridPartType; \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< Yasp ## dim ## dGridType, false >::Type  Yasp ## dim ## dLevelGridPartType;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid< 2, 2 > AluConform2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, false >::Type   AluConform2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluConform2dGridType, false >::Type  AluConform2dLevelGridPartType;
typedef Dune::ALUSimplexGrid< 2, 2 > AluSimplex2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex2dGridType, false >::Type   AluSimplex2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluSimplex2dGridType, false >::Type  AluSimplex2dLevelGridPartType;
typedef Dune::ALUSimplexGrid< 3, 3 > AluSimplex3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex3dGridType, false >::Type   AluSimplex3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluSimplex3dGridType, false >::Type  AluSimplex3dLevelGridPartType;
typedef Dune::ALUCubeGrid< 3, 3 > AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluCube3dGridType, false >::Type  AluCube3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluCube3dGridType, false >::Type AluCube3dLevelGridPartType;
#endif

// +----------------------------------------------------+
// |  * arguments for the product operator test structs |
// +----------------------------------------------------+

typedef testing::Types< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLeafGridPartType, 1, double, 1 >

                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLeafGridPartType, 1, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLeafGridPartType, 1, double, 1 >
#endif
                      > ProductOperatorSpaceTypes;

// +-------------------------------------------------------+
// |  * arguments for the projection operator test structs |
// +-------------------------------------------------------+

#define L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID \
    Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLeafGridPartType, 2, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 2, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 2, double, 1 >

#if HAVE_ALUGRID
typedef testing::Types<
                        L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
                      > L2ProjectionOperatorSpaceTypes;
#endif // HAVE_ALUGRID

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES \
    Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLeafGridPartType, 1, double, 1 > \
  \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLeafGridPartType, 1, double, 1 >

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID \
    Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLeafGridPartType, 1, double, 1 > \
  \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 >


typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > LagrangeProjectionOperatorSpaceTypes;

typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
                      , L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > ProjectionOperatorSpaceTypes;

typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > DirichletProjectionOperatorSpaceTypes;

#undef L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID

// +---------------------------------------------------------+
// |  * arguments for the prolongation operator test structs |
// +---------------------------------------------------------+

#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID \
  /* all combinations which have Spaces::DiscontinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
   /* those below do not work in 3d any more! */ \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > >

#if HAVE_ALUGRID
typedef testing::Types<
                        L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
                      > L2ProlongationOperatorSpaceTypes;
#endif

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLevelGridPartType, 1, double, 1 > >

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLevelGridPartType, 1, double, 1 > > \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLevelGridPartType, 1, double, 1 > > \
  \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLevelGridPartType, 1, double, 1 > >

typedef testing::Types<
                        LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > LagrangeProlongationOperatorSpaceTypes;

typedef testing::Types<
                        LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
                      , L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > ProlongationOperatorSpaceTypes;

#undef L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES

// +--------------------------------------------------------------------------------------+
// | 3rd we combine all test structs with their appropriate arguments to create the tests |
// | (comment out the following lines if you do not want a test to be run)                |
// +--------------------------------------------------------------------------------------+

// +--------------------------------------------------+
// | * those have to come first for error measurement |
// +--------------------------------------------------+

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

// +-----------------------------------------------------------------------+
// | * we need the projection operator tests next for reliable projections |
// +-----------------------------------------------------------------------+

#if HAVE_ALUGRID
TYPED_TEST_CASE(L2ProjectionOperator, L2ProjectionOperatorSpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}
#endif // HAVE_ALUGRID

TYPED_TEST_CASE(LagrangeProjectionOperator, LagrangeProjectionOperatorSpaceTypes);
TYPED_TEST(LagrangeProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

TYPED_TEST_CASE(ProjectionOperator, ProjectionOperatorSpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

TYPED_TEST_CASE(DirichletProjectionOperator, DirichletProjectionOperatorSpaceTypes);
TYPED_TEST(DirichletProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

// +----------------------------------------------------------+
// | * now we can continue with the rest of the product tests |
// +----------------------------------------------------------+

TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(H1SemiLocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiLocalizableProduct, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(H1SemiLocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiLocalizableProduct, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(H1SemiAssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiAssemblableProduct, fulfills_interface) {
  this->fulfills_interface();
}

TYPED_TEST_CASE(H1SemiAssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiAssemblableProduct, produces_correct_results) {
  this->produces_correct_results();
}

// +------------------------------------+
// |  * the prolongation operator tests |
// +------------------------------------+

#if HAVE_ALUGRID
TYPED_TEST_CASE(L2ProlongationOperator, L2ProlongationOperatorSpaceTypes);
TYPED_TEST(L2ProlongationOperator, produces_correct_results) {
  this->produces_correct_results();
}
#endif // HAVE_ALUGRID

TYPED_TEST_CASE(LagrangeProlongationOperator, LagrangeProlongationOperatorSpaceTypes);
TYPED_TEST(LagrangeProlongationOperator, produces_correct_results) {
  this->produces_correct_results();
}

TYPED_TEST_CASE(ProlongationOperator, ProlongationOperatorSpaceTypes);
TYPED_TEST(ProlongationOperator, produces_correct_results) {
  this->produces_correct_results();
}

// +--------------------------------------------------------------------------------------+
// | 4th we run all the tests                                                             |
// | (run the resulting executable with '--gtest_catch_exceptions=0' to see an exception) |
// +--------------------------------------------------------------------------------------+

int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
