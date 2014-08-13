// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products.hh"


template <class SpaceType, class ProductType, class VectorType>
struct AssemblableProduct
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef typename Dune::GDT::AssemblableProductInterface<typename ProductType::Traits> InterfaceType;

  static void fulfills_interface()
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
// static tests
// * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType D_GridViewType;
    typedef typename ProductType::RangeSpaceType D_RangeSpaceType;
    typedef typename ProductType::SourceSpaceType D_SourceSpaceType;
    typedef typename ProductType::MatrixType D_MatrixType;
    typedef typename ProductType::FieldType D_FieldType;
    typedef typename ProductType::PatternType D_PatternType;
    D_PatternType DUNE_UNUSED(d_pattern) = ProductType::pattern(space);
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridViewType I_GridViewType;
    typedef typename InterfaceType::RangeSpaceType I_RangeSpaceType;
    typedef typename InterfaceType::SourceSpaceType I_SourceSpaceType;
    typedef typename InterfaceType::MatrixType I_MatrixType;
    typedef typename InterfaceType::FieldType I_FieldType;
    typedef typename InterfaceType::PatternType I_PatternType;
    static_assert(std::is_same<ProductType, derived_type>::value, "");
    static_assert(std::is_same<D_GridViewType, I_GridViewType>::value, "");
    static_assert(std::is_same<D_RangeSpaceType, I_RangeSpaceType>::value, "");
    static_assert(std::is_same<D_SourceSpaceType, I_SourceSpaceType>::value, "");
    static_assert(std::is_same<D_MatrixType, D_MatrixType>::value, "");
    static_assert(std::is_same<D_FieldType, D_FieldType>::value, "");
    static_assert(std::is_same<D_PatternType, I_PatternType>::value, "");
    I_PatternType DUNE_UNUSED(i_pattern) = InterfaceType::pattern(space);
    // dynamic tests
    D_MatrixType product_matrix(space.mapper().size(), space.mapper().size(), ProductType::pattern(space));
    const VectorType vector(space.mapper().size(), 1.0);
    ProductType product(product_matrix, space, *(space.grid_view()), space);
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view())))
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_RangeSpaceType& d_r = product.range_space();
    if (&d_r != &space)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const D_SourceSpaceType& d_s = product.source_space();
    if (&d_s != &space)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_MatrixType& d_m = product.matrix();
    if (&d_m != &product_matrix)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    D_FieldType d_a = product.apply2(vector, vector);
    // * of the derived type as the interface
    InterfaceType& i_product   = static_cast<InterfaceType&>(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_RangeSpaceType& i_r = i_product.range_space();
    if (&i_r != &d_r)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    const I_SourceSpaceType& i_s = i_product.source_space();
    if (&i_s != &d_s)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_MatrixType& i_m = i_product.matrix();
    if (&i_m != &d_m)
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2(vector, vector);
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a))
      DUNE_THROW_COLORFULLY(Dune::Exception, "");
  }
}; // struct AssemblableProduct


template <class SpaceType>
struct L2AssemblableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<double> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<double> MatrixType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection<GridViewType> ProjectionOperatorType;
  typedef Dune::GDT::Products::L2Assemblable<MatrixType, SpaceType, GridViewType, SpaceType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    VectorType vector(space.mapper().size());
    DiscreteFunctionType discrete_function(space, vector);
    MatrixType product_matrix(space.mapper().size(), space.mapper().size(), ProductType::pattern(space));
    ProductType product(product_matrix, space, *(space.grid_view()), space);
    product.assemble();
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    projection_operator.apply(function_1, discrete_function);
    auto l2_product      = product.apply2(discrete_function.vector(), discrete_function.vector());
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0) << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    projection_operator.apply(function_2, discrete_function);
    l2_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = l2_product - RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 3.0)
                                                   << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    projection_operator.apply(function_3, discrete_function);
    l2_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = l2_product - RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-2)) // <- since the space can be linear only we make quite some projection error
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 5.0)
                                                   << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    AssemblableProduct<SpaceType, ProductType, VectorType>::fulfills_interface();
  }
}; // L2AssemblableProduct


template <class SpaceType>
struct H1SemiAssemblableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<double> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<double> MatrixType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection<GridViewType> ProjectionOperatorType;
  typedef Dune::GDT::Products::H1SemiAssemblable<MatrixType, SpaceType, GridViewType, SpaceType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    VectorType vector(space.mapper().size());
    DiscreteFunctionType discrete_function(space, vector);
    MatrixType product_matrix(space.mapper().size(), space.mapper().size(), ProductType::pattern(space));
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
                                                   << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
    // test 2 (linear)
    const FunctionType function_2(
        "x", "0.5 * x[0] * x[0] - x[0]", 2, "affine gradient", {{"x[0] - 1.0", "0.0", "0.0"}}); // <- this is not as
    // above (b.c. of the
    // projection)
    projection_operator.apply(function_2, discrete_function);
    h1_semi_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = h1_semi_product - RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << RangeFieldType(1.0 / 3.0)
                                                   << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x",
                                  "(1.0/3.0) * x[0] * x[0] * x[0]",
                                  3,
                                  ", quadratic gradient",
                                  {{"x[0]*x[0]", "0.0", "0.0"}}); // <- this is not as above (b.c. of the projection)
    projection_operator.apply(function_3, discrete_function);
    h1_semi_product = product.apply2(discrete_function.vector(), discrete_function.vector());
    error = h1_semi_product - RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << h1_semi_product << " vs. " << RangeFieldType(1.0 / 5.0)
                                                   << " (difference: "
                                                   << std::scientific
                                                   << error
                                                   << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    AssemblableProduct<SpaceType, ProductType, VectorType>::fulfills_interface();
  }
}; // H1SemiAssemblableProduct


TYPED_TEST_CASE(H1SemiAssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiAssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(H1SemiAssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiAssemblableProduct, produces_correct_results)
{
  this->produces_correct_results();
}


TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, produces_correct_results)
{
  this->produces_correct_results();
}

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
